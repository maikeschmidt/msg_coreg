function R = cr_check_geometry(geom, opts)
% cr_check_geometry - Pre-flight validation of a geometry before FEM meshing
%
% TetGen (through iso2mesh surf2mesh) is far stricter than the BEM about
% surface quality. A geometry that produces a perfectly good BEM lead field
% can fail to tetrahedralise, and the failure surfaces deep inside surf2mesh
% with no indication of which mesh is at fault. This runs the checks TetGen
% cares about, cheaply, BEFORE committing hours to a solve.
%
% CALIBRATE AGAINST A GEOMETRY THAT WORKS
%   Absolute thresholds for "too thin a triangle" are guesswork. TetGen
%   accepts far worse than intuition suggests — the coregistration
%   geometries in this study contain faces down to 0.12 degrees and mesh
%   without complaint, because TetGen inserts Steiner points rather than
%   refusing. So pass opts.ref_mesh from a geometry you have ALREADY
%   tetrahedralised successfully, and read the output as "worse than the
%   one that worked", not as an absolute verdict.
%
% WHAT ACTUALLY BREAKS AFTER AN AFFINE WARP
%   cr_build_warp_geometries applies ONE affine transform to every mesh, so
%   nesting is preserved by construction: an affine map cannot move the cord
%   outside the bone if it started inside. What it does do is DEGRADE
%   TRIANGLE QUALITY — anisotropic scaling turns an already-elongated
%   triangle into a sliver, and slivers are what TetGen rejects. So check 1
%   below is the prime suspect, not the nesting checks.
%
%   Conversely, warping ONLY the torso and leaving the internal meshes
%   unwarped DOES break nesting, and would make this worse rather than
%   better. That is what checks 4 and 5 exist to catch.
%
% USAGE:
%   R = cr_check_geometry('geometries_warp01_realistic.mat');
%   R = cr_check_geometry(geom_struct, opts);
%
% INPUT:
%   geom - geometry struct, or path to a geometry .mat
%   opts - optional struct:
%     .meshes        cell of mesh field names to check
%                    (default: every field named mesh_* that is present)
%     .outer         name of the enclosing mesh (default 'mesh_torso')
%     .nested        meshes that must lie inside .outer
%                    (default: every checked mesh except .outer)
%     .min_angle_deg smallest acceptable triangle angle (default 1.0)
%     .min_gap       smallest acceptable gap between surfaces, in the
%                    geometry's own units (default [] = report only)
%     .check_selfint run iso2mesh self-intersection test (default true)
%     .check_nesting run point-in-surface nesting test (default true)
%     .max_probe     vertices sampled per mesh for nesting (default 200)
%     .ref_mesh      per-mesh stats from a geometry KNOWN to tetrahedralise
%                    (the .mesh field of an earlier result). Quality is then
%                    judged relative to it instead of against absolute
%                    thresholds. Strongly recommended — see below.
%     .ref_angle_tol flag if min angle is below ref * this (default 0.5)
%     .ref_sliver_tol flag if sliver count exceeds ref * this (default 3)
%     .verbose       print a report (default true)
%
% OUTPUT:
%   R - struct with:
%     .ok           true if nothing fatal was found
%     .fatal        problems that will genuinely stop TetGen
%     .warnings     things worth looking at, including quality that is much
%                   worse than the reference geometry
%     .quality      sliver counts per mesh — informational, never fatal
%     .mesh         per-mesh struct array: name, n_vert, n_face, min_angle,
%                   n_slivers, n_degenerate, closed, n_open_edges, euler,
%                   volume, selfint
%     .gaps         pairwise minimum vertex-to-vertex distances
%
% NOTES
%   - Pure MATLAB apart from the optional self-intersection test, which uses
%     iso2mesh meshcheckrepair if it is on the path. No toolboxes required.
%   - The point-in-surface test is a vectorised solid-angle (winding number)
%     computation, not a call to tt_is_inside — same maths, but tt_is_inside
%     handles one point per call, which is far too slow here.
%   - Distances are in whatever units the geometry uses; the detected scale
%     is printed so a mm/m mix-up is obvious.
%
% SEE ALSO: cr_check_warp_geometries (batch over a whole warp set)
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of msg_coreg.

if nargin < 2, opts = struct(); end

if ischar(geom) || isstring(geom)
    R.file = char(geom);
    geom   = load(char(geom));
else
    R.file = '<struct>';
end

if ~isfield(opts,'meshes')
    fn = fieldnames(geom);
    opts.meshes = fn(startsWith(fn, 'mesh_') & ...
        cellfun(@(f) isstruct(geom.(f)), fn))';
end
if ~isfield(opts,'outer'),          opts.outer          = 'mesh_torso'; end
if ~isfield(opts,'min_angle_deg'),  opts.min_angle_deg  = 1.0;   end
if ~isfield(opts,'min_gap'),        opts.min_gap        = [];    end
if ~isfield(opts,'check_selfint'),  opts.check_selfint  = true;  end
if ~isfield(opts,'check_nesting'),  opts.check_nesting  = true;  end
if ~isfield(opts,'max_probe'),      opts.max_probe      = 200;   end
if ~isfield(opts,'ref_mesh'),       opts.ref_mesh       = [];    end
if ~isfield(opts,'ref_angle_tol'),  opts.ref_angle_tol  = 0.5;   end
if ~isfield(opts,'ref_sliver_tol'), opts.ref_sliver_tol = 3;     end
if ~isfield(opts,'verbose'),        opts.verbose        = true;  end
if ~isfield(opts,'nested')
    opts.nested = setdiff(opts.meshes, {opts.outer}, 'stable');
end

R.fatal     = {};
R.fatal_cat = {};
R.warnings  = {};
R.quality   = {};

% SELF-CALIBRATION
% opts.tolerated lists problem CATEGORIES that a geometry known to
% tetrahedralise also exhibits. Anything in that list cannot be what stops
% TetGen — it demonstrably does not — so it is reported as a warning rather
% than a blocker. This is what stops an invented threshold from condemning
% geometries that mesh perfectly well.
if ~isfield(opts,'tolerated'), opts.tolerated = {}; end
R.mesh     = struct('name',{},'n_vert',{},'n_face',{},'min_angle',{}, ...
                    'n_slivers',{},'n_degenerate',{},'closed',{}, ...
                    'n_open_edges',{},'euler',{},'volume',{},'selfint',{});

present = opts.meshes(cellfun(@(m) isfield(geom, m), opts.meshes));
if isempty(present)
    R.ok = false;
    R.gaps = struct('a',{},'b',{},'min_dist',{});
    R.fatal{end+1}     = 'No mesh_* fields found in this geometry.';
    R.fatal_cat{end+1} = 'nomesh';
    if opts.verbose, fprintf('%s\n', R.fatal{end}); end
    return;
end


% CHECK 1: PER-MESH TRIANGLE QUALITY AND TOPOLOGY

for m = 1:numel(present)
    name = present{m};
    [V, F] = get_mesh(geom.(name));

    % Non-finite coordinates. TetGen cannot parse these and exits WITHOUT
    % printing anything, which is the empty 'Tetgen command failed:' message.
    % Always fatal — no working geometry can contain one, so this is never
    % demoted by the reference calibration.
    n_bad = sum(any(~isfinite(V), 2));
    if n_bad > 0
        R.fatal{end+1}     = sprintf('%s: %d vertices are NaN or Inf', name, n_bad);
        R.fatal_cat{end+1} = 'nonfinite';
    end
    if any(F(:) < 1) || any(F(:) > size(V,1)) || any(~isfinite(F(:)))
        R.fatal{end+1}     = sprintf('%s: face indices out of range', name);
        R.fatal_cat{end+1} = 'badindex';
    end

    % Triangle angles. Slivers — very small minimum angle — are what makes
    % TetGen give up, and anisotropic scaling manufactures them.
    [ang_min, area] = tri_quality(V, F);
    min_angle    = min(ang_min);
    n_slivers    = sum(ang_min < opts.min_angle_deg);
    n_degenerate = sum(area < eps(max(area)) * 1e3);

    % Edge manifoldness: every edge of a closed surface is in exactly 2 faces.
    % Vertices are merged by position first — meshes imported from STL carry
    % duplicated vertices, and on those the raw index-based edge test reports
    % every edge as open even though the surface is perfectly closed.
    tol  = 1e-6 * norm(max(V,[],1) - min(V,[],1));
    if tol <= 0, tol = 1e-9; end
    [~, ~, vmap] = unique(round(V / tol), 'rows');
    Fm = vmap(F);
    E  = sort([Fm(:,[1 2]); Fm(:,[2 3]); Fm(:,[3 1])], 2);
    E  = E(E(:,1) ~= E(:,2), :);          % drop edges collapsed by the merge
    [~, ~, ic] = unique(E, 'rows');
    cnt = accumarray(ic, 1);
    n_open_edges = sum(cnt ~= 2);
    closed       = n_open_edges == 0;
    euler        = size(V,1) - numel(cnt) + size(F,1);
    volume       = mesh_volume(V, F);

    selfint = NaN;
    if opts.check_selfint && exist('meshcheckrepair','file') == 2
        try
            [~, fi] = evalc('mcr_intersect(V, F)');
            selfint = size(fi,1) < size(F,1);   % faces removed => had them
        catch
            selfint = NaN;
        end
    end

    % Built in ONE call, in the same field order R.mesh was declared with.
    % Assigning a struct whose fields are in a different order into a struct
    % array is an error in MATLAB, not a silent reorder.
    S = struct('name', name, 'n_vert', size(V,1), 'n_face', size(F,1), ...
               'min_angle', min_angle, 'n_slivers', n_slivers, ...
               'n_degenerate', n_degenerate, 'closed', closed, ...
               'n_open_edges', n_open_edges, 'euler', euler, ...
               'volume', volume, 'selfint', selfint);

    % Verdicts
    %
    % FATAL means "TetGen will refuse this". Sliver triangles are NOT in that
    % category: the coregistration geometries carry faces down to 0.12 deg
    % and tetrahedralise perfectly well, because TetGen inserts Steiner
    % points rather than giving up. Treating slivers as fatal flagged every
    % known-good geometry in the set, so they are reported as QUALITY and
    % only escalate when they are much worse than a reference geometry that
    % is known to mesh (opts.ref_mesh).
    if n_degenerate > 0
        R = add_problem(R, opts, 'degenerate', ...
            sprintf('%s: %d zero-area faces', name, n_degenerate));
    end
    if n_slivers > 0
        R.quality{end+1} = sprintf('%s: %d faces below %.2f deg (min %.4f)', ...
            name, n_slivers, opts.min_angle_deg, min_angle);
    end
    if ~isempty(opts.ref_mesh)
        ri = find(strcmp({opts.ref_mesh.name}, name), 1);
        if ~isempty(ri)
            rs = opts.ref_mesh(ri);
            if min_angle < rs.min_angle * opts.ref_angle_tol
                R.warnings{end+1} = sprintf(...
                    '%s: min angle %.4f deg is %.1fx worse than the reference (%.4f)', ...
                    name, min_angle, rs.min_angle / max(min_angle, eps), rs.min_angle);
            end
            if n_slivers > max(rs.n_slivers * opts.ref_sliver_tol, rs.n_slivers + 5)
                R.warnings{end+1} = sprintf(...
                    '%s: %d slivers vs %d in the reference', ...
                    name, n_slivers, rs.n_slivers);
            end
        end
    end
    if ~closed
        R = add_problem(R, opts, 'open', ...
            sprintf('%s: surface not closed, %d non-manifold edges', ...
            name, n_open_edges));
    end
    if euler ~= 2 && closed
        R.warnings{end+1} = sprintf('%s: Euler characteristic %d (a sphere is 2) — genus %g', ...
            name, euler, (2 - euler)/2);
    end
    if isequal(selfint, true)
        R = add_problem(R, opts, 'selfint', ...
            sprintf('%s: surface self-intersects', name));
    end
    if volume <= 0
        R.warnings{end+1} = sprintf('%s: signed volume %.4g — face winding may be inverted', ...
            name, volume);
    end

    R.mesh(end+1) = S; %#ok<AGROW>
end


% CHECK 2: PAIRWISE MINIMUM DISTANCE BETWEEN SURFACES
% Surfaces that touch, or nearly touch, leave TetGen no room to place a
% tetrahedron between them.

R.gaps = struct('a',{},'b',{},'min_dist',{});
for i = 1:numel(present)
    [Vi, ~] = get_mesh(geom.(present{i}));
    for j = i+1:numel(present)
        [Vj, ~] = get_mesh(geom.(present{j}));
        dmin = min_pair_distance(Vi, Vj);
        R.gaps(end+1) = struct('a', present{i}, 'b', present{j}, ...
            'min_dist', dmin); %#ok<AGROW>
        if ~isempty(opts.min_gap) && dmin < opts.min_gap
            R = add_problem(R, opts, 'gap', ...
                sprintf('%s / %s: surfaces %.4g apart (limit %.4g)', ...
                present{i}, present{j}, dmin, opts.min_gap));
        end
    end
end


% CHECK 3: NESTING
% Affine warping preserves this, so a failure means either the base geometry
% was already wrong, or the meshes were NOT all warped by the same transform
% — which is exactly what happens if you warp the torso alone.

if opts.check_nesting && isfield(geom, opts.outer)
    [Vo, Fo] = get_mesh(geom.(opts.outer));
    for k = 1:numel(opts.nested)
        nm = opts.nested{k};
        if ~isfield(geom, nm), continue; end
        [Vn, ~] = get_mesh(geom.(nm));
        P = subsample(Vn, opts.max_probe);
        in = points_inside(P, Vo, Fo);
        n_out = sum(~in);
        if n_out > 0
            R = add_problem(R, opts, 'nesting', sprintf(...
                '%s: %d of %d sampled vertices lie OUTSIDE %s', ...
                nm, n_out, numel(in), opts.outer));
        end
    end
end


% CHECK 4: SOURCES AND SENSORS

if isfield(geom, 'sources_cent') && isfield(geom, 'mesh_wm')
    [Vw, Fw] = get_mesh(geom.mesh_wm);
    SP = coords(geom.sources_cent);
    if any(~isfinite(SP(:)))
        R.fatal{end+1}     = sprintf('sources_cent: %d non-finite coordinates', ...
            sum(any(~isfinite(SP), 2)));
        R.fatal_cat{end+1} = 'nonfinite';
    end
    in = points_inside(SP, Vw, Fw);
    if any(~in)
        R = add_problem(R, opts, 'sources', ...
            sprintf('sources_cent: %d of %d sources outside mesh_wm', ...
            sum(~in), numel(in)));
    end
end

sens_field = '';
for f = {'sens','grad','sensors'}
    if isfield(geom, f{1}), sens_field = f{1}; break; end
end
if ~isempty(sens_field) && isfield(geom, opts.outer)
    sp = [];
    for f = {'coilpos','chanpos','pos'}
        if isfield(geom.(sens_field), f{1})
            sp = geom.(sens_field).(f{1}); break;
        end
    end
    if ~isempty(sp)
        if any(~isfinite(sp(:)))
            R.fatal{end+1}     = sprintf('%s: %d non-finite sensor coordinates', ...
                sens_field, sum(any(~isfinite(sp), 2)));
            R.fatal_cat{end+1} = 'nonfinite';
        end
        % Orientations are renormalised under the inverse-transpose when a
        % geometry is warped, so a degenerate direction divides by zero and
        % silently produces NaN. Checked because it is a known failure mode.
        for oflds = {'coilori','chanori','ori'}
            if isfield(geom.(sens_field), oflds{1})
                O = geom.(sens_field).(oflds{1});
                if isnumeric(O) && any(~isfinite(O(:)))
                    R.fatal{end+1}     = sprintf('%s.%s: %d non-finite orientations', ...
                        sens_field, oflds{1}, sum(any(~isfinite(O), 2)));
                    R.fatal_cat{end+1} = 'nonfinite';
                end
            end
        end
        [Vo, Fo] = get_mesh(geom.(opts.outer));
        in = points_inside(sp, Vo, Fo);
        if any(in)
            R = add_problem(R, opts, 'sensors', sprintf(...
                '%s: %d of %d sensors are INSIDE %s', ...
                sens_field, sum(in), numel(in), opts.outer));
        end
    end
end


R.ok = isempty(R.fatal);

if opts.verbose, print_report(R, opts); end

end


% LOCAL FUNCTIONS

function [V, F] = get_mesh(M)
    if isfield(M,'vertices'), V = M.vertices; else, V = M.v; end
    if isfield(M,'faces'),    F = M.faces;    else, F = M.f; end
    F = double(F(:,1:3));
    V = double(V);
end

function [ang_min, area] = tri_quality(V, F)
% Minimum interior angle (degrees) and area of every triangle.
    A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
    ab = B - A; bc = C - B; ca = A - C;
    la = sqrt(sum(bc.^2, 2));   % side opposite A
    lb = sqrt(sum(ca.^2, 2));
    lc = sqrt(sum(ab.^2, 2));

    area = 0.5 * sqrt(sum(cross(ab, -ca, 2).^2, 2));

    % Law of cosines, clamped against round-off on near-degenerate triangles
    ca_ = clamp((lb.^2 + lc.^2 - la.^2) ./ (2 .* lb .* lc));
    cb_ = clamp((la.^2 + lc.^2 - lb.^2) ./ (2 .* la .* lc));
    cc_ = clamp((la.^2 + lb.^2 - lc.^2) ./ (2 .* la .* lb));

    ang = real(acosd([ca_, cb_, cc_]));
    ang(~isfinite(ang)) = 0;      % zero-length sides => degenerate
    ang_min = min(ang, [], 2);
end

function x = clamp(x)
    x(x >  1) =  1;
    x(x < -1) = -1;
    x(~isfinite(x)) = 1;
end

function v = mesh_volume(V, F)
% Signed volume by the divergence theorem — positive for outward normals.
    A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
    v = sum(dot(A, cross(B, C, 2), 2)) / 6;
end

function d = min_pair_distance(A, B)
% Smallest vertex-to-vertex distance, chunked so a 10k x 10k pair does not
% allocate a 800 MB matrix. Vertex-to-vertex slightly overestimates the true
% surface gap, which is the safe direction for a warning.
    d = inf;
    chunk = max(1, floor(2e6 / max(1, size(B,1))));
    for i = 1:chunk:size(A,1)
        blk = A(i:min(i+chunk-1, size(A,1)), :);
        D2  = sum(blk.^2, 2) + sum(B.^2, 2)' - 2 * (blk * B');
        d   = min(d, sqrt(max(0, min(D2(:)))));
    end
end

function P = subsample(V, n)
    if size(V,1) <= n, P = V; return; end
    P = V(round(linspace(1, size(V,1), n)), :);
end

function R = add_problem(R, opts, category, msg)
% A problem is only FATAL if a geometry known to tetrahedralise does not
% also have it. Otherwise it is demonstrably survivable, and saying
% otherwise would condemn working geometries.
    if any(strcmp(opts.tolerated, category))
        R.warnings{end+1} = sprintf(...
            '%s  [reference has this too — not what blocks TetGen]', msg);
    else
        R.fatal{end+1}     = msg;
        R.fatal_cat{end+1} = category;
    end
end

function P = coords(X)
% Pull an N x 3 coordinate array out of whatever the geometry stores.
% sources_cent is a STRUCT with a .pos field, not a bare array, and passing
% the struct straight through is what made this crash inside tt_is_inside.
    if isnumeric(X)
        P = double(X);
    elseif isstruct(X)
        for f = {'pos','vertices','coilpos','chanpos','p','v'}
            if isfield(X, f{1}) && isnumeric(X.(f{1}))
                P = double(X.(f{1}));
                return;
            end
        end
        error('cr_check_geometry:coords', ...
            'Cannot find coordinates in struct with fields: %s', ...
            strjoin(fieldnames(X)', ', '));
    else
        error('cr_check_geometry:coords', ...
            'Expected numeric coordinates, got %s.', class(X));
    end
    if size(P,2) ~= 3 && size(P,1) == 3, P = P'; end
end

function in = points_inside(P, V, F)
% Point-in-closed-surface by solid angle (winding number): the surface
% subtends 4*pi from an interior point and 0 from an exterior one, so the
% test is |sum| > 2*pi. Robust to concavity, and it does not care about
% face winding because the magnitude is taken.
%
% Vectorised over points AND faces in chunks. tt_is_inside does the same
% maths but one point per call, which is far too slow to probe thousands of
% vertices across several meshes — that loop is what had to be interrupted.
    P = coords(P);
    V = double(V);
    F = double(F(:,1:3));

    nF = size(F,1);
    nP = size(P,1);
    in = false(nP,1);
    if nF == 0 || nP == 0, return; end

    A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
    Ar = reshape(A, 1, nF, 3);
    Br = reshape(B, 1, nF, 3);
    Cr = reshape(C, 1, nF, 3);

    % Keep each chunk's temporaries small. The intermediates are
    % nPoints x nFaces x 3 and there are about ten of them live at once, so
    % a generous budget here thrashes memory rather than going faster.
    chunk = max(1, floor(2e5 / nF));

    for i0 = 1:chunk:nP
        idx = i0:min(i0 + chunk - 1, nP);
        p   = reshape(P(idx,:), numel(idx), 1, 3);

        a = Ar - p;  b = Br - p;  c = Cr - p;
        la = sqrt(sum(a.^2, 3));
        lb = sqrt(sum(b.^2, 3));
        lc = sqrt(sum(c.^2, 3));

        num = sum(a .* cross(b, c, 3), 3);
        den = la .* lb .* lc ...
            + sum(a .* b, 3) .* lc ...
            + sum(b .* c, 3) .* la ...
            + sum(c .* a, 3) .* lb;

        omega = 2 * atan2(num, den);
        in(idx) = abs(sum(omega, 2)) > 2 * pi;
    end
end

function print_report(R, opts)
    fprintf('\n=== Geometry check: %s ===\n', R.file);
    fprintf('%-18s %8s %8s %10s %8s %8s %7s\n', ...
        'mesh', 'verts', 'faces', 'min angle', 'slivers', 'open e', 'closed');
    for m = 1:numel(R.mesh)
        S = R.mesh(m);
        fprintf('%-18s %8d %8d %10.4f %8d %8d %7s\n', ...
            S.name, S.n_vert, S.n_face, S.min_angle, S.n_slivers, ...
            S.n_open_edges, yesno(S.closed));
    end

    if ~isempty(R.gaps)
        fprintf('\nClosest approach between surfaces (geometry units):\n');
        [~, ord] = sort([R.gaps.min_dist]);
        for k = ord(1:min(5, numel(ord)))
            fprintf('  %-16s %-16s %10.4g\n', ...
                R.gaps(k).a, R.gaps(k).b, R.gaps(k).min_dist);
        end
    end

    if isempty(R.fatal)
        fprintf('\nPASS — nothing here should stop TetGen.\n');
    else
        fprintf('\nFAIL — %d blocking problem(s):\n', numel(R.fatal));
        for k = 1:numel(R.fatal), fprintf('  [FATAL] %s\n', R.fatal{k}); end
    end
    for k = 1:numel(R.warnings), fprintf('  [warn ] %s\n', R.warnings{k}); end
    for k = 1:numel(R.quality),  fprintf('  [qual ] %s\n', R.quality{k});  end

    if ~isempty(R.quality) && isempty(opts.ref_mesh)
        fprintf(['\nSliver counts above are informational. TetGen tolerates\n' ...
                 'very thin triangles, so they only matter in comparison:\n' ...
                 'pass opts.ref_mesh from a geometry you have already meshed\n' ...
                 'successfully to see whether these are unusually bad.\n']);
    end
    if ~isempty(R.fatal)
        fprintf(['\nIf the blocking problems are quality-related, try\n' ...
                 'cr_repair_geometry, which remeshes each surface and\n' ...
                 're-runs these checks.\n']);
    end
    fprintf('\n');
end

function [V, F] = mcr_intersect(V, F)
% Wrapper so evalc can suppress meshfix's console chatter with a plain call.
    [V, F] = meshcheckrepair(V, F, 'intersect');
end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end
