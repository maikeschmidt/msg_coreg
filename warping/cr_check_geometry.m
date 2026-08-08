function R = cr_check_geometry(geom, opts)
% cr_check_geometry - Pre-flight validation of a geometry before FEM meshing
%
% TetGen (through iso2mesh surf2mesh) is far stricter than the BEM about
% surface quality. A geometry that produces a perfectly good BEM lead field
% can fail to tetrahedralise, and the failure surfaces deep inside surf2mesh
% with no indication of which mesh is at fault. This runs the checks TetGen
% cares about, cheaply, BEFORE committing hours to a solve.
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
%     .max_probe     vertices sampled per mesh for nesting (default 2000)
%     .verbose       print a report (default true)
%
% OUTPUT:
%   R - struct with:
%     .ok           true if nothing fatal was found
%     .fatal        cell array of problems that will stop TetGen
%     .warnings     cell array of problems worth looking at
%     .mesh         per-mesh struct array: name, n_vert, n_face, min_angle,
%                   n_slivers, n_degenerate, closed, n_open_edges, euler,
%                   volume, selfint
%     .gaps         pairwise minimum vertex-to-vertex distances
%
% NOTES
%   - Pure MATLAB apart from the optional self-intersection test, which uses
%     iso2mesh meshcheckrepair if it is on the path. No toolboxes required.
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
if ~isfield(opts,'max_probe'),      opts.max_probe      = 2000;  end
if ~isfield(opts,'verbose'),        opts.verbose        = true;  end
if ~isfield(opts,'nested')
    opts.nested = setdiff(opts.meshes, {opts.outer}, 'stable');
end

R.fatal    = {};
R.warnings = {};
R.mesh     = struct('name',{},'n_vert',{},'n_face',{},'min_angle',{}, ...
                    'n_slivers',{},'n_degenerate',{},'closed',{}, ...
                    'n_open_edges',{},'euler',{},'volume',{},'selfint',{});

present = opts.meshes(cellfun(@(m) isfield(geom, m), opts.meshes));
if isempty(present)
    R.ok = false;
    R.fatal{end+1} = 'No mesh_* fields found in this geometry.';
    if opts.verbose, fprintf('%s\n', R.fatal{end}); end
    return;
end


% CHECK 1: PER-MESH TRIANGLE QUALITY AND TOPOLOGY

for m = 1:numel(present)
    name = present{m};
    [V, F] = get_mesh(geom.(name));

    S = struct('name', name, 'n_vert', size(V,1), 'n_face', size(F,1));

    % Triangle angles. Slivers — very small minimum angle — are what makes
    % TetGen give up, and anisotropic scaling manufactures them.
    [ang_min, area] = tri_quality(V, F);
    S.min_angle    = min(ang_min);
    S.n_slivers    = sum(ang_min < opts.min_angle_deg);
    S.n_degenerate = sum(area < eps(max(area)) * 1e3);

    % Edge manifoldness: every edge of a closed surface is in exactly 2 faces
    E  = sort([F(:,[1 2]); F(:,[2 3]); F(:,[3 1])], 2);
    [~, ~, ic] = unique(E, 'rows');
    cnt = accumarray(ic, 1);
    S.n_open_edges = sum(cnt ~= 2);
    S.closed       = S.n_open_edges == 0;
    S.euler        = size(V,1) - numel(cnt) + size(F,1);
    S.volume       = mesh_volume(V, F);

    S.selfint = NaN;
    if opts.check_selfint && exist('meshcheckrepair','file') == 2
        try
            evalc('[~, fi] = meshcheckrepair(V, F, ''intersect'');');
            S.selfint = size(fi,1) < size(F,1);   % faces removed => had them
        catch
            S.selfint = NaN;
        end
    end

    % Verdicts
    if S.n_degenerate > 0
        R.fatal{end+1} = sprintf('%s: %d zero-area faces', name, S.n_degenerate);
    end
    if S.n_slivers > 0
        R.fatal{end+1} = sprintf('%s: %d faces with an angle below %.2f deg (min %.4f)', ...
            name, S.n_slivers, opts.min_angle_deg, S.min_angle);
    end
    if ~S.closed
        R.fatal{end+1} = sprintf('%s: surface not closed, %d non-manifold edges', ...
            name, S.n_open_edges);
    end
    if S.euler ~= 2 && S.closed
        R.warnings{end+1} = sprintf('%s: Euler characteristic %d (a sphere is 2) — genus %g', ...
            name, S.euler, (2 - S.euler)/2);
    end
    if isequal(S.selfint, true)
        R.fatal{end+1} = sprintf('%s: surface self-intersects', name);
    end
    if S.volume <= 0
        R.warnings{end+1} = sprintf('%s: signed volume %.4g — face winding may be inverted', ...
            name, S.volume);
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
            R.fatal{end+1} = sprintf('%s / %s: surfaces %.4g apart (limit %.4g)', ...
                present{i}, present{j}, dmin, opts.min_gap);
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
            R.fatal{end+1} = sprintf(...
                '%s: %d of %d sampled vertices lie OUTSIDE %s', ...
                nm, n_out, numel(in), opts.outer);
        end
    end
end


% CHECK 4: SOURCES AND SENSORS

if isfield(geom, 'sources_cent') && isfield(geom, 'mesh_wm')
    [Vw, Fw] = get_mesh(geom.mesh_wm);
    in = points_inside(geom.sources_cent, Vw, Fw);
    if any(~in)
        R.fatal{end+1} = sprintf('sources_cent: %d of %d sources outside mesh_wm', ...
            sum(~in), numel(in));
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
        [Vo, Fo] = get_mesh(geom.(opts.outer));
        in = points_inside(sp, Vo, Fo);
        if any(in)
            R.fatal{end+1} = sprintf(...
                '%s: %d of %d sensors are INSIDE %s', ...
                sens_field, sum(in), numel(in), opts.outer);
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

function in = points_inside(P, V, F)
% Point-in-closed-surface by ray casting along +X, counting crossings.
% Uses tt_is_inside when it is on the path, since that is what the FEM
% pipeline itself uses for seed placement — same answer, same edge cases.
    if exist('tt_is_inside','file') == 2
        in = false(size(P,1),1);
        for i = 1:size(P,1)
            in(i) = logical(tt_is_inside(P(i,:), V, F));
        end
        return;
    end

    A = V(F(:,1),:); B = V(F(:,2),:); C = V(F(:,3),:);
    in = false(size(P,1),1);
    for i = 1:size(P,1)
        p = P(i,:);
        % Triangles whose Y/Z box contains the ray
        cand = min([A(:,2) B(:,2) C(:,2)],[],2) <= p(2) & ...
               max([A(:,2) B(:,2) C(:,2)],[],2) >= p(2) & ...
               min([A(:,3) B(:,3) C(:,3)],[],2) <= p(3) & ...
               max([A(:,3) B(:,3) C(:,3)],[],2) >= p(3) & ...
               max([A(:,1) B(:,1) C(:,1)],[],2) >= p(1);
        if ~any(cand), continue; end

        a = A(cand,:); b = B(cand,:); c = C(cand,:);
        e1 = b - a; e2 = c - a;
        dir = [1 0 0];
        h = [zeros(size(e2,1),1), e2(:,3), -e2(:,2)];    % cross(dir, e2)
        det_ = sum(e1 .* h, 2);
        ok = abs(det_) > 1e-14;
        s  = p - a;
        u  = sum(s .* h, 2) ./ det_;
        q  = cross(s, e1, 2);
        vv = q(:,1) ./ det_;                              % dot(dir, q)
        t  = sum(e2 .* q, 2) ./ det_;
        hit = ok & u >= 0 & u <= 1 & vv >= 0 & (u + vv) <= 1 & t > 1e-12;
        in(i) = mod(sum(hit), 2) == 1;
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
        fprintf('\nFAIL — %d problem(s):\n', numel(R.fatal));
        for k = 1:numel(R.fatal), fprintf('  [FATAL] %s\n', R.fatal{k}); end
    end
    for k = 1:numel(R.warnings), fprintf('  [warn ] %s\n', R.warnings{k}); end

    if ~isempty(R.fatal)
        fprintf(['\nMost slivers after an affine warp come from the warp\n' ...
                 'amplifying triangles that were already elongated. Try\n' ...
                 'cr_repair_geometry, which remeshes each surface at a fixed\n' ...
                 'edge length and re-runs these checks.\n']);
    end
    fprintf('\n');
end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end
