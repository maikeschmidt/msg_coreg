function I = cr_find_intersections(geom, opts)
% cr_find_intersections - Name the meshes that intersect, using TetGen
%
% TetGen reports self-intersections like this:
%
%   PLC Error:  A segment and a facet intersect at point (0.0403,-0.106,-0.0373).
%     Segment: [6849,6840] #-1 (0)
%     Facet:   [16102,15509,14024] #0
%
% which tells you a collision exists but not WHICH anatomical structures are
% involved — the indices are into a merged mesh that no longer knows the
% difference between a vertebra and a lung. This runs the same detection,
% then maps every index back to the mesh it came from, so the output reads
% "mesh_bone component 7 intersects mesh_wm at (x,y,z)".
%
% WHY TETGEN AND NOT AN OWN INTERSECTION TEST
%   TetGen is the tool that will accept or reject the mesh. Anything else
%   risks disagreeing with it — passing a check and still failing to mesh,
%   or vice versa. Using -d asks the actual judge.
%
% USAGE:
%   I = cr_find_intersections('geometries_coreg_rep02_inhomo.mat');
%   I = cr_find_intersections(geom_struct, opts);
%
% INPUT:
%   geom - geometry struct or path to a geometry .mat
%   opts - optional:
%     .meshes    meshes to include (default: all mesh_* fields)
%     .split     split each mesh into connected components before mapping,
%                so a 23-vertebra bone mesh reports WHICH vertebra
%                (default true)
%     .max_report  maximum collisions to list (default 40)
%     .keep_poly   keep the temporary .poly for inspection (default false)
%     .verbose     default true
%
% OUTPUT:
%   I - struct array: mesh_a, comp_a, mesh_b, comp_b, point
%       plus I(1).n_total if more were found than reported
%
% REQUIRES: iso2mesh on the path (for mcpath/mwpath and the TetGen binary).
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
if ~isfield(opts,'split'),      opts.split      = true;  end
if ~isfield(opts,'max_report'), opts.max_report = 40;    end
if ~isfield(opts,'keep_poly'),  opts.keep_poly  = false; end
if ~isfield(opts,'verbose'),    opts.verbose    = true;  end

if ischar(geom) || isstring(geom)
    src  = char(geom);
    geom = load(src);
else
    src = '<struct>';
end

if ~isfield(opts,'meshes')
    fn = fieldnames(geom);
    opts.meshes = fn(startsWith(fn,'mesh_') & ...
        cellfun(@(f) isstruct(geom.(f)), fn))';
end
present = opts.meshes(cellfun(@(m) isfield(geom, m), opts.meshes));

if exist('mwpath','file') ~= 2
    error('iso2mesh is not on the path — mwpath/mcpath are required.');
end


% MERGE, KEEPING PROVENANCE

V = []; F = [];
owner_name = {};      % per-vertex source mesh
owner_comp = [];      % per-vertex component index within that mesh

for m = 1:numel(present)
    [Vm, Fm] = get_mesh(geom.(present{m}));

    if opts.split
        comp = connected_components(Vm, Fm);
    else
        comp = ones(size(Vm,1), 1);
    end

    off = size(V,1);
    V   = [V; Vm];                 %#ok<AGROW>
    F   = [F; Fm + off];           %#ok<AGROW>
    owner_name = [owner_name; repmat(present(m), size(Vm,1), 1)];   %#ok<AGROW>
    owner_comp = [owner_comp; comp];                                %#ok<AGROW>
end

if opts.verbose
    fprintf('=== Intersection search: %s ===\n', src);
    fprintf('Merged %d meshes -> %d vertices, %d faces\n', ...
        numel(present), size(V,1), size(F,1));
end


% WRITE A MINIMAL .POLY AND RUN TETGEN -d
% Written directly rather than through savesurfpoly: no regions, no holes,
% nothing that could itself be malformed. The only question being asked is
% whether the surfaces collide.

poly = mwpath('cr_intersect_check.poly');
fid  = fopen(poly, 'w');
fprintf(fid, '#node list\n%d 3 0 0\n', size(V,1));
fprintf(fid, '%d %.17g %.17g %.17g\n', [(1:size(V,1))', V]');
% Facet section. The header's second number declares whether boundary
% markers are present; it is 0 here, so each facet is exactly
%   <#polygons> <#holes>
%   <#corners> <corner indices>
% Writing a third number on the first line — a marker the header said was
% absent — desynchronises TetGen's parser and it rejects the file, which is
% then indistinguishable from "no intersections found" unless the run is
% validated. See the status check below.
fprintf(fid, '#facet list\n%d 0\n', size(F,1));
fprintf(fid, '1 0\n3 %d %d %d\n', F');
fprintf(fid, '#hole list\n0\n#region list\n0\n');
fclose(fid);

exe = '';
for cand = {'tetgen1.5','tetgen'}
    try
        e = mcpath(cand{1}, getexeext);
        if isfile(e) || isfile([e getexeext]), exe = e; break; end
    catch
    end
end
if isempty(exe)
    error('Could not locate a TetGen binary via mcpath.');
end

[status, out] = system(sprintf(' "%s" -d "%s"', exe, poly));

if ~opts.keep_poly
    delete(poly);
elseif opts.verbose
    fprintf('Kept %s\n', poly);
end

% VALIDATE THAT TETGEN ACTUALLY RAN
%
% A failed invocation produces no parseable collisions, which is exactly
% what a clean geometry produces. Reporting "ok" in that case is worse than
% useless — it is a false all-clear on a geometry that will fail the FEM.
% So a clean verdict requires POSITIVE evidence that TetGen inspected the
% file, not merely the absence of complaints.
found_int = contains(out, 'PLC Error') || ...
            contains(lower(out), 'self-intersection');
ran_ok    = contains(out, 'Delaunizing') || contains(out, 'Opening');

if ~found_int && ~ran_ok
    error('cr_find_intersections:tetgenFailed', ...
        ['TetGen did not produce interpretable output (exit status %d), ' ...
         'so no conclusion can be drawn about intersections.\n' ...
         'Command: "%s" -d "%s"\nOutput:\n%s'], ...
        status, exe, poly, out);
end

R_raw = out;


% PARSE AND MAP BACK

I = struct('mesh_a',{},'comp_a',{},'mesh_b',{},'comp_b',{},'point',{});

% Both wordings TetGen uses
pts  = regexp(out, 'intersect at point \(([^)]*)\)', 'tokens');
segs = regexp(out, 'Segment:\s*\[(\d+),(\d+)\]', 'tokens');
facs = regexp(out, 'Facet:\s*\[(\d+),(\d+),(\d+)\]', 'tokens');
fac2 = regexp(out, 'Facet 1:\s*\[(\d+),(\d+),(\d+)\]', 'tokens');
if isempty(facs) && ~isempty(fac2), facs = fac2; end

n = max([numel(segs), numel(facs)]);
for k = 1:min(n, opts.max_report)
    a_idx = [];
    if k <= numel(segs), a_idx = str2double(segs{k}); end
    b_idx = [];
    if k <= numel(facs), b_idx = str2double(facs{k}); end
    if isempty(a_idx) || isempty(b_idx), continue; end

    [na, ca] = describe(a_idx, owner_name, owner_comp);
    [nb, cb] = describe(b_idx, owner_name, owner_comp);

    if k <= numel(pts)
        pt = sscanf(strrep(pts{k}{1}, ',', ' '), '%f')';
    else
        pt = [NaN NaN NaN];
    end

    I(end+1) = struct('mesh_a', na, 'comp_a', ca, ...
                      'mesh_b', nb, 'comp_b', cb, 'point', pt); %#ok<AGROW>
end

if opts.verbose
    if isempty(I)
        if found_int
            fprintf(['\nTetGen reported a self-intersection but none could ' ...
                     'be parsed. Raw output:\n%s\n'], out);
        else
            fprintf('\nTetGen inspected the file and found no self-intersections.\n');
        end
    else
        fprintf('\n%d collision(s) reported', n);
        if n > opts.max_report
            fprintf(' (showing the first %d)', opts.max_report);
        end
        fprintf(':\n\n');
        for k = 1:numel(I)
            fprintf('  %-22s  x  %-22s  at (%.4g, %.4g, %.4g)\n', ...
                label(I(k).mesh_a, I(k).comp_a), ...
                label(I(k).mesh_b, I(k).comp_b), I(k).point);
        end

        % The pair that matters is usually obvious from a tally
        pairs = arrayfun(@(x) sprintf('%s x %s', ...
            label(x.mesh_a, x.comp_a), label(x.mesh_b, x.comp_b)), ...
            I, 'UniformOutput', false);
        [u, ~, ix] = unique(pairs);
        cnt = accumarray(ix, 1);
        [cnt, o] = sort(cnt, 'descend');
        fprintf('\nBy pair:\n');
        for k = 1:numel(u)
            fprintf('  %4d  %s\n', cnt(k), u{o(k)});
        end
        fprintf(['\nA rigid or affine transform cannot create a collision ' ...
                 'that\nwas not already there, so if these meshes overlap ' ...
                 'here they\noverlap in the base geometry too. Check the ' ...
                 'unmodified model\nbefore blaming the coregistration or ' ...
                 'the warp.\n']);
    end
    fprintf('\n');
end

if ~isempty(I)
    I(1).n_total = n;
    I(1).raw     = R_raw;
end

end


% LOCAL FUNCTIONS

function [V, F] = get_mesh(M)
    if isfield(M,'vertices'), V = M.vertices; else, V = M.v; end
    if isfield(M,'faces'),    F = M.faces;    else, F = M.f; end
    V = double(V); F = double(F(:,1:3));
end

function c = connected_components(V, F)
% Component index per vertex, after merging coincident vertices so that
% surfaces stored with duplicated vertices are not split spuriously.
    n = size(V,1);
    tol = 1e-6 * norm(max(V,[],1) - min(V,[],1));
    if tol <= 0, tol = 1e-9; end
    [~, ~, vmap] = unique(round(V / tol), 'rows');

    nm = max(vmap);
    Fm = vmap(F);
    A  = sparse([Fm(:,1); Fm(:,2); Fm(:,3)], ...
                [Fm(:,2); Fm(:,3); Fm(:,1)], 1, nm, nm);
    A  = A + A';
    [~, cm] = graph_components(A);
    c = cm(vmap);
    c = c(:);
    if numel(c) ~= n, c = ones(n,1); end
end

function [n_comp, lab] = graph_components(A)
% Breadth-first labelling. dmperm would be shorter but is not available in
% every install; this needs nothing beyond sparse.
    n   = size(A,1);
    lab = zeros(n,1);
    n_comp = 0;
    for s = 1:n
        if lab(s), continue; end
        n_comp = n_comp + 1;
        stack  = s;
        lab(s) = n_comp;
        while ~isempty(stack)
            v = stack(end); stack(end) = [];
            nb = find(A(:,v));
            nb = nb(lab(nb) == 0);
            lab(nb) = n_comp;
            stack = [stack; nb]; %#ok<AGROW>
        end
    end
end

function [name, comp] = describe(idx, owner_name, owner_comp)
% TetGen indices are 1-based into the merged node list.
    idx = idx(idx >= 1 & idx <= numel(owner_name));
    if isempty(idx), name = '?'; comp = NaN; return; end
    nm = owner_name(idx);
    name = nm{1};
    cc = owner_comp(idx);
    comp = mode(cc);
    if ~all(strcmp(nm, name)), name = [name ' (+others)']; end
end

function s = label(name, comp)
    if isnan(comp) || comp <= 1
        s = name;
    else
        s = sprintf('%s #%d', name, comp);
    end
end
