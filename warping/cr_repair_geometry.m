function [geom, R] = cr_repair_geometry(geom, opts)
% cr_repair_geometry - Make a warped geometry tetrahedralisable again
%
% An affine warp cannot un-nest the anatomy, but it DOES stretch triangles.
% A triangle that was merely elongated in the base mesh becomes a sliver
% after a 0.85/1.15 anisotropic scale, and slivers are what makes TetGen
% give up. The fix is not to warp less anatomy — it is to restore surface
% quality after warping.
%
% WHAT IT DOES, in order, stopping as soon as the geometry passes:
%   1. remove duplicate and unused vertices        (meshcheckrepair 'dup')
%   2. fix orientation and open edges              (meshcheckrepair 'deep')
%   3. remove self-intersections                   (meshcheckrepair 'meshfix')
%   4. isotropic remesh at a fixed edge length     (remeshsurf / meshresample)
%
% Step 4 is the one that actually cures slivers, and it is also the most
% invasive: it changes the vertex count, so a repaired geometry is no longer
% row-matched to the unrepaired one. That does not affect the BEM-vs-FEM
% comparison — both solvers read the same repaired geometry — but it does
% mean you should repair BEFORE computing either solver's lead field, never
% between them.
%
% USAGE:
%   [geom, R] = cr_repair_geometry(geom);
%   [geom, R] = cr_repair_geometry('geometries_warp07_realistic.mat', opts);
%
% INPUT:
%   geom - geometry struct or path to a geometry .mat
%   opts - optional struct:
%     .meshes        meshes to repair (default: all mesh_* fields)
%     .max_stage     highest repair stage to attempt, 1-4 (default 4)
%     .keepratio     vertex fraction for remeshing (default 0.95). NOT 1.0 —
%                    meshresample at 1.0 is close to a no-op, so stage 4
%                    would not actually rebuild the triangulation. Lower it
%                    further (0.8) if slivers survive.
%     .check_opts    options forwarded to cr_check_geometry
%     .outfile       if set, save the repaired geometry here
%     .verbose       default true
%
% OUTPUT:
%   geom - repaired geometry struct
%   R    - cr_check_geometry result AFTER repair
%
% REQUIRES: iso2mesh on the path (meshcheckrepair, meshresample).
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
if ~isfield(opts,'max_stage'),  opts.max_stage  = 4;         end
if ~isfield(opts,'keepratio'),  opts.keepratio  = 0.95;      end
if ~isfield(opts,'check_opts'), opts.check_opts = struct();  end
if ~isfield(opts,'verbose'),    opts.verbose    = true;      end

if ischar(geom) || isstring(geom)
    src  = char(geom);
    geom = load(src);
    if opts.verbose, fprintf('Loaded %s\n', src); end
end

if ~isfield(opts,'meshes')
    fn = fieldnames(geom);
    opts.meshes = fn(startsWith(fn,'mesh_') & ...
        cellfun(@(f) isstruct(geom.(f)), fn))';
end

if exist('meshcheckrepair','file') ~= 2
    error(['iso2mesh is not on the path — meshcheckrepair is required. ' ...
           'Add iso2mesh before running this.']);
end

co = opts.check_opts;
co.verbose = false;

for stage = 1:opts.max_stage
    R = cr_check_geometry(geom, co);
    if R.ok
        if opts.verbose
            fprintf('Geometry passes after %d repair stage(s).\n', stage-1);
        end
        break;
    end

    if opts.verbose
        fprintf('Stage %d: %s\n', stage, stage_name(stage));
    end

    for m = 1:numel(opts.meshes)
        nm = opts.meshes{m};
        if ~isfield(geom, nm), continue; end
        [V, F] = get_mesh(geom.(nm));
        n0 = size(V,1);

        try
            switch stage
                case 1
                    [V, F] = meshcheckrepair(V, F, 'dup');
                case 2
                    [V, F] = meshcheckrepair(V, F, 'deep');
                case 3
                    [V, F] = meshcheckrepair(V, F, 'meshfix');
                case 4
                    % Isotropic remesh — this is what removes slivers
                    [V, F] = meshresample(V, F, opts.keepratio);
                    [V, F] = meshcheckrepair(V, F, 'deep');
            end
        catch err
            if opts.verbose
                fprintf('  %-18s stage %d failed: %s\n', nm, stage, err.message);
            end
            continue;
        end

        geom.(nm) = set_mesh(geom.(nm), V, F);
        if opts.verbose
            fprintf('  %-18s %6d -> %6d verts\n', nm, n0, size(V,1));
        end
    end
end

R = cr_check_geometry(geom, co);

if opts.verbose
    if R.ok
        fprintf('\nRepaired: geometry should now tetrahedralise.\n');
    else
        fprintf('\nSTILL FAILING after %d stage(s):\n', opts.max_stage);
        for k = 1:numel(R.fatal), fprintf('  %s\n', R.fatal{k}); end
        fprintf(['\nIf this is one of a handful of warps, drop it and report\n' ...
                 'the reduced n. If most warps land here, reduce scale_range\n' ...
                 'in cr_generate_warps — the deformation is too aggressive\n' ...
                 'for the base mesh quality.\n']);
    end
end

if isfield(opts,'outfile') && ~isempty(opts.outfile)
    save(opts.outfile, '-struct', 'geom', '-v7.3');
    if opts.verbose, fprintf('Saved %s\n', opts.outfile); end
end

end


% LOCAL FUNCTIONS

function s = stage_name(k)
    names = {'remove duplicate/unused vertices', ...
             'fix orientation and open edges', ...
             'remove self-intersections', ...
             'isotropic remesh'};
    s = names{k};
end

function [V, F] = get_mesh(M)
    if isfield(M,'vertices'), V = M.vertices; else, V = M.v; end
    if isfield(M,'faces'),    F = M.faces;    else, F = M.f; end
    V = double(V); F = double(F(:,1:3));
end

function M = set_mesh(M, V, F)
    if isfield(M,'vertices'), M.vertices = V; else, M.v = V; end
    if isfield(M,'faces'),    M.faces    = F; else, M.f = F; end
end
