% cr_build_warp_geometries - Build one msg_fwd geometry file per anatomical warp
%
% Applies each warp from cr_generate_warps to the full anatomical model —
% every mesh, the source positions, AND the sensor array — and writes a
% geometry .mat per warp in the format run_bem_leadfields.m and
% run_fem_leadfields.m expect.
%
% SENSOR HANDLING
%   The sensor array is warped by the SAME matrix as the anatomy. This
%   keeps the sensor count and the source-to-sensor correspondence
%   identical across all warps, so leadfields stay row-matched and the
%   group statistics in msg_fwd/stats can pair sources across replicates
%   directly.
%
%   The trade-off, worth stating whenever these geometries are used:
%   warping the array distorts inter-sensor spacing, so these are not
%   physically realisable rigid arrays. That is fine for comparisons made
%   WITHIN each warped geometry, where both models see exactly the same
%   sensors. It is NOT fine for comparing absolute field amplitudes across
%   warps as though they came from a real array.
%
%   Sensor ORIENTATIONS (coilori/chanori) are rotated by the warp's linear
%   part using the inverse-transpose, which is the correct transformation
%   for directions under a non-rigid affine map, then renormalised to unit
%   length. Positions use the full transform.
%
% USAGE:
%   S.geom_file = 'geometries_original_source_original.mat';
%   S.warp_file = 'anatomical_warps.mat';
%   S.outdir    = 'D:\Simulations\Warps\geometries';
%   files = cr_build_warp_geometries(S);
%
% INPUT:
%   S - struct with fields:
%
%   Required:
%     S.geom_file - Path to an EXISTING geometry .mat to warp (the
%                   original anatomical model). Must contain mesh_wm,
%                   mesh_bone, mesh_heart, mesh_lungs, mesh_torso and
%                   sources_cent.
%     S.outdir    - Output folder for warped geometry .mat files
%
%   Optional:
%     S.warp_file - .mat from cr_generate_warps (default: 'anatomical_warps.mat')
%     S.warps     - Warp struct W directly, instead of loading a file
%     S.prefix    - Filename prefix (default: 'geometries_warp')
%     S.variant_tag - Bone-variant tag appended to each filename, so the
%                   output matches the naming the leadfield runners and
%                   msg_fwd/stats expect: geometries_warp01_realistic.
%                   Defaults to the variant inferred from S.geom_file
%                   ('realistic' | 'inhomo' | 'homo' | 'cont'). Pass ''
%                   to omit the tag entirely.
%                   IMPORTANT: a geometry file carries ONE bone mesh, so to
%                   get both a realistic and a continuous model per warp you
%                   must run this script twice, pointing S.geom_file at each
%                   base geometry in turn. The warp set is seeded, so warp k
%                   is the same deformation in both runs.
%     S.overwrite - Rebuild existing files (default: false)
%     S.recentre  - If true, re-derive the warp centre as the centroid of
%                   mesh_torso and rebuild each warp about that point
%                   (default: true). Strongly recommended — warping about
%                   the origin translates the model as well as reshaping it.
%
% OUTPUT:
%   files - cell array of geometry .mat paths written
%
% OUTPUT FILE CONTENTS:
%   Same fields as the input geometry, all warped, plus:
%     warp_matrix        the 4x4 applied
%     warp_index         which warp this came from
%     warp_scales        [1 x 3] per-axis scale factors
%     warp_label         e.g. 'taller / thinner'
%
% NEXT STEPS:
%   Add the printed filenames to run_bem_leadfields.m and
%   run_fem_leadfields.m, run both, then use
%   msg_fwd/stats/st_collect_replicates.
%
% SEE ALSO:
%   cr_generate_warps, cr_plot_warps
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_coreg
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Coregistration Toolbox.

function files = cr_build_warp_geometries(S)

if ~isfield(S, 'geom_file'), error('Please provide S.geom_file.'); end
if ~isfield(S, 'outdir'),    error('Please provide S.outdir.');    end

if ~isfield(S, 'warp_file'), S.warp_file = 'anatomical_warps.mat'; end
if ~isfield(S, 'prefix'),    S.prefix    = 'geometries_warp';      end
if ~isfield(S, 'overwrite'), S.overwrite = false;                  end
if ~isfield(S, 'recentre'),  S.recentre  = true;                   end

% Infer the bone variant from the base geometry filename unless told.
% Order matters: 'cont' is a substring of nothing here, but check the
% longer names first regardless so a rename cannot silently mismatch.
if ~isfield(S, 'variant_tag')
    known = {'realistic', 'inhomo', 'homo', 'cont'};
    [~, base_name] = fileparts(S.geom_file);
    S.variant_tag = '';
    for v = 1:numel(known)
        if contains(base_name, known{v})
            S.variant_tag = known{v};
            break;
        end
    end
    if isempty(S.variant_tag)
        warning(['Could not infer a bone variant from "%s". Output files ' ...
                 'will have no variant tag, which means msg_fwd/stats will ' ...
                 'not be able to tell bone models apart. Set S.variant_tag ' ...
                 'explicitly.'], base_name);
    else
        fprintf('Inferred bone variant "%s" from the base geometry name.\n', ...
            S.variant_tag);
    end
end

% Load warps
if isfield(S, 'warps') && ~isempty(S.warps)
    W = S.warps;
else
    tmp = load(S.warp_file, 'W');
    W   = tmp.W;
end

% Load the base geometry
if ~isfile(S.geom_file)
    error('Geometry file not found: %s', S.geom_file);
end
geom0 = load(S.geom_file);

required = {'mesh_wm','mesh_bone','mesh_heart','mesh_lungs','mesh_torso','sources_cent'};
missing  = required(~isfield(geom0, required));
if ~isempty(missing)
    error('Geometry file is missing required field(s): %s', strjoin(missing, ', '));
end

n = numel(W.matrices);

% Rebuild warps about the torso centroid so they reshape rather than
% translate the model.
if S.recentre
    centre = mean(geom0.mesh_torso.vertices, 1);
    fprintf('Recentring warps about torso centroid [%.1f %.1f %.1f] mm\n', centre);

    P            = W.params;
    P.centre     = centre;
    P.n_warps    = n;
    P.outfile    = fullfile(tempdir, 'anatomical_warps_recentred.mat');
    W            = cr_generate_warps(P);
end

if ~exist(S.outdir, 'dir'); mkdir(S.outdir); end

% Mesh and sensor fields to transform
mesh_fields = {'mesh_wm','mesh_bone','mesh_heart','mesh_lungs','mesh_torso'};
extra_mesh  = {'mesh_spine','mesh_back_muscle'};   % transformed if present
sens_fields = {'front_coils_3axis','back_coils_3axis','experimental_sensors', ...
               'coils_3axis','front_sensors','back_sensors','fullbody_sensors'};

fprintf('\n=== Building %d warped geometries ===\n', n);
fprintf('Base   : %s\n', S.geom_file);
fprintf('Output : %s\n\n', S.outdir);

files = {};

for k = 1:n

    if isempty(S.variant_tag)
        fname = sprintf('%s%02d.mat', S.prefix, k);
    else
        fname = sprintf('%s%02d_%s.mat', S.prefix, k, S.variant_tag);
    end
    outfile = fullfile(S.outdir, fname);

    if isfile(outfile) && ~S.overwrite
        fprintf('  [%2d] exists, skipping: %s\n', k, fname);
        files{end+1} = outfile; %#ok<AGROW>
        continue;
    end

    M    = W.matrices{k};
    geom = geom0;

    % Meshes
    for f = [mesh_fields, extra_mesh]
        fld = f{1};
        if isfield(geom, fld) && isfield(geom.(fld), 'vertices')
            geom.(fld).vertices = apply_T(M, geom.(fld).vertices);
        end
    end

    % Source positions
    geom.sources_cent.pos = apply_T(M, geom.sources_cent.pos);

    % Sensors — positions by M, orientations by inverse-transpose of the
    % linear part (the correct rule for directions under an affine map)
    for f = sens_fields
        fld = f{1};
        if isfield(geom, fld) && ~isempty(geom.(fld))
            geom.(fld) = warp_sensors(geom.(fld), M);
        end
    end

    % Provenance
    geom.warp_matrix = M;
    geom.warp_index  = k;
    geom.warp_scales = W.scales(k, :);
    geom.warp_label  = W.labels{k};

    save(outfile, '-struct', 'geom', '-v7.3');
    files{end+1} = outfile; %#ok<AGROW>

    fprintf('  [%2d] %s  scales [%.3f %.3f %.3f]  %s\n', ...
        k, fname, W.scales(k,:), W.labels{k});
end

fprintf('\n=== %d geometry file(s) written ===\n\n', numel(files));
fprintf('Paste this into the `filenames` list in run_bem_leadfields.m\n');
fprintf('and run_fem_leadfields.m:\n\n');
for i = 1:numel(files)
    [~, nm] = fileparts(files{i});
    fprintf('    ''%s'', ...\n', nm);
end
fprintf('\n');

end


% LOCAL FUNCTIONS

function p = apply_T(T, pts)
% Apply a 4x4 transform to [N x 3] points.
    p = (T * [pts, ones(size(pts,1),1)]')';
    p = p(:, 1:3);
end

function sens = warp_sensors(sens, M)
% Warp a FieldTrip sensor struct. Positions use the full transform;
% orientations use inv(A)' and are renormalised.

    A     = M(1:3, 1:3);
    Aori  = inv(A)';   % direction transform under an affine map

    pos_fields = {'coilpos','chanpos','elecpos','pnt'};
    ori_fields = {'coilori','chanori'};

    for f = pos_fields
        fld = f{1};
        if isfield(sens, fld) && ~isempty(sens.(fld)) && size(sens.(fld),2) == 3
            sens.(fld) = apply_T(M, sens.(fld));
        end
    end

    for f = ori_fields
        fld = f{1};
        if isfield(sens, fld) && ~isempty(sens.(fld)) && size(sens.(fld),2) == 3
            v = (Aori * sens.(fld)')';
            nrm = vecnorm(v, 2, 2);
            nrm(nrm < eps) = 1;          % leave degenerate rows unscaled
            sens.(fld) = v ./ nrm;
        end
    end
end
