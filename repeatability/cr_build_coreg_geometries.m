% cr_build_coreg_geometries - Build one msg_fwd geometry file per saved
%                             coregistration repeat
%
% Replays every transform stored by cr_repeat_coreg through the standard
% mesh registration pipeline and writes a geometry .mat per repeat, in the
% exact format run_bem_leadfields.m and run_fem_leadfields.m expect.
%
% No interactive fiducial selection happens here — the saved transforms
% are supplied via S.transform, so this runs unattended.
%
% USAGE:
%   S.subject     = subject_surface_mesh;
%   S.repeat_file = 'coreg_repeats.mat';
%   S.outdir      = 'D:\Simulations\Coreg\geometries';
%   files = cr_build_coreg_geometries(S);
%
% INPUT:
%   S - struct with fields:
%
%   Required:
%     S.subject      - Subject surface mesh (.vertices, .faces)
%     S.repeat_file  - .mat saved by cr_repeat_coreg (contains R), OR
%                      S.transforms as a {1 x N} cell of 4x4 matrices
%     S.outdir       - Output folder for the geometry .mat files
%
%   Optional:
%     S.bone_modes   - Cell array of bone variants to build per repeat.
%                      Default: {'cont', 'inhomo', 'realistic'} — the three
%                      the paper compares. Note 'realistic' is only
%                      available for anatomical torso mode.
%     S.torso_mode   - 'canonical' (default) | 'anatomical'
%                      Canonical is the point of this analysis.
%     S.spine_mode   - 'full' (default) | 'cervical'
%     S.resolution   - Source spacing along the cord in mm (default: 5)
%     S.sensor_gen   - Generate sensor arrays (default: true). Required if
%                      you do not have an experimental array to reuse.
%     S.sensors      - Experimental sensor array to embed instead of
%                      generating one. When supplied it is saved as
%                      experimental_sensors and sensor_gen is ignored.
%     S.prefix       - Filename prefix (default: 'geometries_coreg')
%     S.overwrite    - Rebuild files that already exist (default: false)
%
% OUTPUT:
%   files - cell array of the geometry .mat paths written
%
% OUTPUT FILE CONTENTS (matches msg_fwd expectations exactly):
%   mesh_wm, mesh_bone, mesh_heart, mesh_lungs, mesh_torso
%   sources_cent                (with .pos and .inside)
%   front_coils_3axis, back_coils_3axis   [or experimental_sensors]
%   coreg_transform             the 4x4 replayed for this file
%   coreg_repeat_index          which repeat this came from
%
% FILENAMES:
%   <prefix>_rep<NN>_<bone_mode>.mat
%   e.g. geometries_coreg_rep01_realistic.mat
%
% NEXT STEPS:
%   Point run_bem_leadfields.m and run_fem_leadfields.m at S.outdir and
%   add the printed filenames to their `filenames` list, then run
%   msg_fwd/stats/st_collect_replicates to pool the results.
%
% SEE ALSO:
%   cr_repeat_coreg, cr_summarise_coreg, cr_check_registration,
%   cr_generate_spine_center
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

function files = cr_build_coreg_geometries(S)

if ~isfield(S, 'subject'), error('Please provide S.subject.'); end
if ~isfield(S, 'outdir'),  error('Please provide S.outdir.');  end

if ~isfield(S, 'bone_modes'), S.bone_modes = {'cont','inhomo','realistic'}; end
if ~isfield(S, 'torso_mode'), S.torso_mode = 'canonical'; end
if ~isfield(S, 'spine_mode'), S.spine_mode = 'full';      end
if ~isfield(S, 'resolution'), S.resolution = 5;           end
if ~isfield(S, 'sensor_gen'), S.sensor_gen = true;        end
if ~isfield(S, 'prefix'),     S.prefix     = 'geometries_coreg'; end
if ~isfield(S, 'overwrite'),  S.overwrite  = false;       end

% Realistic bone meshes do not exist for the canonical torso
if strcmpi(S.torso_mode, 'canonical') && any(strcmpi(S.bone_modes, 'realistic'))
    warning(['Realistic bone is not available for canonical torso mode. ' ...
             'Dropping it from bone_modes.']);
    S.bone_modes = S.bone_modes(~strcmpi(S.bone_modes, 'realistic'));
end

% Gather the transforms
if isfield(S, 'transforms') && ~isempty(S.transforms)
    transforms = S.transforms;
else
    if ~isfield(S, 'repeat_file')
        error('Provide either S.transforms or S.repeat_file.');
    end
    tmp        = load(S.repeat_file, 'R');
    transforms = tmp.R.transforms;
end

n_rep = numel(transforms);
if n_rep == 0
    error('No transforms found to build from.');
end

if ~exist(S.outdir, 'dir'); mkdir(S.outdir); end

fprintf('\n=== Building geometries from %d coregistration repeat(s) ===\n', n_rep);
fprintf('Bone modes : %s\n', strjoin(S.bone_modes, ', '));
fprintf('Output     : %s\n\n', S.outdir);

files = {};

for k = 1:n_rep
    for b = 1:numel(S.bone_modes)

        bone_mode = S.bone_modes{b};
        fname     = sprintf('%s_rep%02d_%s.mat', S.prefix, k, bone_mode);
        outfile   = fullfile(S.outdir, fname);

        if isfile(outfile) && ~S.overwrite
            fprintf('  [%2d/%s] exists, skipping: %s\n', k, bone_mode, fname);
            files{end+1} = outfile; %#ok<AGROW>
            continue;
        end

        fprintf('  [%2d/%s] building %s ...\n', k, bone_mode, fname);

        % Replay the saved transform — no interactive selection
        C            = struct();
        C.subject    = S.subject;
        C.transform  = transforms{k};
        C.torso_mode = S.torso_mode;
        C.spine_mode = S.spine_mode;
        C.bone_mode  = bone_mode;

        if isfield(S, 'sensors') && ~isempty(S.sensors)
            C.sensors    = S.sensors;
            C.sensor_gen = false;
        else
            C.sensor_gen = S.sensor_gen;
        end

        meshes = cr_check_registration(C);
        close(gcf);   % cr_check_registration always draws a check figure

        % Source model along the cord centreline
        y_min = min(meshes.spine.vertices(:, 2));
        y_max = max(meshes.spine.vertices(:, 2));

        G            = struct();
        G.spine      = meshes.spine;
        G.resolution = S.resolution;
        G.ylim       = [y_min y_max];
        G.unit       = 'mm';

        sources_cent = cr_generate_spine_center(G);

        % Assemble in the exact field layout msg_fwd expects.
        % 'wm' is the spinal cord mesh (white matter compartment).
        geom              = struct();
        geom.mesh_wm      = meshes.spine;
        geom.mesh_bone    = meshes.bone;
        geom.mesh_heart   = meshes.heart;
        geom.mesh_lungs   = meshes.lungs;
        geom.mesh_torso   = meshes.torso;
        geom.sources_cent = sources_cent;

        if isfield(S, 'sensors') && ~isempty(S.sensors)
            geom.experimental_sensors = S.sensors;
        else
            if isfield(meshes, 'back_sensors')
                geom.back_coils_3axis = meshes.back_sensors;
            end
            if isfield(meshes, 'front_sensors')
                geom.front_coils_3axis = meshes.front_sensors;
            end
        end

        % Provenance
        geom.coreg_transform    = transforms{k};
        geom.coreg_repeat_index = k;
        geom.bone_mode          = bone_mode;
        geom.torso_mode         = S.torso_mode;

        save(outfile, '-struct', 'geom', '-v7.3');
        files{end+1} = outfile; %#ok<AGROW>
        fprintf('          saved (%d sources)\n', size(sources_cent.pos, 1));
    end
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
