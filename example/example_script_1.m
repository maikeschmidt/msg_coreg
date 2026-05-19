% example_script_1 - Register simulation meshes to experimental OPM sensor space
%
% Demonstrates how to co-register canonical simulation meshes with an
% existing experimental OPM sensor array and an optical body scan.
% Produces all mesh components and a spinal cord source model ready
% for BEM/FEM forward modelling.
%
% WORKFLOW:
%   1. Load and define the optical surface scan of the participant
%   2. Import the experimental OPM sensor array via SPM
%   3. Run cr_check_registration() to load and align all simulation meshes
%   4. Generate a spinal cord centreline source model
%
% USE THIS SCRIPT WHEN:
%   - You already have an experimentally defined OPM sensor layout
%   - You want to run simulations in the same coordinate space as your
%     recorded data
%   - Subject-specific MRI is unavailable (canonical model)
%
% REQUIREMENTS:
%   - SPM (developmental version recommended)
%   - FieldTrip
%   - Helsinki BEM Framework (hbf_lc_p subfolder)
%   - Optical/3D surface scan of participant (.stl or equivalent)
%   - OPM sensor position file (.tsv)
%   - OPM data file (.lvm)
%
% OUTPUTS (workspace variables):
%   all_meshes    - Struct of all registered simulation meshes
%   mesh_torso    - Torso surface mesh
%   mesh_spine    - Spine mesh
%   mesh_bone     - Bone mesh
%   mesh_heart    - Heart mesh
%   mesh_lungs    - Lung mesh
%   spine_sources - Spinal cord centreline source model
%                   (for distributed spinal source simulation)
%
% SEE ALSO:
%   cr_check_registration, cr_generate_spine_center, cr_add_functions
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
% Date:   April 2026
%
% This file is part of the MSG Coregistration Toolbox.


%% Initialise toolbox
cr_add_functions;

%% Load and define the optical surface scan
mesh = struct();
mesh.vertices = ;          % [N x 3] vertex coordinates
mesh.faces    = ;          % [M x 3] face indices
mesh.unit     = 'mm';      % must match the grad structure units

%% Import experimental OPM sensor array
sub         = '';
posfile_back = '.tsv';
savepath     = fullfile('D:\OPM\', ['sub-OP', sub(3:end)], 'nervesim');
cd(savepath)

D            = {};
filetemplate = sprintf('nervesim-right-run_002_array1.lvm');
S            = struct();
S.data       = filetemplate;
S.positions  = posfile_back;
D            = spm_opm_create(S);

exp_sensors  = sensors(D, 'MEG');

%% Load all simulation meshes
S             = [];
S.subject     = mesh;
S.sensors     = exp_sensors;
S.spine_mode  = 'full';       % 'full' or 'cervical'
S.torso_mode  = 'canonical';  % 'canonical' or 'anatomical'
S.bone_mode   = 'cont';       % 'cont', 'homo', 'inhomo', or 'realistic'
% S.brain     = true;         % uncomment to include brain registration

all_meshes  = cr_check_registration(S);
mesh_torso  = all_meshes.torso;
mesh_spine  = all_meshes.spine;
mesh_bone   = all_meshes.bone;
mesh_heart  = all_meshes.heart;
mesh_lungs  = all_meshes.lungs;
% mesh_brain  = all_meshes.brain;   % uncomment if S.brain = true
% mesh_skull  = all_meshes.oskull;  % uncomment if S.brain = true

%% Generate spinal cord centreline source model
y_min = min(mesh_spine.vertices(:,2));
y_max = max(mesh_spine.vertices(:,2));

S            = [];
S.spine      = mesh_spine;
S.resolution = 5;
S.ylim       = [y_min y_max];
S.unit       = 'mm';

spine_sources = cr_generate_spine_center(S);

%% CREATE SHIFTED SOURCE MODELS FOR SENSITIVITY ANALYSIS
% Generates 18 shifted versions of the spinal cord source model by
% translating all source positions independently along each axis (X, Y, Z).
% Together with the original, this gives 19 geometry configurations
% for BEM leadfield computation.
%
% Shift magnitudes (mm), applied independently per axis:
%   X axis (Left-Right):      ±2, ±4, ±6 mm
%   Y axis (Rostral-Caudal):  ±2, ±4, ±6 mm
%   Z axis (Ventral-Dorsal):  ±2, ±4, ±6 mm

% Each row: [dx, dy, dz] — only one axis non-zero per shift
shift_vectors_mm = [
    % X axis shifts (Left-Right)
     2,  0,  0;
     4,  0,  0;
     6,  0,  0;
    -2,  0,  0;
    -4,  0,  0;
    -6,  0,  0;
    % Y axis shifts (Rostral-Caudal)
     0,  2,  0;
     0,  4,  0;
     0,  6,  0;
     0, -2,  0;
     0, -4,  0;
     0, -6,  0;
    % Z axis shifts (Ventral-Dorsal)
     0,  0,  2;
     0,  0,  4;
     0,  0,  6;
     0,  0, -2;
     0,  0, -4;
     0,  0, -6;
];

shift_labels = {
    % X axis
    'shift_x_pos2mm', 'shift_x_pos4mm', 'shift_x_pos6mm', ...
    'shift_x_neg2mm', 'shift_x_neg4mm', 'shift_x_neg6mm'; ...
    % Y axis
    'shift_y_pos2mm', 'shift_y_pos4mm', 'shift_y_pos6mm', ...
    'shift_y_neg2mm', 'shift_y_neg4mm', 'shift_y_neg6mm'; ...
    % Z axis
    'shift_z_pos2mm', 'shift_z_pos4mm', 'shift_z_pos6mm', ...
    'shift_z_neg2mm', 'shift_z_neg4mm', 'shift_z_neg6mm'; ...
};

% Save original geometry as the reference
geom_original = struct();
geom_original.mesh_wm               = mesh_wm;
geom_original.mesh_bone             = mesh_bone;
geom_original.mesh_heart            = mesh_heart;
geom_original.mesh_lungs            = mesh_lungs;
geom_original.mesh_torso            = mesh_torso;
geom_original.sources_cent          = spine_sources;
geom_original.experimental_sensors  = exp_sensors;

outfile_original = fullfile(savepath, 'geometries_original.mat');
save(outfile_original, '-struct', 'geom_original', '-v7.3');
fprintf('Saved: geometries_original.mat\n');

% Generate and save each shifted geometry
for s = 1:size(shift_vectors_mm, 1)

    shift_mm = shift_vectors_mm(s, :);
    label    = shift_labels{s};

    fprintf('Creating shifted geometry: %s (shift = [%d %d %d] mm)...\n', ...
        label, shift_mm(1), shift_mm(2), shift_mm(3));

    % Shift source positions in mm
    shifted_sources     = spine_sources;
    shifted_sources.pos = spine_sources.pos + shift_mm;

    % Package geometry — meshes and sensors unchanged, only sources shift
    geom_shifted                        = struct();
    geom_shifted.mesh_wm                = mesh_wm;
    geom_shifted.mesh_bone              = mesh_bone;
    geom_shifted.mesh_heart             = mesh_heart;
    geom_shifted.mesh_lungs             = mesh_lungs;
    geom_shifted.mesh_torso             = mesh_torso;
    geom_shifted.sources_cent           = shifted_sources;
    geom_shifted.experimental_sensors   = exp_sensors;

    outfile = fullfile(savepath, ['geometries_' label '.mat']);
    save(outfile, '-struct', 'geom_shifted', '-v7.3');
    fprintf('  Saved: geometries_%s.mat\n', label);
end

%% CREATE SHIFTED SENSOR ARRAY MODELS FOR SENSITIVITY ANALYSIS
% Generates sensor-shifted geometry files by translating the entire
% experimental OPM array by a random 3D displacement [dx, dy, dz],
% where each axis component is drawn independently from a uniform
% distribution centred around a target magnitude.
%
% THREE BUNDLES of 8 shifts, one per registration error scale:
%   Bundle 1 — small errors:   each axis shift drawn from U(-3, -1) ∪ U(1, 3) mm
%   Bundle 2 — medium errors:  each axis shift drawn from U(-7, -3) ∪ U(3, 7) mm
%   Bundle 3 — large errors:   each axis shift drawn from U(-13,-7) ∪ U(7,13) mm
%
% Each shift is a unique [dx, dy, dz] vector — all three axes move
% simultaneously but by independently drawn amounts. Signs are random
% so shifts can be in either direction along each axis.
%
% 24 shifted configurations + 1 original = 25 geometry files total.
%
% Sensor orientations (coilori, chanori) and transfer matrix (tra)
% are NOT modified — only coilpos and chanpos are shifted so the
% triaxial orthogonal structure remains intact.
%
% REPRODUCIBILITY:
%   Random seed is set below. Change the seed for a different random
%   realisation, or uncomment the hardcoded block to fix specific values.

% Random seed 
rng(42);   % SET THIS: change seed for a different random realisation

% Bundle definitions
% Each bundle: [lower_bound, upper_bound] for the magnitude of each axis
% shift. Signs are applied randomly after magnitude is drawn.
bundle_names  = {'small_2mm', 'medium_5mm', 'large_10mm'};
bundle_ranges = [1, 3;    % Bundle 1: ~2mm  — uniform U(1,3) mm magnitude
                 3, 7;    % Bundle 2: ~5mm  — uniform U(3,7) mm magnitude
                 7, 13];  % Bundle 3: ~10mm — uniform U(7,13) mm magnitude
n_bundles     = 3;
n_shifts      = 8;   % shifts per bundle

% Generate shift vectors
% shift_vectors{b}(s,:) = [dx, dy, dz] in mm for bundle b, shift s
shift_vectors = cell(1, n_bundles);

for b = 1:n_bundles
    lo  = bundle_ranges(b, 1);
    hi  = bundle_ranges(b, 2);
    vecs = zeros(n_shifts, 3);

    for s = 1:n_shifts
        % Draw magnitude for each axis independently from U(lo, hi)
        magnitudes = lo + (hi - lo) * rand(1, 3);

        % Apply random sign independently to each axis
        signs      = sign(randn(1, 3));   % random +1 or -1 per axis
        signs(signs == 0) = 1;            % handle exact zero (vanishingly rare)

        vecs(s, :) = magnitudes .* signs;
    end

    shift_vectors{b} = vecs;

    fprintf('\nBundle %d (%s) — shift vectors [dx, dy, dz] in mm:\n', ...
        b, bundle_names{b});
    for s = 1:n_shifts
        fprintf('  Shift %d: [%+.2f, %+.2f, %+.2f] mm\n', ...
            s, vecs(s,1), vecs(s,2), vecs(s,3));
    end
end

% Option: hardcode specific shift vectors for exact reproducibility
% Paste the values printed above to lock in specific shifts.
% Uncomment and fill in to bypass random generation entirely.
%
% shift_vectors{1} = [   % Bundle 1 — small (~2mm)
%    +1.23, -2.01, +1.87;
%    -1.56, +1.34, -2.45;
%    +2.11, -1.78, +1.23;
%    -1.89, +2.34, -1.67;
%    +1.45, -1.23, +2.78;
%    -2.34, +1.89, -1.45;
%    +1.67, -2.56, +1.34;
%    -2.01, +1.67, -2.23;
% ];
% shift_vectors{2} = [ ... ];   % Bundle 2 — medium (~5mm)
% shift_vectors{3} = [ ... ];   % Bundle 3 — large (~10mm)

% Save original geometry as reference 
geom_sensor_original = struct();
geom_sensor_original.mesh_wm               = mesh_wm;
geom_sensor_original.mesh_bone             = mesh_bone;
geom_sensor_original.mesh_heart            = mesh_heart;
geom_sensor_original.mesh_lungs            = mesh_lungs;
geom_sensor_original.mesh_torso            = mesh_torso;
geom_sensor_original.sources_cent          = spine_sources;
geom_sensor_original.experimental_sensors  = exp_sensors;
geom_sensor_original.shift_type            = 'sensor';
geom_sensor_original.shift_vectors         = shift_vectors;   % save for reference
geom_sensor_original.bundle_names          = bundle_names;
geom_sensor_original.bundle_ranges_mm      = bundle_ranges;

outfile_sensor_original = fullfile(savepath, 'geometries_sensor_original.mat');
save(outfile_sensor_original, '-struct', 'geom_sensor_original', '-v7.3');
fprintf('\nSaved: geometries_sensor_original.mat\n');

% Generate and save each shifted geometry 
for b = 1:n_bundles
    for s = 1:n_shifts
        shift_vec = shift_vectors{b}(s, :);   % [dx, dy, dz] in mm

        fprintf('Creating: bundle %d (%s), shift %d  [%+.2f, %+.2f, %+.2f] mm...\n', ...
            b, bundle_names{b}, s, shift_vec(1), shift_vec(2), shift_vec(3));

        % Apply shift to sensor positions only
        % coilpos and chanpos shifted by same vector
        % coilori, chanori, tra, balance left completely untouched
        % so triaxial orthogonality is fully preserved
        shifted_sensors         = exp_sensors;
        shifted_sensors.coilpos = exp_sensors.coilpos + shift_vec;
        shifted_sensors.chanpos = exp_sensors.chanpos + shift_vec;

        % Package geometry
        geom_shifted                              = struct();
        geom_shifted.mesh_wm                      = mesh_wm;
        geom_shifted.mesh_bone                    = mesh_bone;
        geom_shifted.mesh_heart                   = mesh_heart;
        geom_shifted.mesh_lungs                   = mesh_lungs;
        geom_shifted.mesh_torso                   = mesh_torso;
        geom_shifted.sources_cent                 = spine_sources;
        geom_shifted.experimental_sensors         = shifted_sensors;
        geom_shifted.shift_type                   = 'sensor';
        geom_shifted.shift_vec_mm                 = shift_vec;
        geom_shifted.bundle_name                  = bundle_names{b};
        geom_shifted.bundle_idx                   = b;
        geom_shifted.shift_idx                    = s;

        label   = sprintf('bundle%d_shift%d', b, s);
        outfile = fullfile(savepath, ['geometries_sensor_' label '.mat']);
        save(outfile, '-struct', 'geom_shifted', '-v7.3');
        fprintf('  Saved: geometries_sensor_%s.mat\n', label);
    end
end

fprintf('\nAll sensor-shifted geometry files saved to: %s\n', savepath);
fprintf('%d configurations total: 1 original + %d shifted (%d bundles x %d shifts)\n', ...
    1 + n_bundles * n_shifts, n_bundles * n_shifts, n_bundles, n_shifts);

% Print filenames for run_bem_leadfields 
fprintf('\nFilenames for run_bem_leadfields.m:\n');
fprintf("  'geometries_sensor_original'\n");
for b = 1:n_bundles
    for s = 1:n_shifts
        fprintf("  'geometries_sensor_bundle%d_shift%d'\n", b, s);
    end
end
