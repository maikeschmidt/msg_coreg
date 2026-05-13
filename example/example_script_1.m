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

fprintf('All geometry files saved to: %s\n', savepath);
fprintf('Ready to run BEM leadfields via run_bem_leadfields in msg_fwd\n');

