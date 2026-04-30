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

%% Ready for BEM/FEM forward modelling
% Pass the following to your forward modelling pipeline (msg_fwd):
%   - all_meshes  (registered simulation meshes)
%   - spine_sources (spinal cord source locations)
%   - all_meshes.transform (transformation matrix)


