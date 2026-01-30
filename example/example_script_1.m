%% Example script to co-register between OPM sensor positions, Optical Scan, and Simulation mesh models
cr_add_functions;

% Read and preprocess the optical scan

mesh = struct();
mesh.vertices = ;
mesh.faces = ;
mesh.unit = 'mm'; % Define unit explicitly - make sure it matches the grad structure!

% call in the experimental opm array
sub = '';
posfile_back = '.tsv';
savepath = fullfile('D:\OPM\', ['sub-OP', sub(3:end)], 'nervesim');
cd(savepath)
 
D = {}; 
filetemplate = sprintf('nervesim-right-run_002_array1.lvm');
S = struct(); 
S.data = filetemplate;
S.positions = posfile_back; 
D = spm_opm_create(S); 
 
exp_sensors = sensors(D, 'MEG');

%% Load all simulation meshes + Create Sensor array

cd('D:\Simulations\Paper_1\but_actualy\geometries')

S = []; 
S.subject = mesh;
S.sensors = exp_sensors;
S.spine_mode = 'full'; %which spine model to load (full or cervical)
S.torso_mode = 'canonical';
S.bone_mode = 'cont'; %type of bone model to load (cont, homo, inhomo, realistic)
% S.brain = true;
all_meshes = cr_check_registration(S);

mesh_torso = all_meshes.torso;
mesh_spine = all_meshes.spine;
mesh_bone = all_meshes.bone;
mesh_heart = all_meshes.heart;
mesh_lungs = all_meshes.lungs;
% mesh_brain = all_meshes.brain;
% mesh_skull = all_meshes.oskull;

%% create spine sources
y_min = min(mesh_spine.vertices(:,2));
y_max = max(mesh_spine.vertices(:,2));

S = [];
S.spine = mesh_spine;
S.resolution = 5; 
S.ylim = [y_min y_max];
S.unit = 'mm';
spine_sources = cr_generate_spine_center(S);

%% this is everything we need for BEM and FEM forward modelling!


