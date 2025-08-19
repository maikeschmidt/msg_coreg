%% Example script to co-register between OPM sensor positions, Optical Scan, and Simulation mesh models

% This pipeline is implemented for a subject dataset where we have the subject scan already in OPM sensor
% space, and hence we can skip registering the subject mesh to the opm sensors and go straight to pulling in
% the mesh models

%%

% First, need to create the spm opm structure to extract the sensor positions
cd('coreg_path\example')
clearvars
close all
clc

cd('D:\');
Metadata;

sub = 'OP00210';
posfile_back = 'D:\OPM\sub-OP00210\scannercast_positions\positions_all_OP00210.tsv';
savepath = fullfile('D:\OPM\', ['sub-OP', sub(3:end)], 'nervesim');
cd(savepath)

D = {}; 
filetemplate = sprintf('nervesim-right-run_002_array1.lvm');
S = struct(); 
S.data = filetemplate;
S.positions = posfile_back; 
D = spm_opm_create(S); 

%% Load and process the optical scan
cd('D:\OPM\scannercast')
optical = ft_read_headshape('surface.stl', 'unit', 'mm');
p = struct();
p.vertices = optical.pos;
p.faces = optical.tri;
p2 = reducepatch(p, 0.01); % Reduce complexity of the mesh
subject.pos = p2.vertices;
subject.tri = p2.faces;

% Prepare mesh structure
mesh = struct();
mesh.vertices = subject.pos;
mesh.faces = subject.tri;
mesh.unit = 'mm'; % Define unit explicitly - make sure it matches the grad structure!

%%
% now we need to register the simulation meshes to the optical scan which
% is in opm sensor space - but we know what the transform matrix is so can
% just go ahead and apply this
% 

S = []; 
S.subject = mesh;
S.sensors = sensors(D, 'MEG'); 
S.spine_mode = 'full'; %full or cervical
S.torso_mode = 'anatomical'; %anatomical or canonical (if choosing canonical - you will need to select three fiducials on the subject mesh (left shoulder, right shoulder, chin))
all_meshes = cr_check_registration(S);

%now we can save individual meshes in the space of opm sensor space
mesh_torso = all_meshes.mesh_torso;
mesh_spine = all_meshes.mesh_spine;
mesh_bone = all_meshes.mesh_vertebrae;
mesh_heart = all_meshes.mesh_heart;
mesh_lungs = all_meshes.mesh_lungs;

%create a source grid along the centerline of the spinal cord

y_min = min(mesh_spine.vertices(:,2));
y_max = max(mesh_spine.vertices(:,2));

S = [];
S.spine = mesh_spine;
S.T = T;
S.resolution = 5; 
S.ylim = [y_min y_max];
S.unit = 'mm';
sources_center_line = cr_generate_spine_center(S);

%from here we now have a transform matrix (T), source grid and all relevent meshes in opm sensor space - this is everything we need for BEM and FEM forward modelling!

transform_matrix = T; %need this for BEM 

mesh_heart.vertices = mesh_heart.vertices; 
mesh_heart.faces = mesh_heart.faces;

mesh_bone.vertices = mesh_bone.vertices;
mesh_bone.faces = mesh_bone.faces;

mesh_lungs.vertices = mesh_lungs.vertices; 
mesh_lungs.faces = mesh_lungs.faces;

mesh_torso.vertices = mesh_torso.vertices;
mesh_torso.faces = mesh_torso.faces;

mesh_wm.vertices = mesh_spine.vertices;
mesh_wm.faces = mesh_spine.faces;

coils_3axis = sensors(D,'meg');


sources_cent = sources_center_line;

% Save to .mat
cd('D:\OPM\sub-OP00210\fwd\');
save('geometries_new_cervical.mat', 'mesh_wm', 'mesh_bone', 'mesh_heart',...
    'mesh_lungs', 'mesh_torso', 'coils_3axis','sources_cent', 'transform_matrix');

