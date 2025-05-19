%% script for spinal cord BEM forward modelling

clearvars
close all
clc

cd('D:\');
Metadata;
proj_init;
cd('D:\Simulations\');
%% Load in all generated meshes/ inputs

geoms_path = ('D:\Simulations\BEM\geometries' );
geoms = load(fullfile(geoms_path, 'geometries_19seg_bem.mat'));
ordering = {'wm','bone', 'heart', 'lungs', 'torso'}; 

clear mesh
mesh_idx = 1;

tt_add_bem; 

for ii = 1:numel(ordering)

    field = ['mesh_' ordering{ii}];

    tmp = [];
    tmp.faces         = geoms.(field).faces;
    tmp.vertices         = geoms.(field).vertices;
    tmp.unit        = 'm'; 
    tmp.name        = ordering{ii};

    orient = hbf_CheckTriangleOrientation(tmp.vertices,tmp.faces);
        if orient == 2
            tmp.faces = tmp.faces(:,[1 3 2]);
        end

    mesh(ii) = tmp;

end

cord_mesh = mesh(1);
bone_mesh = mesh(2); %remove bone here
heart_mesh = mesh(3);
lung_mesh = mesh(4);
torso_mesh = mesh(5);

%% Generate source and grad structures
src = [];
src.pos = geoms.sources_cent.pos;
src.inside = ones(length(src.pos),1);
src.unit = 'm'; 

% Generate the grad structure for the back
grad_back = [];
grad_back.coilpos = geoms.back_coils_3axis.positions;
grad_back.coilori = geoms.back_coils_3axis.orientations;
grad_back.tra = eye(length(grad_back.coilpos));
for ii = 1:length(grad_back.coilpos)
    grad_back.label{ii} = sprintf('Chan-%03d',ii);
end
grad_back.unit = 'm';

% Generate the grad structure for the front
grad_front = [];
grad_front.coilpos = geoms.front_coils_3axis.positions;
grad_front.coilori = geoms.front_coils_3axis.orientations;
grad_front.tra = eye(length(grad_front.coilpos));
for ii = 1:length(grad_front.coilpos)
    grad_front.label{ii} = sprintf('Chan-%03d',ii);
end
grad_front.unit = 'm';

%% lets do some BEM modelling!

i = 19; %which source point from the grid are we modelling for?

cratio = 40;

% Forward model using BEM
S_forward = [];
S_forward.pos = single_source_point;   
S_forward.T = geoms.transform_matrix;                      
S_forward.ori = [1,0,0]; %1 nA/m --> [1,0,0] = left to right, [0,1,0] = anerior to posterior, [0,0,1] = superior to inferior
% S_forward.ori = [0,1,0];
% S_forward.ori = [0,0,1];
S_forward.posunits = 'm';         
S_forward.names = {'spinalcord','vertebrae','blood','lungs','torso'};  
S_forward.ci = [0.33 .33/cratio .62 .05 .23]; 
S_forward.co = [.23 .23 .23 .23 0 ]; 
S_forward.cord = cord_mesh;
S_forward.vertebrae = bone_mesh; 

S_back = S_forward;
S_back.sensors = grad_back; 
Ls_back = tt_fwds_bem5(S_back); 

S_front = S_forward;
S_front.sensors = grad_front;
Ls_front = tt_fwds_bem5(S_front); 

%% split Ls_ structure into its relevent parts
Ls_back_xyz = reshape(Ls_back, [size(grad_back.label, 2)/3 , 3]);
Ls_front_xyz = reshape(Ls_front, [size(grad_front.label, 2)/3 , 3]);

Ls_back_x = Ls_back_xyz(:, 1);
Ls_back_y = Ls_back_xyz(:, 2);
Ls_back_z = Ls_back_xyz(:, 3);

Ls_front_x = Ls_front_xyz(:, 1);
Ls_front_y = Ls_front_xyz(:, 2);
Ls_front_z = Ls_front_xyz(:, 3);

cd('D:\Simulations\BEM\front_v_back\inside_segments');
save('back_20_ori1.mat', 'Ls_back_x', 'Ls_back_y' , 'Ls_back_z');
% save('back_20_ori2.mat', 'Ls_back_x', 'Ls_back_y' , 'Ls_back_z');
% save('back_20_ori3.mat', 'Ls_back_x', 'Ls_back_y' , 'Ls_back_z');
save('front_20_ori1.mat', 'Ls_front_x', 'Ls_front_y', 'Ls_front_z');
% save('front_20_ori2.mat', 'Ls_front_x', 'Ls_front_y', 'Ls_front_z');
% save('front_20_ori3.mat', 'Ls_front_x', 'Ls_front_y', 'Ls_front_z');

