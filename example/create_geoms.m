%% script to generate the meshes and other input for BEM and FEM forward modelling
% just need to comment/ uncomment certain bits to switch between BEM and
% FEM outputs

clearvars
close all
cd('D:\');
Metadata;
cd('D:\Simulations');

%% Import mesh 
mesh2 = ; %import the subject mesh from which you want to create your mesh models from
subject = mesh2;

p = [];
p.vertices = subject.pos;
p.faces = subject.tri;
p2 = reducepatch(p,0.1);
subject.pos = p2.vertices;
subject.tri = p2.faces;

mesh = [];
mesh.vertices = subject.pos;
mesh.faces = subject.tri;
mesh.unit = subject.unit;

% get the fiducuals of the torso - left shoulder, right shoulder, bottom of
% spine
sub_fids_select = spm_mesh_select(mesh);
sub_fids = sub_fids_select';

S = [];
S.subject = mesh; 
S.fiducials = sub_fids;
S.plot = 0;

T = cr_register_torso(S); %creates transform matrix to put canonical torso in the right space
%% load spine stl

%if using one of the provided generic ones skip this section and the next

%if using own spine but need to generate a bone load in cord here and
%follow steps to create bone mesh

%if using own cord and bone meshes - skip this section and next but dont
%forget to apply the transform mation Matrix T to both meshes!

spine = ft_read_headshape('D:\Co-Registration\msg_coreg\cervical_spine.stl', 'unit', 'mm'); %load a custom spine - this needsa to be in the MRI space!!!!

spine_data = [];
spine_data.vertices = spine.pos;
spine_data.faces = spine.tri;
spine_data.unit = spine.unit;

s2 = reducepatch(spine_data, 0.1);
s2.unit = spine_data.unit;

cord_mesh = [];
cord_mesh.faces = s2.faces;
cord_mesh.vertices = s2.vertices;
cord_mesh.unit = spine_data.unit;

y_min = min(cord_mesh.vertices(:,2));
y_max = max(cord_mesh.vertices(:,2));

S = [];
S.spine = cord_mesh;
S.T = T;
S.resolution = 20; 
S.ylim = [y_min y_max];
S.unit = cord_mesh.unit;
src_grid = fem_generate_spine_center(S);

%% Create Bone Mesh

% if using a custom spine - create a toroidal bone model for this

num_cylinders = size(src_grid.pos, 1);
outer_rad = max(src_grid.ydist) + (.1 *  max(src_grid.ydist)); %starting outer radius of torus is max diam of the cord + 10%
height = S.resolution - (0.2 * S.resolution); % height of each bone torus = 20% less than the spacing between centreline points

tt_add_bem; 

temp = {};

for ii = 1:num_cylinders
    origin = src_grid.pos(ii, :); 
    if ii == num_cylinders
        tangent = src_grid.pos(ii, :) - src_grid.pos(ii-1, :); 
    else
        tangent = src_grid.pos(ii+1, :) - src_grid.pos(ii, :); 
    end

    y_dist = src_grid.ydist(ii);
    cord_rad = (y_dist * .05) + y_dist/2;
    thickness = outer_rad - cord_rad;

    M = generate_torus(cord_rad, outer_rad, height, [0, 0, 0]); % Centered at origin initially

    rotated_vertices = rotate_torus(M.vertices, tangent);
    M.vertices = rotated_vertices + origin; 

    temp{ii} = M;

    outer_rad = outer_rad *1.005; %thickness increases by 5% every time
    height = height * 1.005; %height increases by 5% every time
end

combined_vertices = [];
combined_faces = [];

for i = 2:length(temp)
    mesh_temp = temp{i};  
    vertex_offset = size(combined_vertices, 1);  

    combined_vertices = [combined_vertices; mesh_temp.vertices]; 
    combined_faces = [combined_faces; mesh_temp.faces + vertex_offset]; 
end

bone_mesh = struct('vertices', combined_vertices, 'faces', combined_faces);
bone_mesh.unit = cord_mesh.unit;

reduced_bone = reducepatch(bone_mesh, 0.1);
reduced_bone.unit = 'mm';


figure;
hold on;
ft_plot_mesh(cord_mesh, 'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.5);
ft_plot_mesh(reduced_bone, 'facecolor', 'yellow', 'edgecolor', 'none', 'facealpha', 0.7);
axis equal;
grid on;
lighting gouraud;
camlight;  

%% Try and plot the model

S = [];
S.subject = mesh;
S.T = T;
S.spine_mode = 'cervical';  % can be 'default', 'cervical', 'custom

% for 'custom' mode:
% S.spine = user_spine_mesh;
% S.vertebrae = user_bone_mesh;

all_meshes = cr_check_registration(S);
all_meshes.unit = 'mm';


%save for export later 
mesh_heart = all_meshes.mesh_heart;
mesh_lungs = all_meshes.mesh_lungs;
mesh_torso = all_meshes.mesh_torso;
mesh_spine = all_meshes.mesh_spine;
mesh_bone = all_meshes.mesh_vertebrae;

%% source points
y_min = min(mesh_spine.vertices(:,2));
y_max = max(mesh_spine.vertices(:,2));

S = [];
S.spine = mesh_spine;
S.T = T;
S.resolution = 5; 
S.ylim = [y_min y_max];
S.unit = 'mm';
sources_center_line = cr_generate_spine_center(S);


% Plot the spine grid
figure;
hold on;
ft_plot_mesh(mesh_spine, 'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.5)
ft_plot_mesh(mesh_bone, 'facecolor', 'yellow', 'edgecolor', 'none', 'facealpha', 0.5)
scatter3(sources_center_line.pos(:,1),sources_center_line.pos(:,2),sources_center_line.pos(:,3),'black.')
hold off;

%% Save for FEM and BEM

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

back_coils_3axis.positions = back_sensors.coilpos;
back_coils_3axis.orientations = back_sensors.coilori;

front_coils_3axis.positions = front_sensors.coilpos;
front_coils_3axis.orientations = front_sensors.coilori;

sources_cent = sources_center_line;

% Save to .mat
cd('D:\Co-Registration\');
save('geometries_trial.mat', 'mesh_wm', 'mesh_bone', 'mesh_heart',...
    'mesh_lungs', 'mesh_torso', 'back_coils_3axis', 'front_coils_3axis','sources_cent', 'transform_matrix');



