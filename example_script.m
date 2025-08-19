%% Example script to co-register between OPM sensor positions, Optical Scan, and Simulation mesh models

% This is early days, so there may be some snags in the pipeline - let me know!!!

clearvars
close all
cd('D:\');
Metadata;

% First, need to register the halo results to the optical scan OR create
% opm spm struct to extract the relevent reference opms

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

% halo = []; % Placeholder for halo file import if not using positions tsv
marker_opms = {'C2-X', 'D8-X','F5-X'}; % Identify which OPMs were used as the markers
m_opm_idx = find(contains(D.chanlabels,marker_opms)); % Extract the positions of these sensors
all_positions = sensors(D, 'MEG'); 
m_opm_pos = all_positions.chanpos(m_opm_idx, :);

% Read and preprocess the optical scan
optical = ft_read_headshape('D:\OPM\scannercast\sens_pos_fullbody_cast.stl', 'unit', 'mm');
p = struct();
p.vertices = optical.pos;
p.faces = optical.tri;
p2 = reducepatch(p, 0.005); % Reduce complexity of the mesh
subject.pos = p2.vertices;
subject.tri = p2.faces;

% Prepare mesh structure
mesh = struct();
mesh.vertices = subject.pos;
mesh.faces = subject.tri;
mesh.unit = 'mm'; % Define unit explicitly - make sure it matches the grad structure!

% Select the fiducials where the OPMs were placed
fids_select_1 = spm_mesh_select(mesh, 3); % will need to be in the same order that the opm positions are given (in this caseR-side of the neck, L-front of the neck, top of the head)
opm_mesh_fids = fids_select_1'; 

% Procrustes alignment (scaling, reflection, translation)
[d, Z, transform] = procrustes(m_opm_pos, opm_mesh_fids, 'Scaling', true, 'Reflection', true);

% Apply initial Procrustes transformation
initial_transformed_points = mesh.vertices * transform.T + transform.c(1, :);
initial_transformed_fids = opm_mesh_fids * transform.T + transform.c(1, :);

% Implement ICP - Fiducials as Point Clouds
pc_fids_1 = pointCloud(initial_transformed_fids);
pc_matched_fids_1 = pointCloud(m_opm_pos);

% Apply ICP to further refine the transformation
[tform_fiducials, rmse] = pcregrigid(pc_fids_1, pc_matched_fids_1); % RMSE gives alignment error

% Apply the final ICP transformation to the fiducials and the full mesh
transformed_fids1 = transformPointsForward(tform_fiducials, pc_fids_1.Location);
transformed_vertices1 = transformPointsForward(tform_fiducials, pointCloud(initial_transformed_points).Location);

opm_optical = struct();
opm_optical.vertices = transformed_vertices1;
opm_optical.faces = mesh.faces; % Faces remain unchanged
opm_optical.unit = 'm'; 
% Now the optical scan fiducials should be aligned with the OPM sensor positions
%%
% now we need to register the simulation meshes to the optical scan which
% is in opm sensor space

sim_fids_select = spm_mesh_select(opm_optical);
sim_fids = sim_fids_select';

S = [];
S.subject = opm_optical; 
S.fiducials = sim_fids;
S.plot = 'true';

T = cr_register_torso(S);

S = []; 
S.subject = opm_optical;
S.T = T;
S.sensors = sensors(D, 'MEG'); 
all_meshes = cr_check_registration(S);

torso = all_meshes.mesh_torso;
spine = all_meshes.mesh_spine;
bone = all_meshes.mesh_vertebrae;
heart = all_meshes.mesh_heart;
lungs = all_meshes.mesh_lungs;