%% Example script to co-register between OPM sensor positions, Optical Scan, and Simulation mesh models
% the only thing to make sure that is in path is spm - otherwise i *think* this should have it all

% This is early days, so there may be some snags in the pipeline - let me know!!!

% First, need to register the halo results to the optical scan

halo = []; % Placeholder for halo file import
marker_opms = []; % Identify which OPMs were used as the markers
m_opm_pos = []; % Extract the positions of these sensors

% Read and preprocess the optical scan
optical = ft_read_headshape('D:\Optical Scan\17_10_24\OP00180_scan_spineneck.stl', 'unit', 'm');
p = struct();
p.vertices = optical.pos;
p.faces = optical.tri;
p2 = reducepatch(p, 0.05); % Reduce complexity of the mesh
subject.pos = p2.vertices;
subject.tri = p2.faces;

% Prepare mesh structure
mesh = struct();
mesh.vertices = subject.pos;
mesh.faces = subject.tri;
mesh.unit = 'm'; % Define unit explicitly

% Select the fiducials where the OPMs were placed
fids_select_1 = spm_mesh_select(mesh, 3); % will need to be in the same order that the opm positions are given
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
S.sensors = ; %can define the sensors to plot here (i havnt tried this yet with real data only simulated arrays but in theory should work)
all_meshes = cr_check_registration(S);
