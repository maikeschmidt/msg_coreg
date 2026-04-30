% cr_register_torso - Register the canonical torso mesh to a subject body scan
%
% Computes the affine transform required to fit the canonical torso mesh
% (including head and neck) to a subject surface scan using a three-stage
% pipeline: unit normalisation, fiducial-based rigid-body alignment, and
% ICP refinement. The resulting transform matrix can be passed directly to
% cr_load_meshes() and cr_generate_sensor_array_v4().
%
% USAGE:
%   M = cr_register_torso(S)
%
% INPUT:
%   S              - Structure with the following fields:
%
%   Required:
%     S.fiducials    - 3x3 matrix of subject fiducial coordinates [mm]:
%                        Row 1: Left shoulder
%                        Row 2: Right shoulder
%                        Row 3: Chin
%     S.subject      - Subject body scan mesh struct (.vertices, .faces)
%
%   Optional:
%     S.dist         - ICP distance threshold in m for point selection
%                      (default: 0.02)
%     S.plot         - Logical; display registration figure (default: false)
%
% OUTPUT:
%   M              - 4x4 cumulative affine transform matrix
%                    (canonical torso → subject space)
%                    Combines scaling (M0), rigid-body (M1), and ICP (M2):
%                    M = M2 * M1 * M0
%
% DEPENDENCIES:
%   - coreg_path()               : locates canonical_torso.stl
%   - stlread()                  : reads STL files (supports modern and
%                                  legacy MATLAB formats)
%   - determine_body_scan_units(): estimates scale factor from fiducial
%                                  triangle areas (defined in this file)
%   - spm_eeg_inv_rigidreg()     : fiducial-based rigid-body registration
%   - spm_eeg_inv_icp()          : iterative closest point refinement
%   - knnsearch()                : MATLAB Statistics Toolbox k-NN search
%
% NOTES:
%   Registration pipeline:
%     1. Unit normalisation  — scale factor estimated from the area of the
%                              fiducial triangle (shoulder-shoulder-chin)
%                              in each space; rounds to nearest power of 10
%     2. Rigid-body alignment — aligns canonical torso fiducials (vertices
%                              3104, 8807, 858) to subject fiducials
%     3. ICP refinement      — refines using canonical torso vertices
%                              within S.dist * sf of the subject surface
%   - The canonical torso STL is loaded from:
%     <repo_root>/meshes/canonical_torso.stl
%   - Fiducial vertex indices (3104, 8807, 858) correspond to left
%     shoulder, right shoulder, and chin respectively
%
% EXAMPLE:
%   disp('Select fiducials: left shoulder, right shoulder, chin');
%   S.fiducials = spm_mesh_select(subject_mesh)';
%   S.subject   = subject_mesh;
%   S.plot      = true;
%   M = cr_register_torso(S);
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

function M = cr_register_torso(S)

if ~isfield(S,'fiducials'), error('You must provide fiducial locations'); end
if ~isfield(S,'subject'), error('You must provide a mesh to fit to!'); end
if ~isfield(S,'dist'), S.dist = 0.02; end
if ~isfield(S,'plot'), S.plot = false; end

% Load the canonical torso and generate fiducials
torso_file = fullfile(coreg_path,'meshes\canonical_torso.stl'); % Units: m

stl_data = stlread(torso_file);
torso = struct();

if isa(stl_data, 'triangulation')
    torso.vertices = stl_data.Points;
    torso.faces    = stl_data.ConnectivityList;

elseif isstruct(stl_data)

    if isfield(stl_data,'vertices') && isfield(stl_data,'faces')
        torso.vertices = stl_data.vertices;
        torso.faces    = stl_data.faces;

    elseif isfield(stl_data,'Vertices') && isfield(stl_data,'Faces')
        torso.vertices = stl_data.Vertices;
        torso.faces    = stl_data.Faces;

    elseif isfield(stl_data,'points') && isfield(stl_data,'ConnectivityList')
        torso.vertices = stl_data.points;
        torso.faces    = stl_data.ConnectivityList;

    elseif isfield(stl_data,'pts') && isfield(stl_data,'tri')
        torso.vertices = stl_data.pts;
        torso.faces    = stl_data.tri;

    else
        error('Unsupported STL format in %s: cannot find vertices/faces.', torso_file);
    end

else
    error('Unknown data type returned by stlread for %s.', torso_file);
end


% Fiducials of the torso
% - Left shoulder (point 3107)
% - Right shoulder (point 8838)
% - Chin (point 860)
% - Lower Spine (point 5568)

torso_fids = [torso.vertices(3104,:);
              torso.vertices(8807,:);
              torso.vertices(858,:)];
              % torso.vertices(5556,:)];
% Step 1: Normalize units between subject and canonical torso
sf = determine_body_scan_units(S.fiducials, torso_fids);
M0 = diag([sf, sf, sf, 1]); % Scaling matrix

torso.vertices = (M0 * [torso.vertices, ones(size(torso.vertices,1),1)]')'; 
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

torso_fids = [torso.vertices(3104,:);
              torso.vertices(8807,:);
              torso.vertices(858,:)];
              % torso.vertices(5556,:)]; % Update fiducials after scaling

% Step 2: Rigid body transform based on fiducial alignment
M1 = spm_eeg_inv_rigidreg(S.fiducials', torso_fids');
torso.vertices = (M1 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

torso_fids = [torso.vertices(3104,:);
              torso.vertices(8807,:);
              torso.vertices(858,:)];
              % torso.vertices(5556,:)];% Update fiducials

% Step 3: ICP Refinement
[~, D] = knnsearch(S.subject.vertices, torso.vertices);
id = find(D <= S.dist * sf); % Select closest points

M2 = spm_eeg_inv_icp(S.subject.vertices', torso.vertices(id,:)', ...
                     S.fiducials', torso_fids', [], [], 1);

torso.vertices = (M2 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

% Step 4: Visualization (Plotting)
if S.plot
    figure; clf;
    hold on;
    patch('Vertices', S.subject.vertices, 'Faces', S.subject.faces, ...
          'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', 0.3);
    patch('Vertices', torso.vertices, 'Faces', torso.faces, ...
          'FaceColor', 'none', 'EdgeColor', 'b', 'EdgeAlpha', 0.3);

    % Plot fiducials
    plot3(S.fiducials(:,1), S.fiducials(:,2), S.fiducials(:,3), 'co', ...
          'MarkerFaceColor', 'c', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    plot3(torso_fids(:,1), torso_fids(:,2), torso_fids(:,3), 'mo', ...
          'MarkerFaceColor', 'm', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    
    axis equal;
    axis off;
    set(gcf, 'Color', 'w');
    hold off;
end

% Compute final transformation matrix
M = M2 * M1 * M0;

end


function sf = determine_body_scan_units(body_fids,torso_fids)
% Determine if the units of the torso and body are the same by looking at
% the triangle which is made between the fiducials

body_vec = (body_fids([1 2],:) - body_fids(3,:));
thorax_vec = torso_fids([1 2],:) - torso_fids(3,:);
    
body_area = norm(cross(body_vec(1,:),body_vec(2,:)));
thorax_area = norm(cross(thorax_vec(1,:),thorax_vec(2,:)));

pow = round(log10(sqrt(body_area/thorax_area)));

sf = 10^pow;

end
