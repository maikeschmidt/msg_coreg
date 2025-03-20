function M = cr_register_torso(S)
% Generates the affine transform required to fit the canonical torso 
% (including head and neck) to a subject scan.

if ~isfield(S,'fiducials'), error('You must provide fiducial locations'); end
if ~isfield(S,'subject'), error('You must provide a mesh to fit to!'); end
if ~isfield(S,'dist'), S.dist = 0.06; end
if ~isfield(S,'plot'), S.plot = false; end

% Load the canonical torso and generate fiducials
%----------------------------------------------
torso_file = fullfile(coreg_path,'torso.stl'); % Units: m
torso_temp = stlread(torso_file);
torso = [];
torso.vertices = torso_temp.vertices;
torso.faces = torso_temp.faces;

% Fiducials of the torso
% - Left shoulder (point 3093)
% - Right shoulder (point 8774)
% - Chin (point 887)
torso_fids = [torso.vertices(3093,:);
              torso.vertices(8774,:);
              torso.vertices(887,:)];

% Step 1: Normalize units between subject and canonical torso
%--------------------------------------------------------
sf = determine_body_scan_units(S.fiducials, torso_fids);
M0 = diag([sf, sf, sf, 1]); % Scaling matrix

torso.vertices = (M0 * [torso.vertices, ones(size(torso.vertices,1),1)]')'; 
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates
torso_fids = [torso.vertices(3093,:);
              torso.vertices(8774,:);
              torso.vertices(887,:)]; % Update fiducials after scaling

% Step 2: Rigid body transform based on fiducial alignment
%--------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(S.fiducials', torso_fids');
torso.vertices = (M1 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates
torso_fids = [torso.vertices(3093,:);
              torso.vertices(8774,:);
              torso.vertices(887,:)]; % Update fiducials

% Step 3: ICP Refinement
%--------------------------------------------------------
[~, D] = knnsearch(S.subject.vertices, torso.vertices);
id = find(D <= S.dist * sf); % Select closest points

M2 = spm_eeg_inv_icp(S.subject.vertices', torso.vertices(id,:)', ...
                     S.fiducials', torso_fids', [], [], 1);

torso.vertices = (M2 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

% Step 4: Visualization (Plotting)
%--------------------------------------------------------
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
% Determine if the units of the thorax and body are the same by looking at
% the triangle which is made between the fiducials

body_vec = (body_fids([1 2],:) - body_fids(3,:));
thorax_vec = torso_fids([1 2],:) - torso_fids(3,:);
    
body_area = norm(cross(body_vec(1,:),body_vec(2,:)));
thorax_area = norm(cross(thorax_vec(1,:),thorax_vec(2,:)));

pow = round(log10(sqrt(body_area/thorax_area)));

sf = 10^pow;

end

