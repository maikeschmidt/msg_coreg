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
% NOTE: 'meshes' and the filename are separate fullfile arguments so the
% separator is chosen by the platform. Writing 'meshes\canonical_torso.stl'
% as one string works on Windows but fails on macOS and Linux, where the
% backslash is a literal filename character rather than a separator.
torso_file = fullfile(coreg_path, 'meshes', 'canonical_torso.stl'); % Units: m

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


% Fiducials of the canonical torso: left shoulder, right shoulder, chin.
%
% THESE ARE COORDINATES, NOT VERTEX INDICES — AND THAT IS DELIBERATE.
%
% This code previously used the hardcoded vertex indices 3104, 8807 and 858.
% Those index the UNMERGED vertex list of canonical_torso.stl, which has
% 3 vertices per face = 12,390 entries. MATLAB's built-in stlread returns a
% triangulation with duplicate points MERGED — 2,067 vertices — so the
% indices ran off the end and the function failed with
%   "Index in position 1 exceeds array bounds. Index must not exceed 2067."
% The original author must have used an stlread that did not merge points.
%
% The coordinates below were recovered from the unmerged list at exactly
% those three indices, so they are the original intended landmarks, not a
% re-pick. They are checked as anatomically consistent: the shoulders sit
% on opposite sides in X at near-equal height, and the chin is near the
% midline, above the shoulders, and at the anterior extreme.
%
% Storing coordinates and resolving them to the nearest vertex at run time
% makes this robust to vertex merging, vertex reordering, and any change of
% STL reader — none of which preserve indices, all of which preserve
% geometry.
%
% DO NOT REPLACE THIS WITH INDICES FOR THE MERGED LIST. Merging does not
% define an ordering, so the index depends on the reader: for these same
% three points MATLAB's stlread yields 574, 1544 and 178, while a
% sorted-unique merge (e.g. numpy) yields 2031, 16 and 995. Only the
% coordinates are stable.

torso_fid_ref = [ ...
     2.0898,   0.2269,  -0.4374;   % Left shoulder
    -1.8253,   0.0529,  -0.2542;   % Right shoulder
     0.1836,   1.1304,   1.2111];  % Chin

% Resolve to vertex indices ONCE, on the untransformed mesh. The indices
% are then reused after each transform, since a transform moves vertices
% but never reorders them.
fid_idx = zeros(3,1);
fid_err = zeros(3,1);
for f = 1:3
    d = vecnorm(torso.vertices - torso_fid_ref(f,:), 2, 2);
    [fid_err(f), fid_idx(f)] = min(d);
end

% Guard against a swapped-in mesh. The reference points are vertices of the
% shipped canonical_torso.stl, quoted to 4 decimal places against float32
% STL storage, so the residual is ~5e-5 on a mesh of extent ~8 (relative
% ~9e-6). Anything materially larger means the mesh has changed and the
% landmarks are no longer the intended ones.
mesh_scale = max(max(torso.vertices) - min(torso.vertices));
if max(fid_err) > 1e-4 * mesh_scale
    warning('cr_register_torso:fiducialMismatch', ...
        ['Canonical torso fiducials did not match a mesh vertex exactly ' ...
         '(worst offset %.3g, mesh extent %.3g). canonical_torso.stl may ' ...
         'have been replaced. Re-derive torso_fid_ref before trusting the ' ...
         'registration.'], max(fid_err), mesh_scale);
end

torso_fids = torso.vertices(fid_idx, :);

% Step 1: Normalize units between subject and canonical torso
sf = determine_body_scan_units(S.fiducials, torso_fids);
M0 = diag([sf, sf, sf, 1]); % Scaling matrix

torso.vertices = (M0 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

torso_fids = torso.vertices(fid_idx, :); % Update fiducials after scaling

% Step 2: Rigid body transform based on fiducial alignment
M1 = spm_eeg_inv_rigidreg(S.fiducials', torso_fids');
torso.vertices = (M1 * [torso.vertices, ones(size(torso.vertices,1),1)]')';
torso.vertices = torso.vertices(:,1:3); % Remove homogenous coordinates

torso_fids = torso.vertices(fid_idx, :); % Update fiducials

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
