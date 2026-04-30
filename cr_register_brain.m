% cr_register_brain - Register SPM template brain meshes to a subject head surface
%
% Fits the SPM canonical brain, scalp, inner skull, and outer skull meshes
% to a subject head surface using a three-stage pipeline: unit normalisation,
% fiducial-based rigid-body alignment, and ICP refinement. Includes automatic
% detection and correction of mirror-flip ambiguity.
%
% USAGE:
%   [temp_brain, T] = cr_register_brain(S)
%
% INPUT:
%   S              - Structure with the following fields:
%
%   Required:
%     S.fiducials    - 3x3 matrix of subject fiducial coordinates [mm]:
%                        Row 1: Nasion (NAS)
%                        Row 2: Left pre-auricular point (LPA)
%                        Row 3: Right pre-auricular point (RPA)
%     S.head         - Subject head/torso mesh struct (.vertices, .faces)
%                      used to verify brain placement after registration
%
%   Optional:
%     S.dist         - ICP distance threshold in m (default: 0.02)
%     S.plot         - Logical; display registration figure (default: false)
%
% OUTPUT:
%   temp_brain     - Struct containing registered SPM canonical meshes:
%                      .brain   - Cortical surface mesh
%                      .scalp   - Scalp surface mesh
%                      .iskull  - Inner skull mesh
%                      .oskull  - Outer skull mesh
%                    Each mesh is a struct with .vertices [N x 3] and
%                    .faces [M x 3]
%   T              - 4x4 cumulative transform matrix (template → subject)
%
% DEPENDENCIES:
%   - spm()                    : SPM installation for canonical meshes and
%                                fiducial template
%   - spm_eeg_fixpnt()         : reads fiducial point file
%   - ft_read_headshape()      : FieldTrip headshape reader
%   - determine_body_scan_units() : estimates scale factor between spaces
%   - spm_eeg_inv_rigidreg()   : fiducial-based rigid-body registration
%   - spm_eeg_inv_icp()        : iterative closest point refinement
%   - tt_is_inside()           : point-in-mesh test for flip detection
%
% NOTES:
%   Registration pipeline:
%     1. Unit normalisation  — scale template to match subject units
%     2. Rigid-body alignment — align template fiducials to subject fiducials
%     3. ICP refinement      — refine using Gaussian point clouds (~20mm
%                              radius, 300 points) seeded around fiducials
%     4. Mirror-flip check   — tests 20% of brain vertices against S.head;
%                              if misplaced, reflects across the fiducial
%                              plane and updates T accordingly
%
% EXAMPLE:
%   S.fiducials = [nas; lpa; rpa];   % 3x3 matrix of coordinates
%   S.head      = torso_mesh;
%   S.plot      = true;
%   [temp_brain, T] = cr_register_brain(S);
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

function [temp_brain, T] = cr_register_brain(S)

if ~isfield(S, 'fiducials'), error('You must provide fiducial locations (nas, LP and RP)'); end
if ~isfield(S, 'head'), error('You must provide a head shape to fit to!'); end
if ~isfield(S, 'dist'), S.dist = 0.02; end
if ~isfield(S, 'plot'), S.plot = false; end

% Load the spm template brain and associated fiducials

brain_gii = gifti(fullfile(spm('dir'),'canonical','cortex_5124.surf.gii'));
scalp_gii = gifti(fullfile(spm('dir'), 'canonical','scalp_2562.surf.gii'));
iskull_gii = gifti(fullfile(spm('dir'), 'canonical','iskull_2562.surf.gii'));
oskull_gii = gifti(fullfile(spm('dir'), 'canonical','oskull_2562.surf.gii'));

brain.vertices = double(brain_gii.vertices);
brain.faces    = double(brain_gii.faces);
scalp.vertices = double(scalp_gii.vertices);
scalp.faces = double(scalp_gii.faces);
iskull.vertices = double(iskull_gii.vertices);
iskull.faces = double(iskull_gii.faces);
oskull.vertices = double(oskull_gii.vertices);
oskull.faces = double(oskull_gii.faces);

fid_template = spm_eeg_fixpnt(ft_read_headshape( ...
                fullfile(spm('dir'),'EEGtemplates','fiducials.sfp')));
fid_template = fid_template.fid.pnt(1:3,:);

% Normalise units between head and template brain

sf = determine_body_scan_units(S.fiducials, fid_template);
M0 = diag([sf,sf,sf,1]);

brain.vertices = (M0 * [brain.vertices, ones(size(brain.vertices,1),1)]')'; 
brain.vertices = brain.vertices(:,1:3); % Remove homogenous coordinates
fid_template   = (M0 * [fid_template, ones(3,1)]')';
fid_template   = fid_template(:,1:3); %update fiducials after scaling

% Rigid body transform based on fiducial alignment

M1 = spm_eeg_inv_rigidreg(S.fiducials, fid_template);

brain.vertices = (M1 * [brain.vertices, ones(size(brain.vertices,1),1)]')';
brain.vertices = brain.vertices(:,1:3); % Remove homogenous coordinates
fid_template   = (M1 * [fid_template, ones(3,1)]')';
fid_template   = fid_template(:,1:3); %update fiducials after scaling

% ICP Refinement
% need to create point clouds around the fiducials due to the brain mesh
% and head mesh being such different shapes - tehn run ICP on the clouds
% formed
cloud_template = make_fiducial_cloud(fid_template);
cloud_subject  = make_fiducial_cloud(S.fiducials);

M2 = spm_eeg_inv_icp(cloud_subject', cloud_template', ...
                         S.fiducials', fid_template', [], [], 1);

brain.vertices = (M2 * [brain.vertices, ones(size(brain.vertices,1),1)]')';
brain.vertices = brain.vertices(:,1:3); % Remove homogenous coordinates
fid_template   = (M2 * [fid_template, ones(3,1)]')';
fid_template   = fid_template(:,1:3); %update fiducials after scaling

%final transform and meshes
T = M2 * M1 * M0;
temp_brain = brain;

% check if brain is inside the mesh - mirror ambiguity issue
nverts = size(brain.vertices,1);
sample_idx = randperm(nverts, round(0.2 * nverts));
sample_points = brain.vertices(sample_idx, :);

inside_flags = false(size(sample_points,1),1);
for i = 1:size(sample_points,1)
    inside_flags(i) = tt_is_inside(sample_points(i,:), S.head.vertices, S.head.faces);
end

if ~all(inside_flags)
    warning('Brain appears flipped — applying mirroring across fiducial plane');

    % Define plane from 3 fiducials: Nasion (p0), LPA (p1), RPA (p2)
    p0 = S.fiducials(1,:); % NAS
    p1 = S.fiducials(2,:); % LPA
    p2 = S.fiducials(3,:); % RPA

    % Plane normal: n = (p1 - p0) x (p2 - p0)
    n = cross(p1 - p0, p2 - p0);
    n = n / norm(n);

    % Construct mirror (reflection) matrix
    R = eye(3) - 2 * (n' * n);  

    % Apply mirroring about the plane
    verts = brain.vertices;
    brain.vertices = (R * (verts - p0)' )' + p0;

    % Update final transform
    T_mirror = eye(4);
    T_mirror(1:3,1:3) = R;
    T_mirror(1:3,4) = p0' - R*p0';
    T = T_mirror * T;
end

scalp.vertices = (T * [scalp.vertices, ones(size(scalp.vertices,1),1)]')';
scalp.vertices = scalp.vertices(:,1:3); % Remove homogenous coordinates
iskull.vertices = (T * [iskull.vertices, ones(size(iskull.vertices,1),1)]')';
iskull.vertices = iskull.vertices(:,1:3); % Remove homogenous coordinates
oskull.vertices = (T * [oskull.vertices, ones(size(oskull.vertices,1),1)]')';
oskull.vertices = oskull.vertices(:,1:3); % Remove homogenous coordinates

temp_brain.brain = brain;
temp_brain.scalp = scalp;
temp_brain.iskull = iskull;
temp_brain.oskull = oskull;

% Visualisation
if S.plot
    figure; clf;
    hold on;
    patch('Vertices', S.head.vertices, 'Faces', S.head.faces,...
        'FaceColor', 'none', 'EdgeColor', 'k', 'EdgeAlpha', 0.3);
    patch('Vertices', brain.vertices, 'Faces', brain.faces,...
        'FaceColor', 'none', 'Edgecolor', 'b', 'EdgeAlpha', 0.3);
    patch('Vertices', iskull.vertices, 'Faces', iskull.faces,...
    'FaceColor', 'none', 'Edgecolor', 'r', 'EdgeAlpha', 0.3);
    patch('Vertices', oskull.vertices, 'Faces', oskull.faces,...
    'FaceColor', 'none', 'Edgecolor', 'g', 'EdgeAlpha', 0.3);

    plot3(S.fiducials(:,1), S.fiducials(:,2), S.fiducials(:,3), 'co', ...
          'MarkerFaceColor', 'c', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    plot3(fid_template(:,1), fid_template(:,2), fid_template(:,3), 'mo', ...
          'MarkerFaceColor', 'm', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');
    
    axis equal;
    axis off;
    set(gcf, 'Color', 'w');
    hold off;
end

end


function cloud = make_fiducial_cloud(fid_pts)
    % random Gaussian cloud around fiducials 
    n = 300; % number of samples per fid
    cloud = [];
    for i = 1:size(fid_pts,1)
        r = randn(n,3);              
        r = r ./ vecnorm(r,2,2);     
        rad = 20 + randn(n,1)*0.5;    % radius 20mm +/- 0.5
        pts = fid_pts(i,:) + r .* rad;
        cloud = [cloud; pts];
    end
end

