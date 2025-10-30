function grad = cr_generate_sensor_array_v4(S)
% Generate triaxial MEG/OPM sensor array across torso surface with coverage
%
% Required fields in S:
%   S.subject   - subject mesh (struct with .vertices and .faces) OR path to mesh
%   S.T         - 4x4 transform matrix (canonical torso → subject)
%
% Optional fields:
%   S.resolution - grid spacing in mm (default 10)
%   S.depth      - depth offset from skin in mm (default 10)
%   S.frontflag  - 0=back (default), 1=front
%   S.zlim       - [min_z, max_z] for cropping
%   S.triaxial   - 1 (default)
%   S.coverage   - scalar (0..1) OR [top bottom left right] (0..1)

if ~isfield(S,'subject'), error('Please provide subject mesh!'); end
if ~isfield(S,'T'), error('Please provide transform matrix S.T!'); end
if ~isfield(S,'resolution'), S.resolution = 10; end
if ~isfield(S,'depth'), S.depth = 10; end
if ~isfield(S,'frontflag'), S.frontflag = 0; end
if ~isfield(S,'zlim'), S.zlim = []; end
if ~isfield(S,'triaxial'), S.triaxial = 1; end
if ~isfield(S,'coverage'), S.coverage = 0.6; end

% first load the torso
fprintf('Loading canonical torso and fiducials...\n')

if ~isfield(S,'torsotype'), S.torsotype = 'anatomical'; end
switch lower(S.torsotype)
    case 'anatomical'
        torso_file = fullfile(coreg_path,'mri_torso.stl');
    case 'canonical'
        torso_file = fullfile(coreg_path,'canonical_torso.stl');
    otherwise
        error('S.torsotype must be ''anatomical'' or ''canonical''');
end

if ~isfile(torso_file)
    error('Torso mesh not found: %s', torso_file);
end

torso = stlread(torso_file);


% Load fiducials corresponding to torso type
fids = cr_get_fids(S.torsotype);  
if size(fids,1) ~=4 || size(fids,2)~=3
    error('cr_get_fids(%s) must return 4x3 [Lsh;Rsh;chin;lumbar]', S.torsotype);
end


% canonical fiducials 
fids = cr_get_fids();   % rows: [L_sh,R_sh,chin,lumbar]
if size(fids,1) < 4 || size(fids,2) ~= 3
    error('cr_get_fids must return a 4x3 matrix of fiducials [Lsh; Rsh; chin; lumbar]');
end

% transform canonical -> subject
torso = spm_mesh_transform(torso, S.T);
fids_h = (S.T * [fids, ones(size(fids,1),1)]')';
fids = fids_h(:,1:3);

% shoulder midpoint & spine direction 
hp = 0.5*(fids(1,:) + fids(2,:));    % midpoint between left/right shoulder
spine_vec = fids(4,:) - hp;          % shoulder->lumbar
spine_positive = (spine_vec(2) > 0); % lumbar is +Y (true if lumbar is larger Y)
ydir = double(spine_positive);

% compute bounds (after transform)
min_x = min(torso.vertices(:,1)); max_x = max(torso.vertices(:,1));
min_y = min(torso.vertices(:,2)); max_y = max(torso.vertices(:,2));
min_z = min(torso.vertices(:,3)); max_z = max(torso.vertices(:,3));
if ~isempty(S.zlim)
    min_z = max(min_z, S.zlim(1));
    max_z = min(max_z, S.zlim(2));
end

% choose start plane along Z and ray direction 
z_range = max_z - min_z;
offset_dist = max(80, 0.5*z_range);

if S.frontflag == 1
    % FRONT (anterior: +Z direction)
    z_start = max_z + offset_dist;
    ray = [0 0 -1];  % shoot toward -Z into torso
else
    % BACK (posterior: -Z direction)
    z_start = min_z - offset_dist;
    ray = [0 0 1];   % shoot toward +Z into torso
end

% create grid (X,Y vary, Z fixed at z_start) 
[xgrid, ygrid, zgrid] = meshgrid(min_x:S.resolution:max_x, ...
                                 min_y:S.resolution:max_y, ...
                                 z_start);

% Raycast torso (pass 1) 
fprintf('Raycasting torso surface (Pass 1/2)...\n')
Gtor = gifti; Gtor.vertices = torso.vertices; Gtor.faces = torso.faces;
[~, nrms] = spm_mesh_normals(Gtor);

plane_xyz = zeros(0,3);
grid_id = [];

for ii = 1:numel(xgrid)
    R = struct('orig',[xgrid(ii) ygrid(ii) zgrid(ii)]', 'vec', ray');
    [I,P] = spm_mesh_ray_intersect(Gtor,R);
    if ~isempty(P)
        hits = find(I);
        if size(P,1) > 1
            switch S.frontflag
                case 1, [~,pid] = min(P(:,3)); % front: smallest Z hit along ray inward
                case 0, [~,pid] = max(P(:,3)); % back: largest Z
            end
        else
            pid = 1;
        end

        nrm = nrms(hits(pid),:); nrm = nrm./norm(nrm);
        ang = abs(acosd(dot(nrm, ray))); ang = min(ang,180-abs(ang));

        if ang < 90
            plane_xyz(end+1,:) = P(pid,:);
            grid_id(end+1) = ii;
        end
    end
end

% Raycast subject (pass 2) 
fprintf('Raycasting subject surface (Pass 2/2)...\n')
% prepare subject mesh
if isstruct(S.subject) && isfield(S.subject,'vertices') && isfield(S.subject,'faces')
    Gsub = gifti; Gsub.vertices = S.subject.vertices; Gsub.faces = S.subject.faces;
else
    try
        Gsub = gifti(S.subject);
    catch
        error('S.subject must be a mesh struct with vertices/faces or a valid gifti filename');
    end
end
[~, nrms_sub] = spm_mesh_normals(Gsub);

nGrid = numel(grid_id);
plane_xyz2 = nan(nGrid,3);
plane_n2 = nan(nGrid,3);

for kk = 1:nGrid
    ii = grid_id(kk);
    R = struct('orig',[xgrid(ii) ygrid(ii) zgrid(ii)]', 'vec', ray');
    [I,P] = spm_mesh_ray_intersect(Gsub, R);
    if ~isempty(P)
        hits = find(I);
        if size(P,1) > 1
            switch S.frontflag
                case 1, [~,pid] = min(P(:,3));
                case 0, [~,pid] = max(P(:,3));
            end
        else
            pid = 1;
        end
        plane_xyz2(kk,:) = P(pid,:);
        plane_n2(kk,:)   = nrms_sub(hits(pid),:);
    end
end

% remove rows with no hit in second pass
valid = all(~isnan(plane_xyz2),2);
plane_xyz = plane_xyz(valid,:);
plane_xyz2 = plane_xyz2(valid,:);
plane_n2 = plane_n2(valid,:);

% choose final Z per point 
switch S.frontflag
    case 1, plane_z = min([plane_xyz(:,3), plane_xyz2(:,3)],[],2); % front
    case 0, plane_z = max([plane_xyz(:,3), plane_xyz2(:,3)],[],2); % back
end

% assemble final positions (before depth)
grid_pre_depth = [plane_xyz(:,1), plane_xyz(:,2), plane_z];

% apply coverage 
% head->lumbar axis is Y (ht_axis=2); left-right axis is X (lr_axis=1)
ht_axis = 2;
lr_axis = 1;
ht_min_all = min(grid_pre_depth(:,ht_axis));
ht_max_all = max(grid_pre_depth(:,ht_axis));

if numel(S.coverage) == 1
    coverage = min(max(S.coverage,0),1);
    cut_amount = (1 - coverage) * (ht_max_all - ht_min_all);
    ht_min = ht_min_all + cut_amount/2;
    ht_max = ht_max_all - cut_amount/2;
    lr_min = min(grid_pre_depth(:,lr_axis));
    lr_max = max(grid_pre_depth(:,lr_axis));
elseif numel(S.coverage) == 4
    cov = min(max(S.coverage,0),1);
    ht_len = ht_max_all - ht_min_all;
    ht_min = ht_min_all + cov(2) * ht_len; % bottom cut
    ht_max = ht_max_all - cov(1) * ht_len; % top cut
    lr_min_all = min(grid_pre_depth(:,lr_axis));
    lr_max_all = max(grid_pre_depth(:,lr_axis));
    lr_len = lr_max_all - lr_min_all;
    lr_min = lr_min_all + cov(3) * lr_len; % left cut
    lr_max = lr_max_all - cov(4) * lr_len; % right cut
else
    error('S.coverage must be scalar or 4-element vector');
end

% apply cuts
in_ht = grid_pre_depth(:,ht_axis) >= ht_min & grid_pre_depth(:,ht_axis) <= ht_max;
in_lr = grid_pre_depth(:,lr_axis) >= lr_min & grid_pre_depth(:,lr_axis) <= lr_max;
keep = in_ht & in_lr;

grid_kept = grid_pre_depth(keep,:);
normals_kept = plane_n2(keep,:);   % corresponding normals (if needed)
if isempty(grid_kept)
    warning('Coverage settings removed all sensors. Returning empty grad.');
    grad = [];
    return;
end

% apply depth shift along ray
grid_all = grid_kept - ray * S.depth;

% build grad (triaxial)
fprintf('Creating triaxial sensors...\n')
nSensors = size(grid_all,1);
grad = [];
unit = 'mm';

% Option A: keep fixed XYZ axes 
if S.triaxial
    grad.coilpos = repmat(grid_all, 3, 1);
    grad.coilori = [repmat([1 0 0], nSensors, 1);
                    repmat([0 1 0], nSensors, 1);
                    repmat([0 0 1], nSensors, 1)];
    grad.label = cell(3*nSensors,1);
    for ii = 1:nSensors
        idx = (ii-1)*3 + (1:3);
        grad.label{idx(1)} = sprintf('mag-%04d-R',ii);
        grad.label{idx(2)} = sprintf('mag-%04d-T1',ii);
        grad.label{idx(3)} = sprintf('mag-%04d-T2',ii);
    end
else
    grad.coilpos = grid_all;
    grad.coilori = repmat(ray, nSensors, 1);
    grad.label = arrayfun(@(ii) sprintf('mag-%04d-R',ii), 1:nSensors, 'uni', 0)';
end

% If you'd rather have radial = surface normal and tangents from null,
% replace the coilori block above with the code below :
%{
if S.triaxial
    grad.coilpos = repmat(grid_all, 3, 1);
    grad.coilori = zeros(3*nSensors,3);
    grad.label = cell(3*nSensors,1);
    for i = 1:nSensors
        nrm = normals_kept(i,:)';
        if norm(nrm) < eps
            nrm = [0 1 0]';
        else
            nrm = nrm / norm(nrm);
        end
        T = null(nrm'); tang1 = T(:,1)'; tang2 = T(:,2)';
        idx = (i-1)*3 + (1:3);
        grad.coilori(idx(1),:) = nrm';
        grad.coilori(idx(2),:) = tang1;
        grad.coilori(idx(3),:) = tang2;
        grad.label{idx(1)} = sprintf('mag-%04d-R',i);
        grad.label{idx(2)} = sprintf('mag-%04d-T1',i);
        grad.label{idx(3)} = sprintf('mag-%04d-T2',i);
    end
end
%}

grad.tra = speye(numel(grad.label));
[grad.chanunit{1:numel(grad.label)}] = deal('T');
[grad.chantype{1:numel(grad.label)}] = deal('megmag');
grad.unit = unit;
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', unit);

fprintf('Sensor array generation complete! (%d sensors)\n', nSensors);
end