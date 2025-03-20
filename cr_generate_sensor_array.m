function grad = cr_generate_sensor_array(S)

if ~isfield(S,'subject'); error('please provide subject mesh!'); end
if ~isfield(S,'T'); error('please provide the transformation matrix!'); end
if ~isfield(S,'resolution');    S.resolution = 10;      end
if ~isfield(S,'depth');         S.depth = 10;      end
if ~isfield(S,'frontflag');         S.frontflag = 0;      end
if ~isfield(S,'zlim');         S.zlim = [];      end
if ~isfield(S,'triaxial');     S.triaxial = 1; end
if ~isfield(S,'unit');  S.unit = 'm'; end

%-rotate the body scan for easier grid generation later
%-----------------------------------------------------
fids = S.fids;
torso_path = fullfile(coreg_path, 'torso.stl');
torso = [];
temp = stlread(torso_path);
torso.vertices = temp.vertices; %Points;
torso.faces = temp.faces; %ConnectivityList;

% halfway point between shoulder fids
hp = 0.5*(fids(2,:) + fids(1,:));
% vector between shoulder fids
hpv = fids(2,:) - fids(1,:);

% find the rotation in z which aligns the shoulders in y (is this a special
% case with just this scan?)
angD = -atand(hpv(2)/hpv(1));
R1 = rotmatZ(angD,hp);

% rotate the meshes for the grid generation
torso_rot = spm_mesh_transform(torso,R1);
sub_rot = spm_mesh_transform(S.subject,R1);
fids = spm_mesh_transform(S.fids,R1*S.T);

% work out if the spine is in the +/-ve direction along y.
tmp = fids(3,:) - hp;
ydir = heaviside(tmp(2));

if S.frontflag,
    ydir=~ydir;
end;

% get units to scale space with later;
unit = S.unit;

%-Generate a plane of the back and raycast to get shape
%--------------------------------------------------------------------
min_x = min(torso_rot.vertices(:,1));
max_x = max(torso_rot.vertices(:,1));
min_z = min(torso_rot.vertices(:,3));
max_z = max(torso_rot.vertices(:,3));

if ~isempty(S.zlim),
    min_z=max(min_z,S.zlim(1));
    max_z=min(max_z,S.zlim(2));
end;
if S.frontflag,
    y_start = hp(2) - 10*tmp(2);
else
    y_start = hp(2) + 10*tmp(2);
end;

[xgrid,ygrid,zgrid] = meshgrid(min_x:S.resolution:max_x...
    ,y_start,min_z:S.resolution:max_z);
t1=torso_rot.vertices;
% figure;
% plot3(t1(:,1),t1(:,2),t1(:,3),'r.');axis equal; hold on;
% plot3(t1(:,1),ones(size(t1(:,2))).*y_start,t1(:,3),'b.');axis equal



switch ydir
    case 1
        ray = [0 -1 0];
    case 0
        ray = [0 1 0];
end


% First pass based on the torso
G = gifti(torso_rot);
[~,nrms] = spm_mesh_normals(G);

fprintf('Generating sensor positions: Pass 1/2\n')

plane_x1 = [];
plane_y1 = [];
plane_z1 = [];
grid_id = [];
for ii = 1:numel(xgrid)
    R = struct('orig',[xgrid(ii) ygrid(ii) zgrid(ii)]',...
        'vec',ray');
    [I,P] = spm_mesh_ray_intersect(G,R);
    if ~isempty(P)
        hits = find(I);
        if size(P,1) > 1
            switch ydir
                case 1
                    [~,pid] = max(P(:,2));
                case 0
                    [~,pid] = min(P(:,2));
            end
        else
            pid = 1;
        end

        % filter by angle, to not get side-on sensors
        nrm = nrms(hits(pid),:);
        nrm = nrm./norm(nrm);
        ang = abs(acosd(dot(nrm,ray)));
        ang = min(ang,180-abs(ang));

        if ang < 90 %% keep more

            plane_x1(end+1) = P(pid,1);
            plane_y1(end+1) = P(pid,2);
            plane_z1(end+1) = P(pid,3);
            grid_id(end+1) = ii;
        end
    end
end

% second pass based on subject scan.
fprintf('Generating sensor positions: Pass 2/2\n')
G = gifti(sub_rot);
plane_x2 = [];
plane_y2 = [];
plane_z2 = [];
for ii = 1:numel(grid_id)
    R = struct('orig',[xgrid(grid_id(ii)) ygrid(grid_id(ii)) zgrid(grid_id(ii))]',...
        'vec',ray');
    [I,P] = spm_mesh_ray_intersect(G,R);
    if ~isempty(P)
        hits = find(I);
        if size(P,1) > 1
            switch ydir
                case 1
                    [~,pid] = max(P(:,2));
                case 0
                    [~,pid] = min(P(:,2));
            end
        else
            pid = 1;
        end
        plane_x2(end+1) = P(pid,1);
        plane_y2(end+1) = P(pid,2);
        plane_z2(end+1) = P(pid,3);
    else
        plane_x2(end+1) = NaN;
        plane_y2(end+1) = NaN;
        plane_z2(end+1) = NaN;
    end
end

% find nans from the second pass
id = find(isnan(plane_z2));
plane_x1(id) = [];
plane_x2(id) = [];
plane_y1(id) = [];
plane_y2(id) = [];
plane_z1(id) = [];
plane_z2(id) = [];

switch ydir %% 1 for sensors on back in example subj
    case 1
        plane_y = max([plane_y1; plane_y2]);
    case 0
        plane_y = min([plane_y1; plane_y2]);
end

fprintf('COMPLETE!\n')

%-shift across depths
%--------------------
grid_all = [plane_x1; plane_y; plane_z1]'- ray*S.depth;


%-make FT style source grad container, transform back to orig space
%-------------------------------------------------------------------
grad = [];
if S.triaxial
    grad.coilpos = repmat(grid_all,3,1);
    grad.coilori = [repmat([1 0 0],length(grid_all),1);...
        repmat([0 1 0],length(grid_all),1);...
        repmat([0 0 1],length(grid_all),1)];
    grad.label = [];
    for prefix = {'X','Y','Z'}
        for ii = 1:length(grid_all)
            grad.label{end+1} = sprintf('mag-%04d-%s',ii,prefix{:});
        end
    end
    grad.label = grad.label';
else % make it just one axis
    grad.coilpos = grid_all;
    grad.coilori = repmat([0 1 0],length(grid_all),1);
    grad.label = [];
    for prefix = {'Y'}
        for ii = 1:length(grid_all)
            grad.label{end+1} = sprintf('mag-%04d-%s',ii,prefix{:});
        end
    end
end

grad.tra = speye(numel(grad.label));
[grad.chanunit{1:numel(grad.label)}] = deal('T');
[grad.chantype{1:numel(grad.label)}] = deal('megmag');
grad.unit = unit;
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', unit);
grad = ft_transform_geometry(inv(R1),grad);

end
function M = rotmatZ(deg,cp)
% make a rotation matrix around the Z axis centered at a point of choice
MT = [1 0 0 -cp(1);
    0 1 0 -cp(2);
    0 0 1 -cp(3);
    0 0 0 1];

MR = [cosd(deg) -sind(deg) 0 0;
    sind(deg) cosd(deg) 0 0
    0 0 1 0;
    0 0 0 1];

iMT = [1 0 0 cp(1);
    0 1 0 cp(2);
    0 0 1 cp(3);
    0 0 0 1];

M = iMT*MR*MT;
end