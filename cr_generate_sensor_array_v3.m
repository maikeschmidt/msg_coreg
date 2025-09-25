function grad = cr_generate_sensor_array_v3(S)
% CR_GENERATE_SENSOR_ARRAY_V3
% Generate torso sensors with automatic front/back axis detection
% and optional percentage coverage along head-to-lumbar and left-right axes
%
% Inputs:
%   S.torso        : torso mesh struct with .vertices (Nx3) and .faces (Mx3)
%   S.resolution   : spacing of sensors (mm, default 10)
%   S.depth        : inward shift from surface (mm, default 10)
%   S.frontflag    : [] (both), 1 (front), 0 (back)
%   S.triaxial     : 1 (default) or 0
%   S.margin       : margin to start rays outside mesh (default 50)
%   S.debug        : true to plot mesh + hits (default false)
%   S.coverage     : scalar (0–1, default 0.6) or
%                    [top bottom left right] fractions (0–1 each)

if ~isfield(S,'torso'); error('Please provide torso mesh!'); end
if ~isfield(S,'resolution'); S.resolution = 10; end
if ~isfield(S,'depth');      S.depth = 10; end
if ~isfield(S,'frontflag');  S.frontflag = []; end
if ~isfield(S,'triaxial');   S.triaxial = 1; end
if ~isfield(S,'margin');     S.margin = 50; end
if ~isfield(S,'debug');      S.debug = false; end
if ~isfield(S,'coverage');   S.coverage = 0.6; end

torso = S.torso;
axisNames = {'x-axis','y-axis','z-axis'};

% Get normals
G = gifti(torso);
[~, nrms] = spm_mesh_normals(G);

% flip normals if pointing inward
cent = mean(G.vertices,1);
if dot(nrms(1,:), G.vertices(1,:) - cent) < 0
    nrms = -nrms;
end

% ---- automatic front/back axis detection ----
ranges = max(torso.vertices) - min(torso.vertices);  % 1x3 vector
[~, fb_axis] = min(ranges); % shortest axis = front/back
other_axes = setdiff(1:3, fb_axis);
fprintf('Automatic front/back axis: %s (shortest axis)\n', axisNames{fb_axis});

% ---- head-to-lumbar axis detection ----
[~, ht_axis] = max(ranges);  % longest axis = vertical
fprintf('Automatic head–lumbar axis: %s (longest axis)\n', axisNames{ht_axis});
ht_min = min(torso.vertices(:,ht_axis));
ht_max = max(torso.vertices(:,ht_axis));

% ---- coverage option ----
if numel(S.coverage) == 1
    coverage = min(max(S.coverage,0),1);
    cut_amount = (1 - coverage) * (ht_max - ht_min);
    ht_min = ht_min + cut_amount;
    fprintf('Applying %.0f%% coverage along %s\n', coverage*100, axisNames{ht_axis});
elseif numel(S.coverage) == 4
    cov = min(max(S.coverage,0),1);

    % top/bottom (head-foot)
    ht_len = ht_max - ht_min;
    ht_min = ht_min + cov(2) * ht_len; % bottom cut
    ht_max = ht_max - cov(1) * ht_len; % top cut
    fprintf('Coverage along %s: top %.0f%%, bottom %.0f%%\n', ...
        axisNames{ht_axis}, cov(1)*100, cov(2)*100);

    % left/right (the remaining axis)
    lr_axis = setdiff(1:3,[ht_axis fb_axis]);
    lr_min = min(torso.vertices(:,lr_axis));
    lr_max = max(torso.vertices(:,lr_axis));
    lr_len = lr_max - lr_min;
    lr_min = lr_min + cov(3) * lr_len; % left cut
    lr_max = lr_max - cov(4) * lr_len; % right cut
    fprintf('Coverage along %s: left %.0f%%, right %.0f%%\n', ...
        axisNames{lr_axis}, cov(3)*100, cov(4)*100);
else
    error('S.coverage must be either 1 value or 4 values.');
end

% ---- bounding box on plane perpendicular to front/back ----
range1 = max(torso.vertices(:,other_axes(1))) - min(torso.vertices(:,other_axes(1)));
range2 = max(torso.vertices(:,other_axes(2))) - min(torso.vertices(:,other_axes(2)));

tol1 = 0.1 * range1 / 2; 
tol2 = 0.1 * range2 / 2;

min1 = min(torso.vertices(:,other_axes(1))) + tol1;
max1 = max(torso.vertices(:,other_axes(1))) - tol1;
min2 = min(torso.vertices(:,other_axes(2))) + tol2;
max2 = max(torso.vertices(:,other_axes(2))) - tol2;

[grid1, grid2] = meshgrid(min1:S.resolution:max1, min2:S.resolution:max2);
list1 = grid1(:);
list2 = grid2(:);

hit_points = zeros(0,3);
hit_normals = zeros(0,3);
hit_side = [];

% frontflag automatic if empty
if isempty(S.frontflag)
    frontflag_auto = 1; % front = min coord
else
    frontflag_auto = S.frontflag;
end

% ---- front/back passes ----
for which = 1:2
    if which==1 % front
        if ~isempty(S.frontflag) && S.frontflag==0, continue; end
        if isempty(S.frontflag) && frontflag_auto==0, continue; end
        ray = zeros(1,3); ray(fb_axis) = 1;
        start_val = min(torso.vertices(:,fb_axis)) - S.margin;
    else % back
        if ~isempty(S.frontflag) && S.frontflag==1, continue; end
        if isempty(S.frontflag) && frontflag_auto==1, continue; end
        ray = zeros(1,3); ray(fb_axis) = -1;
        start_val = max(torso.vertices(:,fb_axis)) + S.margin;
    end

    for ii = 1:numel(list1)
        origin = zeros(1,3);
        origin(other_axes(1)) = list1(ii);
        origin(other_axes(2)) = list2(ii);
        origin(fb_axis) = start_val;

        R = struct('orig', origin', 'vec', ray');
        [I,P] = spm_mesh_ray_intersect(G,R);
        if ~isempty(P)
            hits = find(I);
            [~,pid] = min(abs(P(:,fb_axis)-start_val));
            nrm = nrms(hits(pid),:); nrm=nrm/norm(nrm);

            if dot(nrm, ray) < 0
                pt = P(pid,:);
                sensor_pt = pt + S.depth*ray;

                % filter by ht_axis (coverage) and optionally lr_axis
                pass_ht = (sensor_pt(ht_axis) >= ht_min && sensor_pt(ht_axis) <= ht_max);
                if numel(S.coverage) == 4
                    pass_lr = (sensor_pt(lr_axis) >= lr_min && sensor_pt(lr_axis) <= lr_max);
                else
                    pass_lr = true;
                end

                if pass_ht && pass_lr
                    hit_points(end+1,:) = sensor_pt;
                    hit_normals(end+1,:) = nrm;
                    hit_side(end+1,1) = which;
                end
            end
        end
    end
end

nSensors = size(hit_points,1);
if nSensors==0
    warning('No sensors found. Check mesh and resolution.');
    grad = [];
    return;
end

% ---- build FT-style grad struct ----
if S.triaxial
    grad.coilpos = repmat(hit_points,3,1);
    grad.coilori = zeros(3*nSensors,3);
    grad.label = cell(3*nSensors,1);
    for i = 1:nSensors
        nrm = hit_normals(i,:)';
        T = null(nrm'); tang1=T(:,1)'; tang2=T(:,2)';
        idx = (i-1)*3 + (1:3);
        grad.coilori(idx(1),:) = nrm';
        grad.coilori(idx(2),:) = tang1;
        grad.coilori(idx(3),:) = tang2;
        grad.label{idx(1)} = sprintf('mag-%04d-R',i);
        grad.label{idx(2)} = sprintf('mag-%04d-T1',i);
        grad.label{idx(3)} = sprintf('mag-%04d-T2',i);
    end
else
    grad.coilpos = hit_points;
    grad.coilori = hit_normals;
    grad.label = cell(nSensors,1);
    for i = 1:nSensors
        grad.label{i} = sprintf('mag-%04d-R',i);
    end
end

grad.tra = speye(numel(grad.label));
[grad.chanunit{1:numel(grad.label)}] = deal('T');
[grad.chantype{1:numel(grad.label)}] = deal('megmag');
grad.unit = 'mm';
grad = ft_datatype_sens(grad,'amplitude','T','distance',grad.unit);

% ---- debug plot ----
if S.debug
    figure('Name','Torso sensors debug','Color','w'); hold on;
    patch('Faces',torso.faces,'Vertices',torso.vertices,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.35,'EdgeColor','none');
    scatter3(hit_points(:,1),hit_points(:,2),hit_points(:,3),36,hit_side,'filled');
    quiver3(hit_points(:,1),hit_points(:,2),hit_points(:,3),hit_normals(:,1),hit_normals(:,2),hit_normals(:,3),10,'k');
    axis equal; camlight headlight; lighting gouraud;
    xlabel('X'); ylabel('Y'); zlabel('Z'); title('Sensors');
end

% ---- reporting ----
fprintf('COMPLETE! %d sensors generated (triaxial=%d)\n', nSensors, S.triaxial);
nFront = sum(hit_side==1);
nBack  = sum(hit_side==2);
fprintf('Front sensors: %d, Back sensors: %d\n', nFront, nBack);

% ---- warnings for extreme coverage settings ----
if (isempty(S.frontflag) || S.frontflag==1) && nFront == 0
    warning('No FRONT sensors generated – check S.coverage values (top/bottom/left/right).');
end
if (isempty(S.frontflag) || S.frontflag==0) && nBack == 0
    warning('No BACK sensors generated – check S.coverage values (top/bottom/left/right).');
end


end

