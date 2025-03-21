function M = cr_register_thorax(S)
% Generates the affine transform required to fit the canonical thorax to a
% subject scan

if ~isfield(S,'fiducials'); error('you must provide fiducial locations'); end
if ~isfield(S,'subject'); error('you must provide a mesh to fit to!'); end
if ~isfield(S,'dist'); S.dist = 0.06; end
if ~isfield(S,'plot'); S.plot = false; end



% Load the canonical torso, generate fiducials
%----------------------------------------------
thorax = export(gifti(fullfile(coreg_path,'thorax.gii')),'patch'); % Units: m

% fiducuals of the torso
% - left shoulder   (pt 230)
% - right shoulder  (pt 1)
% - low spine       (pt 409)
thorax_fids = [thorax.vertices(230,:);
    thorax.vertices(1,:);
    thorax.vertices(409,:)];

% Register part 0: get the two meshes into similar units
%--------------------------------------------------------
sf = determine_body_scan_units(S.fiducials,thorax_fids);
M0 = [sf 0 0 0
    0 sf 0 0
    0 0 sf 0
    0 0 0 1];

thorax = spm_mesh_transform(thorax,M0);
% Update fiducials post scaling
thorax_fids = [thorax.vertices(2,:);
    thorax.vertices(1,:);
    thorax.vertices(409,:)];

if S.plot
    figure;clf
    patch(S.subject,'facecolor','none','edgecolor','k',...
        'edgealpha',0.3,'clipping','off');
    axis equal
    axis off
    set(gcf,'color','w');
    hold on
    h_fid   = plot3(S.fiducials(:,1),S.fiducials(:,2),...
        S.fiducials(:,3),'oc');
    set(h_fid,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');
end

% Register part 1: rigid body transform between fiducials
%--------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(S.fiducials',thorax_fids');
thorax = spm_mesh_transform(thorax,M1);
% Update fiducials post RBT
thorax_fids = [thorax.vertices(2,:);
    thorax.vertices(1,:);
    thorax.vertices(409,:)];

% Register part 2: ICP
%--------------------------------------------------------
[~,D] = knnsearch(S.subject.vertices,thorax.vertices);
id = find(D<=S.dist*sf);
thorax_fids = [thorax.vertices(2,:);
    thorax.vertices(1,:);
    thorax.vertices(409,:)];

M2 = spm_eeg_inv_icp(S.subject.vertices',thorax.vertices(id,:)',...
    S.fiducials',thorax_fids',[],[],1);

if S.plot
    thorax = spm_mesh_transform(thorax, M2);
    thorax_fids = [thorax.vertices(2,:);
                    thorax.vertices(1,:);
                    thorax.vertices(409,:)];
    patch(thorax,'facecolor','none','edgecolor','b','edgealpha',0.3)
    h_fidmr = plot3(thorax_fids(:,1),thorax_fids(:,2),...
        thorax_fids(:,3),'dm');
    set(h_fidmr,'MarkerFaceColor','m',...
        'MarkerSize',12,'MarkerEdgeColor','k');
end

M = M2*M1*M0;

end

function sf = determine_body_scan_units(body_fids,thorax_fids)
% Determine if the units of the thorax and body are the same by looking at
% the triangle which is made between the fiducials

body_vec = (body_fids([1 2],:) - body_fids(3,:));
thorax_vec = thorax_fids([1 2],:) - thorax_fids(3,:);
    
body_area = norm(cross(body_vec(1,:),body_vec(2,:)));
thorax_area = norm(cross(thorax_vec(1,:),thorax_vec(2,:)));

pow = round(log10(sqrt(body_area/thorax_area)));

sf = 10^pow;

end

