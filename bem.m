
clearvars
close all
clc

cd('D:\');
Metadata;
proj_init;
cd('D:\OPM\sub-OP00210\fwd');

%% Load geometries
geoms_path = ('D:\OPM\sub-OP00210\fwd' );
geoms = load(fullfile(geoms_path, 'geometries_realistic_bone_cervical.mat'));

ordering_cord = {'wm','bone','heart','lungs','torso'};
ordering_brain = {'brain', 'oskull', 'torso'};

tt_add_bem;
reduction_factor = 0.2;  
reduction_torso = 0.5;

for ii = 1:numel(ordering_cord)
    field = ['mesh_' ordering_cord{ii}];
    mesh_tmp = geoms.(field);

    if ii <= 2
        patch_in.vertices = mesh_tmp.vertices;
        patch_in.faces = mesh_tmp.faces;

        patch_out = reducepatch(patch_in, reduction_factor);

        pos = patch_out.vertices;
        tri = patch_out.faces;
    elseif ii == 5
        patch_in.vertices = mesh_tmp.vertices;
        patch_in.faces = mesh_tmp.faces;

        patch_out = reducepatch(patch_in, reduction_torso);

        pos = patch_out.vertices;
        tri = patch_out.faces;
    else
        pos = mesh_tmp.vertices;
        tri = mesh_tmp.faces;
    end

    bnd_cord(ii).pos = pos;
    bnd_cord(ii).tri = tri;
    bnd_cord(ii).unit = 'mm';

    orient = hbf_CheckTriangleOrientation(bnd_cord(ii).pos, bnd_cord(ii).tri);
    if orient == 2
        bnd_cord(ii).tri = bnd_cord(ii).tri(:, [1 3 2]);
    end

    bnd_cord(ii) = ft_convert_units(bnd_cord(ii), 'm');
end

for ii = 1:numel(ordering_brain)
    field = ['mesh_' ordering_brain{ii}];
    mesh_tmp = geoms.(field);

    patch_in.vertices = mesh_tmp.vertices;
    patch_in.faces = mesh_tmp.faces;

    patch_out = reducepatch(patch_in, reduction_torso);

    pos = patch_out.vertices;
    tri = patch_out.faces;

    bnd_brain(ii).pos = pos;
    bnd_brain(ii).tri = tri;
    bnd_brain(ii).unit = 'mm';

    orient = hbf_CheckTriangleOrientation(bnd_brain(ii).pos, bnd_brain(ii).tri);
    if orient == 2
        bnd_brain(ii).tri = bnd_brain(ii).tri(:, [1 3 2]);
    end

    bnd_brain(ii) = ft_convert_units(bnd_brain(ii), 'm');
end
%% Conductivities
cratio = 40;
ci_cord = [0.33 (0.33/cratio) .62 .05 .23];
co_cord = [(0.33/cratio) .23 .23 .23 0];

ci_brain = [0.33 0.0042 0.23];
co_brain = [0.0042 0.33 0];

cfg = [];
cfg.method = 'bem_hbf';
cfg.conductivity = [ci_cord;co_cord];
vol_cord = ft_prepare_headmodel(cfg, bnd_cord);

cfg = [];
cfg.method = 'bem_hbf';
cfg.conductivity = [ci_brain; co_brain];
vol_brain = ft_prepare_headmodel(cfg, bnd_brain);

%% Spinal cord sources
sources_spine = [];
sources_spine.pos = geoms.sources_cent.pos;
sources_spine.inside = true(size(sources_spine.pos, 1), 1);
sources_spine.unit = 'mm';
sources_spine = ft_convert_units(sources_spine, 'm');

%% Brain sources
sources_brain = [];
sources_brain.pos = geoms.sources_brain.pos;
sources_brain.inside = true(size(sources_brain.pos, 1), 1);
sources_brain.unit = 'mm';
sources_brain = ft_convert_units(sources_brain, 'm');

%% Sensors
grad = geoms.coils_3axis;
grad = ft_convert_units(grad, 'm');
grad = ft_datatype_sens(grad);

%% Compute leadfields 
cfg = [];
cfg.grad       = grad;
cfg.headmodel  = vol_cord;
cfg.sourcemodel = sources_spine;
cfg.reducerank = 'no';
cfg.channel    = 'all';
cfg.normalize  = 'no';
leadfield_cord = ft_prepare_leadfield(cfg);

% Scale to nAm
leadfield_cord.leadfield = cellfun(@(x) x * 1e-9, leadfield_cord.leadfield, ...
                              'UniformOutput', false);
leadfield_cord.unit = 'T/nAm';

cfg = [];
cfg.grad       = grad;
cfg.headmodel  = vol_brain;
cfg.sourcemodel = sources_brain;
cfg.reducerank = 'no';
cfg.channel    = 'all';
cfg.normalize  = 'no';
leadfield_brain = ft_prepare_leadfield(cfg);

% Scale to nAm
leadfield_brain.leadfield = cellfun(@(x) x * 1e-9, leadfield_brain.leadfield, ...
                              'UniformOutput', false);
leadfield_brain.unit = 'T/nAm';

%% Save each leadfield component (X, Y, Z)
outdir = 'D:\OPM\sub-OP00210\fwd\cervical_complex_brainspine_sources';
if ~exist(outdir, 'dir'); mkdir(outdir); end
cd(outdir)

nsensors = size(grad.chanori,1)/3;

% Spinal sources 
for idx = 1:size(sources_spine.pos, 1)
    Lmat = leadfield_cord.leadfield{idx}; 

    for ori_idx = 1:3
        L_ori = Lmat(:, ori_idx);  

        Ls_x = L_ori(1:nsensors);                        % Sensor X-axis
        Ls_y = L_ori(nsensors+1:2*nsensors);             % Sensor Y-axis
        Ls_z = L_ori((2*nsensors)+1:3*nsensors);         % Sensor Z-axis

        fname = sprintf('spine_%d_ori%d.mat', idx, ori_idx);
        save(fname, 'Ls_x', 'Ls_y', 'Ls_z');
    end
end


% Brain sources 
for idx = 1:size(sources_brain.pos, 1)
    Lmat = leadfield_brain.leadfield{idx}; 

    for ori_idx = 1:3
        L_ori = Lmat(:, ori_idx);  

        Ls_x = L_ori(1:nsensors);                        % Sensor X-axis
        Ls_y = L_ori(nsensors+1:2*nsensors);             % Sensor Y-axis
        Ls_z = L_ori((2*nsensors)+1:3*nsensors);         % Sensor Z-axis

        fname = sprintf('brain_%d_ori%d.mat', idx, ori_idx);
        save(fname, 'Ls_x', 'Ls_y', 'Ls_z');
    end
end