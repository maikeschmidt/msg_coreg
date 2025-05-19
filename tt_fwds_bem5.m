function L = tt_fwds_bem5(S)
% Generate the lead fields for a 5-shell BEM, adding spinal cord and vertebrae

% Input validation
if ~isfield(S,'pos'); error('please specify the source positions!'); end
if ~isfield(S,'posunits'); error('please specify the current positions units!'); end
if ~isfield(S,'ori'); S.ori = []; end
if ~isfield(S,'sensors'); error('please specify the sensor structure!'); end
if ~isfield(S,'T'); S.T = eye(4); end
if ~isfield(S,'names'); S.names = {'spinalcord','vertebrae','blood','lungs','torso'}; end
if ~isfield(S, 'cord'); S.cord = []; end
if ~isfield(S, 'vertebrae'); S.vertebrae = []; end
if ~isfield(S,'ci'); S.ci = [.33 .007 .62 .05 .23]; end
if ~isfield(S,'co'); S.co = [.23 .23 .23 .23 0]; end
if ~isfield(S,'isa'); S.isa = []; end

% Ensure necessary functions are loaded
if isempty(which('hbf_BEMOperatorsPhi_LC'))
    tt_add_bem;
end

%load meshes
meshes = fem_load_meshes(S.T);
[~, sf] = tt_determine_mesh_units(meshes);
bmeshes = {}; 

if ~isempty(S.cord)
    cord_tmp = [];
    cord_tmp.p = S.cord.vertices / sf; 
    cord_tmp.e = S.cord.faces; 
    bmeshes{end+1} = cord_tmp; 
end

if ~isempty(S.vertebrae)
    vert_tmp = [];
    vert_tmp.p = S.vertebrae.vertices / sf; 
    vert_tmp.e = S.vertebrae.faces; 
    bmeshes{end+1} = vert_tmp; 
end


for ii = 1:numel(meshes)
    bmeshes{end+1}.p = meshes{ii}.vertices / sf; % Convert to meters
    bmeshes{end}.e = meshes{ii}.faces; % Mesh faces
end

% Convert sensor positions to meters
S.sensors = ft_convert_units(S.sensors, 'm');

% Prepare the source model and convert units
cfg = [];
cfg.method = 'basedonpos';
cfg.sourcemodel.pos = S.pos;
cfg.sourcemodel.unit = S.posunits;
src = ft_prepare_sourcemodel(cfg);
src = ft_convert_units(src, 'm');

% Coil positions and orientations
coils = [];
coils.p = S.sensors.coilpos;
coils.n = S.sensors.coilori;

% Conductivity parameters for each shell
ci = S.ci; 
co = S.co; 

% Generate the transfer matrix for BEM
D = hbf_BEMOperatorsPhi_LC(bmeshes);
if isempty(S.isa)
    Tphi_full = hbf_TM_Phi_LC(D, ci, co);
else
    fprintf('%-40s: %30s\n', 'Applying ISA', S.names{S.isa});
    Tphi_full = hbf_TM_Phi_LC_ISA2(D, ci, co, S.isa);
end

% Generate the B-field transfer matrix
DB = hbf_BEMOperatorsB_Linear(bmeshes, coils);
TB = hbf_TM_Bvol_Linear(DB, Tphi_full, ci, co);

% Compute the lead field
if isempty(S.ori)
    L = S.sensors.tra * hbf_LFM_B_LC(bmeshes, coils, TB, src.pos);
else
    L = S.sensors.tra * hbf_LFM_B_LC(bmeshes, coils, TB, S.pos, S.ori);
end
