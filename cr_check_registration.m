function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'T'); warning('No transformation matrix found'); S.T = eye(4); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'spine'); S.spine = []; end
if ~isfield(S, 'vertebrae'); S.vertebrae = []; end

% Determine spine loading mode
if ~isfield(S, 'spine_mode')
    S.spine_mode = 'default';
end

% If custom meshes are provided, override spine_mode
if ~isempty(S.spine) || ~isempty(S.vertebrae)
    S.spine_mode = 'custom';
end

% Set spine and bone types based on mode
switch lower(S.spine_mode)
    case 'default'
        loadSpine = true;
        spineType = 'spine';
        boneType = 'bone';
    case 'cervical'
        loadSpine = true;
        spineType = 'cervical_spine';
        boneType = 'cervical_bone';
    case 'custom'
        loadSpine = false;
        spineType = '';
        boneType = '';
    otherwise
        error('Unknown spine_mode: %s', S.spine_mode);
end

% Load canonical or default meshes
meshes = cr_load_meshes(S.T, loadSpine, spineType, boneType);
unit = tt_determine_mesh_units(meshes);

% Plot subject mesh (torso/body)
figure
ft_plot_mesh(S.subject, 'facecolor', 'none', 'edgecolor', 'k', ...
    'clipping', 'off', 'edgealpha', 0.05); % Grey subject
hold on

% Plot canonical meshes (heart, lungs, torso)
structure_names = {'heart', 'lungs', 'torso'};
structure_colors = {[0 0 1], [0 1 0], [0.5 0 0.5]}; % blue, green, purple
structure_alphas = [0.3, 0.3, 0.1];

for ii = 1:3
    tmp.pnt = meshes{ii}.vertices;
    tmp.tri = meshes{ii}.faces;
    tmp.unit = unit;

    ft_plot_mesh(tmp, 'facecolor', structure_colors{ii}, ...
        'edgecolor', 'none', 'facealpha', structure_alphas(ii));

    switch ii
        case 1, output_meshes.mesh_heart = tmp;
        case 2, output_meshes.mesh_lungs = tmp;
        case 3, output_meshes.mesh_torso = tmp;
    end
end

% Plot spine
if ~isempty(S.spine)
    spine_tmp.pnt = S.spine.vertices;
    spine_tmp.tri = S.spine.faces;
    spine_tmp.unit = unit;
    ft_plot_mesh(spine_tmp, 'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.5);
    output_meshes.mesh_spine = spine_tmp;
elseif loadSpine && numel(meshes) > 3
    spine_tmp.pnt = meshes{4}.vertices;
    spine_tmp.tri = meshes{4}.faces;
    spine_tmp.unit = unit;
    ft_plot_mesh(spine_tmp, 'facecolor', 'red', 'edgecolor', 'none', 'facealpha', 0.5);
    output_meshes.mesh_spine = spine_tmp;
end

% Plot bone/vertebrae
if ~isempty(S.vertebrae)
    bone_tmp.pnt = S.vertebrae.vertices;
    bone_tmp.tri = S.vertebrae.faces;
    bone_tmp.unit = unit;
    ft_plot_mesh(bone_tmp, 'facecolor', 'yellow', 'edgecolor', 'none', 'facealpha', 0.2);
    output_meshes.mesh_vertebrae = bone_tmp;
elseif loadSpine && numel(meshes) > 4
    bone_tmp.pnt = meshes{5}.vertices;
    bone_tmp.tri = meshes{5}.faces;
    bone_tmp.unit = unit;
    ft_plot_mesh(bone_tmp, 'facecolor', 'yellow', 'edgecolor', 'none', 'facealpha', 0.2);
    output_meshes.mesh_vertebrae = bone_tmp;
end

xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
grid on

if ~isempty(S.sensors)
    ft_plot_sens(ft_convert_units(S.sensors, unit))
end

end

