function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'T'); warning('No transformation matrix found'); S.T = eye(4); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'spine'); S.spine = []; end
if ~isfield(S, 'vertebrae'); S.vertebrae = []; end

% Plot the subject mesh (torso/body)
figure
ft_plot_mesh(S.subject, 'facecolor', 'none', 'edgecolor', 'y', ...
    'clipping', 'off', 'edgealpha', 0.2);
hold on

% Determine whether to load the canonical spine
loadCanonicalSpine = isempty(S.spine) && isempty(S.vertebrae);
meshes = cr_load_meshes(S.T, loadCanonicalSpine);
unit = tt_determine_mesh_units(meshes);

% Define colors and alpha values
colors = {'r', 'g', 'b'}; % Heart, lungs, torso
alphas = [0.3, 0.3, 0.1];

% Assign heart, lungs, and torso
for ii = 1:3
    tmp.pnt = meshes{ii}.vertices;
    tmp.tri = meshes{ii}.faces;
    tmp.unit = unit;

    % Plot mesh
    ft_plot_mesh(tmp, 'facecolor', colors{ii}, 'edgecolor', 'none', 'facealpha', alphas(ii));

    % Save mesh data
    switch ii
        case 1
            output_meshes.mesh_heart = tmp;
        case 2
            output_meshes.mesh_lungs = tmp;
        case 3
            output_meshes.mesh_torso = tmp;
    end
end

% Handle spine and vertebrae cases
if ~isempty(S.spine)
    % Use user-provided spine (no transformation)
    spine_tmp.pnt = S.spine.vertices;
    spine_tmp.tri = S.spine.faces;
    spine_tmp.unit = unit;
    ft_plot_mesh(spine_tmp, 'facecolor', 'blue', 'edgecolor', 'none', 'facealpha', 0.5);
    output_meshes.mesh_spine = spine_tmp;
elseif loadCanonicalSpine && numel(meshes) > 3
    % Use canonical spine (transformed)
    spine_tmp.pnt = meshes{4}.vertices;
    spine_tmp.tri = meshes{4}.faces;
    spine_tmp.unit = unit;
    ft_plot_mesh(spine_tmp, 'facecolor', 'blue', 'edgecolor', 'none', 'facealpha', 0.5);
    output_meshes.mesh_spine = spine_tmp;
end

if ~isempty(S.vertebrae)
    % Use user-provided vertebrae (no transformation)
    bone_tmp.pnt = S.vertebrae.vertices;
    bone_tmp.tri = S.vertebrae.faces;
    bone_tmp.unit = unit;
    ft_plot_mesh(bone_tmp, 'facecolor', 'black', 'edgecolor', 'none', 'facealpha', 0.2);
    output_meshes.mesh_vertebrae = bone_tmp;
end

xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
grid on

% Plot sensors if available
if ~isempty(S.sensors)
    ft_plot_sens(ft_convert_units(S.sensors, unit))
end
end
