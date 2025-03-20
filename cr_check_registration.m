function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'T'); warning('No transformation matrix found'); S.T = eye(4); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'spine'); S.spine = []; end
if ~isfield(S, 'vertebrae'); S.vertebrae = []; end

% Plot the subject mesh (torso/body)
figure
ft_plot_mesh(S.subject, 'facecolor', 'none', 'edgecolor', 'k', ...
    'clipping', 'off', 'edgealpha', 0.2);
hold on

% Load heart and lungs with transformation, torso without transformation
meshes = cr_load_meshes(S.T);
unit = cr_determine_mesh_units(meshes);

% Ensure the meshes are correctly assigned
if numel(meshes) < 3
    error('Expected 3 meshes (heart, lungs, torso), but found only %d.', numel(meshes));
end

% Define colors for visualization
colors = {'r', 'g', 'b', 'y'};
alphas = [0.3, 0.3, 0.3, 0.8];

for ii = 1:numel(meshes)
    tmp.pnt = meshes{ii}.vertices;
    tmp.tri = meshes{ii}.faces;
    tmp.unit = unit;

    % Plot mesh with color and alpha
    ft_plot_mesh(tmp, 'facecolor', colors{ii}, 'edgecolor', 'none', 'facealpha', alphas(ii));

    % Save the mesh data
    switch ii
        case 1
            output_meshes.mesh_heart.vertices = tmp.pnt;
            output_meshes.mesh_heart.faces = tmp.tri;
        case 2
            output_meshes.mesh_lungs.vertices = tmp.pnt;
            output_meshes.mesh_lungs.faces = tmp.tri;
        case 3
            output_meshes.mesh_torso.vertices = tmp.pnt;
            output_meshes.mesh_torso.faces = tmp.tri;
    end
end
end

