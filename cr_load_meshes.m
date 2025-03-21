function [meshes] = cr_load_meshes(T, loadSpine)
    % Load meshes with transformations applied
    % T: Transformation matrix
    % loadSpine: Boolean flag to determine whether to load the canonical spine

    if nargin < 1 || isempty(T)
        T = eye(4); 
    end
    if nargin < 2
        loadSpine = true; % Default: Load spine unless overridden
    end

    names = {'heart', 'lungs', 'torso'}; % Always load these
    if loadSpine
        names{end + 1} = 'spine'; % Add spine only if needed
    end

    % Load and transform all except spine (if user-specified)
    applyTransform = true; % Always apply transform to heart, lungs, torso
    meshes = load_and_transform_meshes(names, T, applyTransform);
end



function meshes = load_and_transform_meshes(fileNames, transformMatrix, applyTransform)
    meshes = {};
    for i = 1:numel(fileNames)
        stlFile = fullfile(coreg_path, [fileNames{i}, '.stl']);
        mesh = [];

        if exist(stlFile, 'file')
            tempMesh = stlread(stlFile);
            mesh.vertices = tempMesh.vertices;
            mesh.faces = double(tempMesh.faces);
        else
            warning('File %s not found. Skipping.', stlFile);
            continue;
        end

        % Apply transformation only if specified
        if applyTransform
            mesh = spm_mesh_transform(mesh, transformMatrix);
        end

        % Append mesh
        meshes{end + 1} = mesh;
    end
end
