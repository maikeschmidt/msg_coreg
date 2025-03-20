function [meshes] = cr_load_meshes(T)
    % Load meshes with transformations from predefined file names
    % T: Transformation matrix for 'heart' and 'lungs' only

    if nargin < 1 || isempty(T)
        T = eye(4); 
    end

    names_transformed = {'heart', 'lungs', 'torso', 'spine'};

    meshes = {};

    % Load and transform heart and lungs
    meshes = load_and_transform_meshes(names_transformed, T, meshes, true);

end


function meshes = load_and_transform_meshes(fileNames, transformMatrix, meshes, applyTransform)
    for i = 1:numel(fileNames)
        stlFile = fullfile(coreg_path, [fileNames{i}, '.stl']);
        mesh = [];

        if exist(stlFile, 'file')
            % Load STL (.stl) file using stlread
            tempMesh = stlread(stlFile);
            mesh.vertices = tempMesh.vertices;
            mesh.faces = tempMesh.faces;
        else
            warning('File %s not found as either .gii or .stl. Skipping.', fileNames{i});
            continue;
        end

        % Apply transformation if needed
        if applyTransform
            mesh = spm_mesh_transform(mesh, transformMatrix);
        end

        % Ensure faces are double
        mesh.faces = double(mesh.faces);

        % Append to meshes
        meshes{end + 1} = mesh;
    end
end
