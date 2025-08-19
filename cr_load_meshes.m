function meshes = cr_load_meshes(T, loadSpine, spineType, boneType, torsoType, lungType, heartType)

if nargin < 1 || isempty(T), T = eye(4); end
if nargin < 2, loadSpine = true; end

names = {heartType, lungType, torsoType};

if loadSpine
    names{end+1} = spineType;
    names{end+1} = boneType;
end

meshes = load_and_transform_meshes(names, T, true);

end


function meshes = load_and_transform_meshes(fileNames, transformMatrix, applyTransform)

    tolerance = 1e-6;
    meshes = struct();  

    for i = 1:numel(fileNames)
        meshName = fileNames{i};
        stlFile = fullfile(coreg_path, [meshName, '.stl']);
        mesh = [];

        if exist(stlFile, 'file')
            tempMesh = stlread(stlFile);

            V = tempMesh.vertices;
            F = double(tempMesh.faces);

            % Deduplicate vertices
            [V_unique, ~, ic] = unique(round(V * (1 / tolerance)) * tolerance, 'rows', 'stable');
            F_fixed = ic(F);

            mesh.vertices = V_unique;
            mesh.faces = F_fixed;
            mesh.unit = 'mm';
        else
            warning('File %s not found. Skipping.', stlFile);
            continue;
        end

        if applyTransform
            mesh = spm_mesh_transform(mesh, transformMatrix);
        end

        % Save under name
        meshes.(meshName) = mesh;
    end
end

