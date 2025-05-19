function meshes = cr_load_meshes(T, loadSpine, spineType, boneType)

if nargin < 1 || isempty(T)
    T = eye(4);
end
if nargin < 2
    loadSpine = true;
end
if nargin < 3
    spineType = 'spine';
end
if nargin < 4
    boneType = 'bone';
end

names = {'heart', 'lungs', 'torso'};

if loadSpine
    names{end + 1} = spineType;
    names{end + 1} = boneType;
end

applyTransform = true;
meshes = load_and_transform_meshes(names, T, applyTransform);

end



function meshes = load_and_transform_meshes(fileNames, transformMatrix, applyTransform)

    tolerance = 1e-6;

    meshes = {};
    for i = 1:numel(fileNames)
        stlFile = fullfile(coreg_path, [fileNames{i}, '.stl']);
        mesh = [];

        if exist(stlFile, 'file')
            tempMesh = stlread(stlFile);

            V = tempMesh.vertices;
            F = double(tempMesh.faces);
            [V_unique, ~, ic] = unique(round(V * (1 / tolerance)) * tolerance, 'rows', 'stable');
            F_fixed = ic(F);

            mesh.vertices = V_unique;
            mesh.faces = F_fixed;
        else
            warning('File %s not found. Skipping.', stlFile);
            continue;
        end

        if applyTransform
            mesh = spm_mesh_transform(mesh, transformMatrix);
        end

        meshes{end + 1} = mesh;
    end
end
