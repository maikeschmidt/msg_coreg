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

    % --- Mapping from variant filenames to core names ---
    nameMap = containers.Map( ...
        {'mri_full_spine', 'mri_cervical_spine', 'spine', 'cervical_spine', ...
         'mri_full_bone', 'mri_cervical_bone', 'canonical_bone', 'cervical_bone', ...
         'mri_torso', 'canonical_torso', ...
         'mri_lungs', 'canonical_lungs', ...
         'heart', 'canonical_heart'}, ...
        {'spine', 'spine', 'spine', 'spine', ...
         'bone', 'bone', 'bone', 'bone', ...
         'torso', 'torso', ...
         'lungs', 'lungs', ...
         'heart', 'heart'} ...
    );

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

        % Map variant name to core name
        if isKey(nameMap, meshName)
            coreName = nameMap(meshName);
        else
            coreName = meshName; % fallback
        end

        meshes.(coreName) = mesh;
    end
end
