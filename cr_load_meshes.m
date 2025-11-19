function meshes = cr_load_meshes(T, loadSpine, spineType, boneType, torsoType, lungType, heartType, includeMuscle)

if nargin < 1 || isempty(T), T = eye(4); end
if nargin < 2, loadSpine = true; end
if nargin < 8, includeMuscle = false; end

names = {heartType, lungType, torsoType};

if loadSpine
    names{end+1} = spineType;
    names{end+1} = boneType;
end

if includeMuscle
    names{end+1} = 'back_muscle_temp';
end

meshes = load_and_transform_meshes(names, T, true);

end

function meshes = load_and_transform_meshes(fileNames, transformMatrix, applyTransform)

    tolerance = 1e-6;
    meshes = struct();  

    % Mapping from variant filenames to core names
    nameMap = containers.Map( ...
        {'realistic_full_bone','realistic_cervical_bone', ...
         'mri_full_inhomo','mri_cervical_inhomo','canonical_full_inhomo','canonical_cervical_inhomo', ...
         'mri_full_homo','mri_cervical_homo','canonical_full_homo','canonical_cervical_homo', ...
         'mri_full_cont','mri_cervical_cont','canonical_full_cont','canonical_cervical_cont', ...
         'mri_torso','canonical_torso', ...
         'mri_lungs','canonical_lungs', ...
         'heart','canonical_heart', ...
         'back_muscle_temp'}, ...
        {'bone','bone', ...
         'bone','bone','bone','bone', ...
         'bone','bone','bone','bone', ...
         'bone','bone','bone','bone', ...
         'torso','torso', ...
         'lungs','lungs', ...
         'heart','heart', ...
         'back_muscle'} ...
    );

    for i = 1:numel(fileNames)
        meshName = fileNames{i};
        stlFile = fullfile(coreg_path, [meshName, '.stl']);
        mesh = [];

        % Load STL
        if ~exist(stlFile, 'file')
            warning('File %s not found. Skipping.', stlFile);
            continue;
        end

        tempMesh = stlread(stlFile);

        % Detect STL structure format
        if isstruct(tempMesh)

            if isfield(tempMesh, 'vertices') && isfield(tempMesh, 'faces')
                V = tempMesh.vertices;
                F = double(tempMesh.faces);

            elseif isfield(tempMesh, 'pos') && isfield(tempMesh, 'tri')
                V = tempMesh.pos;
                F = double(tempMesh.tri);

            else
                error('stlread: Unsupported STL structure format. Fields found: %s', strjoin(fieldnames(tempMesh)));
            end

        elseif isnumeric(tempMesh)
            % Very old stlread versions return only vertices
            if size(tempMesh,2) == 3
                V = tempMesh;
                error('stlread returned vertex-only matrix. Faces must be generated manually.');
            else
                error('stlread: Unrecognized numeric matrix output.');
            end
        else
            error('stlread: Unexpected output type: %s', class(tempMesh));
        end

        % Deduplicate vertices
        V_rounded = round(V * (1 / tolerance)) * tolerance;
        [V_unique, ~, ic] = unique(V_rounded, 'rows', 'stable');
        F_fixed = ic(F);

        % Assign mesh
        mesh.vertices = V_unique;
        mesh.faces    = F_fixed;
        mesh.unit     = 'mm';

        % Optional transform
        if applyTransform
            mesh = spm_mesh_transform(mesh, transformMatrix);
        end

        % Map variant name → core name
        if isKey(nameMap, meshName)
            coreName = nameMap(meshName);
        else
            coreName = meshName; % fallback
        end

        meshes.(coreName) = mesh;
    end
end
