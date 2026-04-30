% cr_load_meshes - Load and transform anatomical mesh components for BEM simulation
%
% Loads a set of anatomical STL meshes (torso, lungs, heart, spine, bone,
% and optionally back muscle) from the repository meshes folder, applies a
% 4x4 rigid-body transform to register them to subject space, and returns
% them as a struct of mesh components ready for BEM forward modelling.
%
% USAGE:
%   meshes = cr_load_meshes(T, loadSpine, spineType, boneType, ...
%                           torsoType, lungType, heartType, includeMuscle)
%
% INPUT:
%   T             - 4x4 rigid-body transform matrix (canonical → subject
%                   space). Pass [] or omit to use identity (default: eye(4))
%   loadSpine     - Logical; include spine and bone meshes (default: true)
%   spineType     - Spine mesh filename string, e.g.:
%                     'spine' | 'mri_full_spine' | 'mri_cervical_spine' |
%                     'cervical_spine'
%   boneType      - Bone mesh filename string, e.g.:
%                     'realistic_full_bone' | 'mri_full_inhomo' |
%                     'canonical_full_homo' | 'mri_cervical_cont'
%   torsoType     - Torso mesh filename string:
%                     'mri_torso' | 'canonical_torso'
%   lungType      - Lung mesh filename string:
%                     'mri_lungs' | 'canonical_lungs'
%   heartType     - Heart mesh filename string:
%                     'heart' | 'canonical_heart'
%   includeMuscle - Logical; include back muscle mesh (default: false)
%                   Only available for anatomical torso mode
%
% OUTPUT:
%   meshes        - Struct with fields for each loaded mesh component.
%                   Field names are normalised to core names regardless of
%                   variant filename used:
%                     .torso       - Torso surface mesh
%                     .lungs       - Lung mesh
%                     .heart       - Heart mesh
%                     .bone        - Bone/spine bone mesh
%                     .spine       - Spine centreline mesh
%                     .back_muscle - Back muscle mesh (if includeMuscle)
%                   Each mesh is a struct with:
%                     .vertices    - [N x 3] vertex coordinates (mm)
%                     .faces       - [M x 3] face indices
%                     .unit        - 'mm'
%
% DEPENDENCIES:
%   - coreg_path()           : locates the repository meshes folder
%   - spm_mesh_transform()   : applies the 4x4 transform to each mesh
%   - stlread()              : reads STL files (supports modern and legacy
%                              MATLAB formats)
%
% NOTES:
%   - All STL files are expected in <repo_root>/meshes/
%   - Vertices are deduplicated with a tolerance of 1e-6 before output
%   - Variant mesh filenames (e.g. 'mri_full_inhomo') are mapped to core
%     struct field names (e.g. 'bone') automatically
%   - Missing STL files produce a warning and are skipped rather than
%     causing an error
%
% EXAMPLE:
%   T = cr_register_torso(regS);
%   meshes = cr_load_meshes(T, true, 'mri_full_spine', 'mri_full_homo', ...
%                           'mri_torso', 'mri_lungs', 'heart', false);
%
% REPOSITORY:
%   https://github.com/maikeschmidt/msg_coreg
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
% Date:   April 2026
%
% This file is part of the MSG Coregistration Toolbox.

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
        stlFile = fullfile(coreg_path, 'meshes', [meshName, '.stl']);
        mesh = [];

        % Load STL
        if ~exist(stlFile, 'file')
            warning('File %s not found. Skipping.', stlFile);
            continue;
        end

        tempMesh = stlread(stlFile);

        % Detect STL structure format
                tempMesh = stlread(stlFile);
        
                % Detect STL output format
                if isa(tempMesh, 'triangulation')
                    % Modern MATLAB stlread → triangulation object
                    V = tempMesh.Points;
                    F = tempMesh.ConnectivityList;
        
                elseif isstruct(tempMesh)
        
                    % Struct-based STL readers
                    if isfield(tempMesh, 'vertices') && isfield(tempMesh, 'faces')
                        V = tempMesh.vertices;
                        F = double(tempMesh.faces);
        
                    elseif isfield(tempMesh, 'pos') && isfield(tempMesh, 'tri')
                        V = tempMesh.pos;
                        F = double(tempMesh.tri);
        
                    else
                        error(['Unsupported STL structure format. Fields found: ', ...
                               strjoin(fieldnames(tempMesh), ', ')]);
                    end
        
                elseif isnumeric(tempMesh)
        
                    % Very old STL readers returning only vertices
                    if size(tempMesh,2) == 3
                        error(['stlread returned only vertices. Faces are missing and ', ...
                               'must be reconstructed or a different stlread must be used.']);
                    else
                        error('Unrecognized numeric STL output format.');
                    end
        
                else
                    error('Unsupported output type from stlread: %s', class(tempMesh));
                end
        
                % Deduplicate vertices
                V_rounded = round(V * (1 / tolerance)) * tolerance;
                [V_unique, ~, ic] = unique(V_rounded, 'rows', 'stable');
                F_fixed = ic(F);
        
                % Build mesh
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
