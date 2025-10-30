function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'torso_mode'); error('Please set torso_mode = canonical or anatomical'); end
if ~isfield(S, 'spine_mode'); S.spine_mode = 'full'; end
if ~isfield(S, 'bone_mode'); S.bone_mode = 'default'; end
%paramteres needed if wanting to generate a sensor array - defualt = false
if ~isfield(S, 'sensor_gen'); S.sensor_gen = 'false'; end
if ~isfield(S,'resolution'); S.resolution = 30; end
if ~isfield(S,'depth');      S.depth = -10; end
if ~isfield(S,'margin');     S.margin = 50; end
if ~isfield(S, 'coverage'); S.coverage = 0.6; end

% Step 1: get transform
switch lower(S.torso_mode)
    case 'canonical'
        % Prompt fiducials
        disp('please select three fiducials: left shoulder, right shoulder and chin');
        sim_fids_select = spm_mesh_select(S.subject);
        sim_fids = sim_fids_select';
        regS.subject   = S.subject;
        regS.fiducials = sim_fids;
        regS.plot      = true;
        T = cr_register_torso(regS);
    case 'anatomical'
        T = [0.0000, -0.0000,  1.5385,   18.2670;
                    -0.0000,  1.5385, -0.0000, -173.5329;
                    -1.5385, -0.0000,  0.0000,   46.4977;
                    -0.0000,  0.0000, -0.0000,    1.0000];
    otherwise
        error('Unknown torso_mode %s', S.torso_mode);
end

% Step 2: decide filenames
switch lower(S.torso_mode)
    case 'canonical'
        heartType = 'canonical_heart';
        lungType  = 'canonical_lungs';
        torsoType = 'canonical_torso';
        if strcmpi(S.spine_mode, 'cervical')
            spineType = 'cervical_spine';
        else
            spineType = 'spine';
        end

    case 'anatomical'
        heartType = 'heart';
        lungType  = 'mri_lungs';
        torsoType = 'mri_torso';
        if strcmpi(S.spine_mode, 'cervical')
            spineType = 'mri_cervical_spine';
        else
            spineType = 'mri_full_spine';
        end
end

% Step 3: bone type logic 
switch lower(S.bone_mode)
    case 'realistic'
        if strcmpi(S.torso_mode, 'canonical')
            error('Realistic bone meshes are not available for canonical mode.');
        end
        if strcmpi(S.spine_mode, 'cervical')
            boneType = 'realistic_cervical_bone';
        else
            boneType = 'realistic_full_bone';
        end

    case 'inhomo'
        if strcmpi(S.torso_mode, 'anatomical')
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'mri_cervical_inhomo';
            else
                boneType = 'mri_full_inhomo';
            end
        else % canonical
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'canonical_cervical_inhomo';
            else
                boneType = 'canonical_full_inhomo';
            end
        end

    case 'homo'
        if strcmpi(S.torso_mode, 'anatomical')
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'mri_cervical_homo';
            else
                boneType = 'mri_full_homo';
            end
        else % canonical
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'canonical_cervical_homo';
            else
                boneType = 'canonical_full_homo';
            end
        end

    case 'cont'
        if strcmpi(S.torso_mode, 'anatomical')
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'mri_cervical_cont';
            else
                boneType = 'mri_full_cont';
            end
        else % canonical
            if strcmpi(S.spine_mode, 'cervical')
                boneType = 'canonical_cervical_cont';
            else
                boneType = 'canonical_full_cont';
            end
        end

    otherwise
        error('Unknown bone_mode "%s". Valid options: realistic, inhomo, homo, cont.', S.bone_mode);
end


% Step 4: load meshes
% include vagus only for anatomical torso_mode
includeVagus = strcmpi(S.torso_mode, 'anatomical');
meshes = cr_load_meshes(T, true, spineType, boneType, torsoType, lungType, heartType, includeVagus);
if isfield(meshes, spineType)
    meshes.spine = meshes.(spineType);
    if ~strcmp(spineType, 'spine')
        meshes = rmfield(meshes, spineType); % remove old field
    end
end
torso = meshes.torso;

% Step 4b: optional sensor generation
if isfield(S,'sensor_gen') && (islogical(S.sensor_gen) && S.sensor_gen || ischar(S.sensor_gen) && strcmpi(S.sensor_gen,'true'))
    % Generate back sensors
    S_v3 = [];
    S_v3.subject = S.subject;
    S_v3.T = T;
    S_v3.resolution = S.resolution;
    S_v3.depth      = S.depth;
    S_v3.margin     = S.margin;
    S_v3.coverage = S.coverage;
    S_v3.frontflag  = 1;   % back sensors (flag=1)
    S_v3.triaxial   = 1;
    S_v3.torsotype = S.torso_mode;
    back_sensors = cr_generate_sensor_array_v4(S_v3);

    % Generate front sensors
    S_v3.frontflag = 0; % front sensors (flag=0)
    front_sensors = cr_generate_sensor_array_v4(S_v3);

    % Add to output
    meshes.back_sensors  = back_sensors;
    meshes.front_sensors = front_sensors;
end


% Step 5: brain registration
if isfield(S, 'brain') && S.brain
    disp('Please select three fiducials on the subject head: NAS, LPA, RPA');
    brain_fids_select = spm_mesh_select(torso); 
    brain_fids = brain_fids_select';  

    regS.fiducials = brain_fids;
    regS.head      = torso; 
    regS.plot      = true;

    temp_brain = cr_register_brain(regS);

    meshes.brain  = temp_brain.brain;
    meshes.scalp  = temp_brain.scalp;
    meshes.iskull = temp_brain.iskull;
    meshes.oskull = temp_brain.oskull;
end

% Step 6: plotting
figure('Name','Registration check','Color','w'); hold on;

% === Subject surface (grey outline only) ===
if isfield(S.subject, 'vertices') && isfield(S.subject, 'faces')
    patch('Vertices', S.subject.vertices, ...
          'Faces', S.subject.faces, ...
          'EdgeAlpha', 0.1, ...
          'EdgeColor', [0.8 0.8 0.8], ...
          'FaceColor', 'none');
else
    warning('Subject mesh missing vertices/faces; skipping subject surface plot.');
end
legendEntries = {'Subject'};

% === Simulation meshes ===
meshNames = fieldnames(meshes);
sensorFields = {'front_sensors','back_sensors'};

for i = 1:numel(meshNames)
    if ismember(meshNames{i}, sensorFields)
        continue; % skip sensors
    end

    mesh_i = meshes.(meshNames{i});
    if isempty(mesh_i) || ~isfield(mesh_i,'vertices') || ~isfield(mesh_i,'faces') ...
            || isempty(mesh_i.vertices) || isempty(mesh_i.faces)
        warning('Mesh "%s" is empty or malformed — skipping.', meshNames{i});
        continue;
    end

    % Default colour and transparency
    c = [0.5 0.5 0.5];
    alphaVal = 0.3;

    % Assign per-mesh colour
    name = lower(meshNames{i});
    if contains(name, 'vagus_nerve')
        c = [1.0 0.2 0.8];   % bright magenta
        alphaVal = 1.0;      % fully opaque so it’s clearly visible
    elseif contains(name, 'torso')
        c = [0.7 0.6 0.9];   % light purple
        alphaVal = 0.3;
    elseif contains(name, 'spine')
        c = [1.0 0.0 0.0];   % red
        alphaVal = 0.4;
    elseif contains(name, 'bone')
        c = [1.0 1.0 0.0];   % yellow
        alphaVal = 0.3;
    elseif contains(name, 'lung')
        c = [0.0 0.0 1.0];   % blue
        alphaVal = 0.3;
    elseif contains(name, 'heart')
        c = [0.0 1.0 0.0];   % green
        alphaVal = 0.4;
    elseif contains(name, 'brain')
        c = [0.3 0.3 0.9];   % deep blue
        alphaVal = 0.35;
    elseif contains(name, 'iskull')
        c = [0.8 0.8 0.8];
        alphaVal = 0.35;
    elseif contains(name, 'oskull')
        c = [0.5 0.5 0.5];
        alphaVal = 0.35;
    end

    % Plot this mesh
    patch('Vertices', mesh_i.vertices, ...
          'Faces', mesh_i.faces, ...
          'FaceAlpha', alphaVal, ...
          'EdgeColor', 'none', ...
          'FaceColor', c);

    legendEntries{end+1} = meshNames{i};
end

% === Sensors (if provided) ===
if ~isempty(S.sensors)
    ft_plot_sens(S.sensors)
    legendEntries{end+1} = 'Sensors';
end
if isfield(meshes,'back_sensors') && ~isempty(meshes.back_sensors)
    ft_plot_sens(meshes.back_sensors);
end
if isfield(meshes,'front_sensors') && ~isempty(meshes.front_sensors)
    ft_plot_sens(meshes.front_sensors);
end

% === Final plot settings ===
axis equal; grid on; view(3);
lighting gouraud; camlight;
legend(legendEntries, 'Interpreter','none');
title(sprintf('Registration check (%s torso, %s spine, %s bone)', ...
    S.torso_mode, S.spine_mode, S.bone_mode));

output_meshes = meshes;
output_meshes.transform = T;

end

