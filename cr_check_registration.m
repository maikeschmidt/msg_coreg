function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'torso_mode'); error('Please set torso_mode = canonical or anatomical'); end
if ~isfield(S, 'spine_mode'); S.spine_mode = 'full'; end

% --- Step 1: get transform
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

% --- Step 2: decide filenames
switch lower(S.torso_mode)
    case 'canonical'
        heartType = 'canonical_heart';
        lungType  = 'canonical_lungs';
        torsoType = 'canonical_torso';
        if strcmpi(S.spine_mode, 'cervical')
            spineType = 'cervical_spine';
            boneType  = 'cervical_bone';
        else
            spineType = 'spine';
            boneType  = 'canonical_bone';
        end
    case 'anatomical'
        heartType = 'heart';
        lungType  = 'mri_lungs';
        torsoType = 'mri_torso';
        if strcmpi(S.spine_mode, 'cervical')
            spineType = 'mri_cervical_spine';
            boneType  = 'mri_cervical_bone';
        else
            spineType = 'mri_full_spine';
            boneType  = 'mri_full_bone';
        end
end


% load meshes
meshes = cr_load_meshes(T, true, spineType, boneType, torsoType, lungType, heartType);

% brain registration
if isfield(S, 'brain') && S.brain
    disp('Please select three fiducials on the subject head: NAS, LPA, RPA');
    brain_fids_select = spm_mesh_select(S.subject); 
    brain_fids = brain_fids_select';  
    
    % regS.subject   = S.subject;
    regS.fiducials = brain_fids;
    regS.head      = S.subject; 
    regS.plot      = false;
    
    temp_brain = cr_register_brain(regS);
    
    meshes.brain  = temp_brain.brain;
    meshes.scalp  = temp_brain.scalp;
    meshes.iskull = temp_brain.iskull;
    meshes.oskull = temp_brain.oskull;
end


figure('Name','Registration check','Color','w'); hold on;

% Subject surface (grey)
patch('Vertices', S.subject.vertices, ...
      'Faces', S.subject.faces, ...
      'EdgeAlpha',0.1,'EdgeColor','[0.8 0.8 0.8]','FaceColor','none');
legendEntries = {'Subject'};

% Simulation meshes
meshNames = fieldnames(meshes);
for i = 1:numel(meshNames)
    if contains(lower(meshNames{i}), 'scalp')
        continue;  
    end
    
    mesh_i = meshes.(meshNames{i});
    
    % Set colors based on mesh type using contains
    name = lower(meshNames{i});
    if contains(name, 'torso')
        c = [0.7 0.6 0.9];       % light purple
    elseif contains(name, 'spine')
        c = [1.0 0.0 0.0];       % red
    elseif contains(name, 'bone')
        c = [1.0 1.0 0.0];       % yellow
    elseif contains(name, 'lung')
        c = [0.0 0.0 1.0];       % blue
    elseif contains(name, 'heart')
        c = [0.0 1.0 0.0];       % green
    elseif contains(name, 'brain')
        c = [0.3 0.3 0.9];       % deep blue
    elseif contains(name, 'iskull')
        c = [0.8 0.8 0.8];       % grey
    elseif contains(name, 'oskull')
        c = [0.5 0.5 0.5];       % darker grey
    else
        c = [0.5 0.5 0.5];       % default grey
    end
    
    patch('Vertices', mesh_i.vertices, ...
          'Faces', mesh_i.faces, ...
          'FaceAlpha',0.3,'EdgeColor','none','FaceColor',c);
    
    legendEntries{end+1} = meshNames{i}; 
end

% Sensors (if provided)
if ~isempty(S.sensors)
    ft_plot_sens(S.sensors)
    legendEntries{end+1} = 'Sensors';
end

axis equal; grid on; view(3);
legend(legendEntries, 'Interpreter','none');
title(sprintf('Registration check (%s, %s spine)', S.torso_mode, S.spine_mode));

output_meshes = meshes;
output_meshes.transform = T;

end