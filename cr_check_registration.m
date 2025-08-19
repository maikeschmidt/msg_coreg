function output_meshes = cr_check_registration(S)

if ~isfield(S, 'subject'); error('Please provide a subject mesh'); end
if ~isfield(S, 'sensors'); S.sensors = []; end
if ~isfield(S, 'torso_mode'); error('Please set torso_mode = canonical or anatomical'); end
if ~isfield(S, 'spine_mode'); S.spine_mode = 'full'; end

% --- Step 1: get transform
switch lower(S.torso_mode)
    case 'canonical'
        % Prompt fiducials
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
% meshes = cr_load_meshes(T, true, 'mri_full_spine', 'mri_full_bone', 'mri_torso', 'mri_lungs', 'heart');

meshes = cr_load_meshes(T, true, spineType, boneType, torsoType, lungType, heartType);


figure('Name','Registration check','Color','w'); hold on;
meshNames = fieldnames(meshes);
colors = lines(numel(meshNames)+1);

% Subject surface (grey)
patch('Vertices', S.subject.vertices, ...
      'Faces', S.subject.faces, ...
      'FaceAlpha',0.1,'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
legendEntries = {'Subject'};

% Simulation meshes
for i = 1:numel(meshNames)
    mesh_i = meshes.(meshNames{i});
    patch('Vertices', mesh_i.vertices, ...
          'Faces', mesh_i.faces, ...
          'FaceAlpha',0.3,'EdgeColor','none','FaceColor',colors(i,:));
    legendEntries{end+1} = meshNames{i}; %#ok<AGROW>
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