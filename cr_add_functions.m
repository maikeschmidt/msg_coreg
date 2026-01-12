function cr_add_functions
% Add Helsinki BEM library and ensure HBF helper functions are accessible
% Add fieldtrp and private functions to path that are needed to run these
% simulations

% Paths
hbf_path = fullfile(coreg_path, 'hbf_lc_p');
addpath(hbf_path);

if isempty(which('hbf_SetPaths'))
    disp('Cloning BEM library to repository!');
    !git submodule update --init
end
hbf_SetPaths;

% Destination for private HBF/FT functions (normal folder, not private)
dir_out = fullfile(coreg_path, 'hbf_private');
if ~exist(dir_out,'dir'), mkdir(dir_out); end

% Functions to copy
fnames = { ...
    'hbf_LFM_B_LC_xyz', ...
    'hbf_Phiinf_xyz', ...
    'hbf_Binf_xyz', ...
    'bmesh2bnd', ...
    'fixbalance', ...
    'filetype_check_uri' ...
    'filetype_check_header',...
};

exts   = {'m','p'};

% HBF private folder
dir_in = fullfile(hbf_path, 'hbf_calc', 'private');

% FieldTrip path
ft_path = fileparts(which('ft_defaults'));
ft_private = dir(fullfile(ft_path, '**', 'private'));

% Loop over each function
for ii = 1:numel(fnames)
    copied = false;
    for jj = 1:numel(exts)
        % First, try HBF folder
        fin = fullfile(dir_in, [fnames{ii} '.' exts{jj}]);
        if exist(fin, 'file')
            copyfile(fin, fullfile(dir_out, [fnames{ii} '.' exts{jj}]));
            copied = true;
            break;
        end

        % Then, check FieldTrip private folders
        if ~copied
            for k = 1:numel(ft_private)
                fin_ft = fullfile(ft_private(k).folder, [fnames{ii} '.' exts{jj}]);
                if exist(fin_ft, 'file')
                    copyfile(fin_ft, fullfile(dir_out, [fnames{ii} '.' exts{jj}]));
                    copied = true;
                    break;
                end
            end
        end
    end
    if ~copied
        warning('Could not find %s in HBF or FieldTrip private folders.', fnames{ii});
    end
end

% Finally, add the new folder to MATLAB path
addpath(dir_out);

fprintf('helper functions copied to %s and added to path.\n', dir_out);

end

