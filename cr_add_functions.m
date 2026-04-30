% cr_add_functions - Check and configure all required toolbox dependencies
%
% Detects and adds all necessary functions and libraries to the MATLAB path
% for the MSG Coregistration Toolbox to run. Specifically:
%   - Adds the Helsinki BEM (HBF) library (cloning via git submodule if absent)
%   - Creates lowercase wrapper functions for FieldTrip compatibility
%   - Copies required private HBF/FieldTrip functions to an accessible folder
%
% This function should be called once at the start of a session before
% using any other toolbox functions.
%
% USAGE:
%   cr_add_functions()
%
% INPUT:
%   None
%
% OUTPUT:
%   None — modifies the MATLAB path in place and copies files as needed
%
% DEPENDENCIES:
%   - coreg_path()         : locates the repository root directory
%   - hbf_SetPaths()       : sets up Helsinki BEM library paths
%   - FieldTrip            : ft_defaults must be on the MATLAB path
%   - Git                  : required if HBF submodule has not been initialised
%
% NOTES:
%   - The HBF library is expected at: <repo_root>/hbf_lc_p
%   - FieldTrip wrapper functions are written to: <repo_root>/hbf_ft_wrappers
%   - Private helper functions are copied to:     <repo_root>/hbf_private
%   - Wrappers are only created where lowercase versions are missing,
%     to avoid conflicts with existing FieldTrip installations
%
% EXAMPLE:
%   cr_add_functions()   % run once at the start of your session
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
%
function cr_add_functions
% Add Helsinki BEM library and ensure HBF helper functions are accessible

% Paths
hbf_path = fullfile(coreg_path, 'hbf_lc_p');
addpath(hbf_path);

if isempty(which('hbf_SetPaths'))
    disp('Cloning BEM library to repository!');
    !git submodule update --init
end
hbf_SetPaths;

% Patch FieldTrip HBF detection (case-sensitivity fix)
% ft_hastoolbox lower()s function names, but HBF uses uppercase filenames
% We create lowercase wrappers that forward to the uppercase originals

wrapper_dir = fullfile(coreg_path, 'hbf_ft_wrappers');
if ~exist(wrapper_dir,'dir'), mkdir(wrapper_dir); end

% Functions FieldTrip checks for HBF
hbf_ft_deps = {
    'hbf_tm_phi_lc',        'HBF_TM_PHI_LC'
    'hbf_lfm_b_lc',         'HBF_LFM_B_LC'
    'hbf_bemoperatorsb_linear', 'HBF_BEMOperatorsB_Linear'
};

for i = 1:size(hbf_ft_deps,1)
    lower_name = hbf_ft_deps{i,1};
    upper_name = hbf_ft_deps{i,2};

    % Only create wrapper if lowercase version is missing
    if isempty(which(lower_name)) && ~isempty(which(upper_name))
        wrapper_file = fullfile(wrapper_dir, [lower_name '.m']);

        fid = fopen(wrapper_file,'w');
        fprintf(fid, ...
            'function varargout = %s(varargin)\n' + ...
            '%% Auto-generated FieldTrip compatibility wrapper\n' + ...
            '[varargout{1:nargout}] = %s(varargin{:});\n' + ...
            'end\n', ...
            lower_name, upper_name);
        fclose(fid);

        fprintf('Created HBF wrapper: %s -> %s\n', lower_name, upper_name);
    end
end

addpath(wrapper_dir);


% Destination for private HBF/FT functions (normal folder, not private)
dir_out = fullfile(coreg_path, 'hbf_private');
if ~exist(dir_out,'dir'), mkdir(dir_out); end

% Functions to copy
fnames = {'hbf_LFM_B_LC_xyz', 'hbf_Phiinf_xyz', 'hbf_Binf_xyz', 'bmesh2bnd', 'fixbalance'};
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

fprintf('HBF helper functions copied to %s and added to path.\n', dir_out);

end

