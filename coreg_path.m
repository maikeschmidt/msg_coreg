% Automatically detects the file path for where the msg_coreg
% repository is saved
%
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Date: April 2026
%

function folderpath = coreg_path
folderpath = fileparts(mfilename('fullpath'));
end
