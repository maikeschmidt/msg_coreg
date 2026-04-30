% coreg_path - Automatically detect the file path of the msg_coreg repository
%
% Returns the absolute path to the directory where this function file is
% saved. Used internally by the toolbox to locate repository resources
% without requiring manual path configuration.
%
% USAGE:
%   folderpath = coreg_path()
%
% INPUT:
%   None
%
% OUTPUT:
%   folderpath - Character vector containing the absolute path to the
%                directory where coreg_path.m is located
%
% EXAMPLE:
%   path = coreg_path();
%   disp(path)
%   % e.g. '/home/user/toolboxes/msg_coreg'
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

function folderpath = coreg_path
folderpath = fileparts(mfilename('fullpath'));
end
