% cr_get_fids - Return anatomical fiducial positions for a given torso type
%
% Provides a set of four landmark positions used to register the canonical
% or anatomical torso mesh to a subject surface. Fiducials are defined in
% the native space of each mesh type.
%
% USAGE:
%   fids = cr_get_fids()
%   fids = cr_get_fids(torsotype)
%
% INPUT:
%   torsotype  - (optional) String specifying the torso mesh type:
%                  'anatomical' - MRI-derived torso mesh (default)
%                  'canonical'  - Template/canonical torso mesh
%
% OUTPUT:
%   fids       - 4x3 matrix of fiducial coordinates [mm]:
%                  Row 1: Left shoulder
%                  Row 2: Right shoulder
%                  Row 3: Chin
%                  Row 4: Lumbar
%
% NOTES:
%   - Fiducial coordinates are hardcoded in the native space of each
%     mesh type and should not be changed unless the reference meshes
%     are updated
%   - These fiducials are used by cr_register_torso() and
%     cr_generate_sensor_array_v4() to compute the torso-to-subject
%     transform
%
% EXAMPLE:
%   fids = cr_get_fids('canonical');
%   fids = cr_get_fids('anatomical');  % default
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


function fids = cr_get_fids(torsotype)
% Returns 4x3 fiducials for canonical or anatomical torso
% Rows: [L_shoulder; R_shoulder; chin; lumbar]

if nargin<1, torsotype='anatomical'; end

switch lower(torsotype)
    case 'anatomical'
        % load anatomical fiducials
        fids = [
    30.114 141.137 127.07 
    39.799 142.93 -125.9 
    -49.1 188.98 -12.184 
    93.691 -163.67 -7.542 
];
    case 'canonical'
        % load canonical fiducials
        fids = [
  -47.1279  -13.0101   40.7700
   76.5116 -132.5392   95.2911
   80.2956 -170.6540   55.9314
  103.4446  -52.8995  -16.8010
        ];
    otherwise
        error('Unknown torsotype: %s', torsotype);
end

end
