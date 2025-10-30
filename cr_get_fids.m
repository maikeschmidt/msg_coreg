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
