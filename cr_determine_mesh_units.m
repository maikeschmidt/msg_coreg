function [units,sf] = cr_determine_mesh_units(meshes)
% Function from George - copied from tt_dertimine_mesh_units into here so
% no external tool boxes are needed

% try and guess what the units of a registered thorax mesh are based on
% surface area compared to its template state;

for ii = 1:length(meshes)
    A(ii) = spm_mesh_area(meshes{ii});
end
% thorax is assumed to be the mesh with the largest surface area
area_thorax = max(A);
area_orignal = 0.59; % 0.59 m^2 

log_ratio = abs(round(log10(sqrt(area_thorax/area_orignal))));

switch log_ratio
    case 0
        units = 'm';
        sf = 1;
    case 1
        units = 'dm';
        sf = 10;
    case 2
        units = 'cm';
        sf = 100;
    case 3
        units = 'mm';
        sf = 1000;
    otherwise
        error('thorax mesh units could not be identified!');
end