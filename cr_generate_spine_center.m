% cr_generate_spine_center - Generate a centre-line source model along the spine mesh
%
% Computes a series of centre points along the spine mesh by slicing at
% regular Y intervals, finding the cross-sectional midpoint and local
% diameter at each level. The result is a source struct compatible with
% downstream BEM/leadfield computations.
%
% USAGE:
%   src = cr_generate_spine_center(S)
%
% INPUT:
%   S                  - Structure with the following fields:
%
%   Required:
%     S.spine            - Spine mesh struct with .vertices and .faces
%     S.ylim             - [y_min, y_max] range along Y axis to sample (mm)
%     S.unit             - Units string, e.g. 'mm'
%
%   Optional:
%     S.resolution       - Initial slice spacing in mm (default: 1)
%     S.spacing_increase - Fractional increase in spacing per step, e.g.
%                          0.02 for 2% progressive spacing (default: 0)
%
% OUTPUT:
%   src           - Source model struct containing:
%                     .pos    - Centre-line positions [N x 3] (mm)
%                     .zdist  - Local cross-sectional diameter at each
%                               centre point [N x 1]
%                     .inside - Logical array of ones [N x 1]
%                     .unit   - Units string (inherited from S.unit)
%                     .zlim   - [min_z, max_z] range of collected Z values
%
% NOTES:
%   - Slicing proceeds from y_max down to y_min along the Y axis
%   - At each Y level, vertices within the current resolution window are
%     projected onto a local cross-sectional plane to estimate diameter
%     via convex hull; levels with fewer than 3 points use a nominal
%     diameter of 0.001
%   - If S.spacing_increase > 0, slice spacing grows progressively,
%     producing coarser sampling toward the lower spine
%
% EXAMPLE:
%   S.spine            = spine_mesh;
%   S.ylim             = [y_min, y_max];
%   S.unit             = 'mm';
%   S.resolution       = 1;
%   S.spacing_increase = 0.02;   % 2% progressive spacing
%   src = cr_generate_spine_center(S);
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


function src = cr_generate_spine_center(S)
    % Ensure required fields are present
    if ~isfield(S, 'resolution'); S.resolution = 1; end
    if ~isfield(S, 'spine'); error('Please provide the spine mesh!'); end
    if ~isfield(S, 'ylim'); error('Please provide the ylim (y_min, y_max)!'); end
    if ~isfield(S, 'unit'); error('Please provide the units!'); end
    if ~isfield(S, 'spacing_increase'); S.spacing_increase = 0; end  % e.g., 0.02 for 2%

    spine_mesh = S.spine;
    unit = S.unit;

    y_min = S.ylim(1);
    y_max = S.ylim(2);

    center_line = [];
    y_dist = [];
    z_collected = [];

    y = y_max;
    res = S.resolution;

    while y >= y_min
        points_at_y = spine_mesh.vertices(abs(spine_mesh.vertices(:,2) - y) < res, :);

        if ~isempty(points_at_y)
            x_mid = mean(points_at_y(:,1));
            z_mid = mean(points_at_y(:,3));

            if size(center_line, 1) > 1
                tangent = center_line(end, :) - center_line(end-1, :);
            else
                tangent = [0, 1, 0];
            end
            tangent = tangent / norm(tangent);

            up = [0, 0, 1]; 
            if abs(dot(up, tangent)) > 0.9
                up = [1, 0, 0]; 
            end
            normal1 = cross(tangent, up); normal1 = normal1 / norm(normal1);
            normal2 = cross(tangent, normal1);

            center_point = repmat([x_mid, y, z_mid], size(points_at_y, 1), 1);
            projected_points = [(points_at_y - center_point) * normal1', ...
                                (points_at_y - center_point) * normal2'];

            if size(projected_points, 1) > 2
                k = convhull(projected_points(:,1), projected_points(:,2));
                boundary_points = projected_points(k, :);
                radii = sqrt(sum(boundary_points.^2, 2));
                min_diam = 2 * mean(radii);
            else
                min_diam = 0.001;
            end

            center_line = [center_line; x_mid, y, z_mid];
            y_dist = [y_dist; min_diam];
            z_collected = [z_collected; z_mid];
        end

        % Increase the spacing dynamically
        res = res * (1 + S.spacing_increase);
        y = y - res;
    end

    src.pos = center_line;
    src.zdist = y_dist;
    src.inside = ones(size(src.pos, 1), 1);
    src.unit = unit;

    if ~isempty(z_collected)
        src.zlim = [min(z_collected), max(z_collected)];
    else
        src.zlim = [0, 0];
    end
end
