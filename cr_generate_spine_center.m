function src = cr_generate_spine_center(S)
    % Ensure required fields are present
    if ~isfield(S, 'T'); error('Please provide the transformation matrix!'); end
    if ~isfield(S, 'resolution'); S.resolution = 1; end  % Default resolution
    if ~isfield(S, 'spine'); error('Please provide the spine mesh!'); end
    if ~isfield(S, 'zlim'); error('Please provide the zlim (z_min, z_max)!'); end
    if ~isfield(S, 'unit'); error('Please provide the units!'); end

    spine_mesh = S.spine;
    unit = S.unit;

    z_min = S.zlim(1);
    z_max = S.zlim(2);

    center_line = [];
    y_dist = [];

    z_values = z_max:-S.resolution:z_min;
    
    for z = z_values
        points_at_z = spine_mesh.vertices(abs(spine_mesh.vertices(:,3) - z) < S.resolution, :);
        
        if ~isempty(points_at_z)
            x_mid = mean(points_at_z(:,1));
            y_mid = mean(points_at_z(:,2));
            
            if size(center_line,1) > 1
                tangent = center_line(end, :) - center_line(end-1, :);
            else
                tangent = [0, 0, 1]; %
            end
            tangent = tangent / norm(tangent);

            % Find two orthonormal vectors to define the cross-section plane
            up = [0, 0, 1]; 
            if abs(dot(up, tangent)) > 0.9
                up = [1, 0, 0]; 
            end
            normal1 = cross(tangent, up);
            normal1 = normal1 / norm(normal1);
            normal2 = cross(tangent, normal1);

            center_point = repmat([x_mid, y_mid, z], size(points_at_z, 1), 1);
            projected_points = [(points_at_z - center_point) * normal1', ...
                                (points_at_z - center_point) * normal2'];
            
            if size(projected_points,1) > 2
                k = convhull(projected_points(:,1), projected_points(:,2));
                boundary_points = projected_points(k, :);
                radii = sqrt(sum(boundary_points.^2, 2));
                min_diam = 2 * mean(radii); % Averaging distances
            else
                min_diam = 0.001; % Default small diameter if too few points
            end

            center_line = [center_line; x_mid, y_mid, z];
            y_dist = [y_dist; min_diam]; 
        end
    end
    
    src.pos = center_line;
    src.ydist = y_dist;
    src.inside = ones(size(src.pos, 1), 1);
    src.unit = unit;
end
