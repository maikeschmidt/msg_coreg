% cr_summarise_coreg - Quantify the spread across repeated manual
%                      coregistrations
%
% Decomposes each saved 4x4 transform into interpretable components
% (translation, rotation, scale) and reports how much they vary across
% repeats. Produces the numbers you need to state in the manuscript when
% describing what the coregistration repeats actually bound.
%
% USAGE:
%   Stats = cr_summarise_coreg(matfile)
%   Stats = cr_summarise_coreg(matfile, do_plot)
%   Stats = cr_summarise_coreg(R)
%
% INPUT:
%   matfile - path to the .mat saved by cr_repeat_coreg (containing R),
%             or the R struct itself
%   do_plot - logical, draw summary figure (default: true)
%
% OUTPUT:
%   Stats - struct with fields:
%     .n              number of repeats
%     .translation    [n x 3] translation component of each transform (mm)
%     .scale          [n x 3] per-axis scale factor of each transform
%     .euler_deg      [n x 3] ZYX Euler angles of the rotation part (deg)
%     .fid_coords     [n x 3 x 3] selected fiducials per repeat
%     .fid_spread_mm  [3 x 1] SD of each fiducial's position across repeats
%     .fid_rms_mm     scalar RMS fiducial displacement from the mean
%                     selection — the single headline number for
%                     "how repeatable was the manual selection?"
%     .translation_sd [1 x 3] SD of translation across repeats (mm)
%     .scale_sd       [1 x 3] SD of per-axis scale across repeats
%     .euler_sd_deg   [1 x 3] SD of Euler angles across repeats (deg)
%     .max_vertex_disp_mm  worst-case displacement, across repeats, of a
%                     unit-cube corner under each transform relative to
%                     the mean transform. A geometry-level summary of how
%                     far the models actually move.
%
% INTERPRETATION FOR THE MANUSCRIPT:
%   fid_rms_mm answers "how precisely can a user place the fiducials?"
%   translation_sd / euler_sd_deg answer "how much does that move the
%   model?" Report both. If the BEM-FEM agreement is stable across
%   repeats despite a non-trivial fid_rms_mm, that is the result: the
%   BEM-FEM comparison is robust to realistic coregistration error.
%
% SEE ALSO:
%   cr_repeat_coreg, cr_build_coreg_geometries
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
%
% This file is part of the MSG Coregistration Toolbox.

function Stats = cr_summarise_coreg(src, do_plot)

if nargin < 2, do_plot = true; end

% Accept a path or the struct directly
if ischar(src) || isstring(src)
    tmp = load(char(src), 'R');
    R   = tmp.R;
else
    R = src;
end

n = R.n_repeats;
if n < 2
    error('Need at least 2 repeats to summarise variability (have %d).', n);
end

Stats.n           = n;
Stats.translation = zeros(n, 3);
Stats.scale       = zeros(n, 3);
Stats.euler_deg   = zeros(n, 3);
Stats.fid_coords  = zeros(n, 3, 3);

for k = 1:n
    T = R.transforms{k};
    A = T(1:3, 1:3);

    Stats.translation(k, :) = T(1:3, 4)';

    % Per-axis scale = column norms of the linear part
    sc = vecnorm(A, 2, 1);
    Stats.scale(k, :) = sc;

    % Rotation = linear part with scale divided out
    Rot = A ./ sc;
    Stats.euler_deg(k, :) = rad2deg(rotm_to_euler_zyx(Rot));

    Stats.fid_coords(k, :, :) = R.fiducials{k};
end

% Fiducial selection repeatability
mean_fids = squeeze(mean(Stats.fid_coords, 1));      % [3 x 3]
dev       = zeros(n, 3);
for k = 1:n
    f        = squeeze(Stats.fid_coords(k, :, :));   % [3 x 3]
    dev(k,:) = vecnorm(f - mean_fids, 2, 2)';        % distance per fiducial
end
Stats.fid_spread_mm = std(dev, 0, 1)';
Stats.fid_rms_mm    = sqrt(mean(dev(:).^2));

Stats.translation_sd = std(Stats.translation, 0, 1);
Stats.scale_sd       = std(Stats.scale,       0, 1);
Stats.euler_sd_deg   = std(Stats.euler_deg,   0, 1);

% Geometry-level displacement: how far do the corners of a 200 mm cube
% centred on the model move, repeat vs mean transform?
probe = 100 * [ 1  1  1; 1  1 -1; 1 -1  1; 1 -1 -1;
               -1  1  1; -1  1 -1; -1 -1  1; -1 -1 -1];
T_mean = mean(cat(3, R.transforms{:}), 3);
disp_mm = zeros(n, 1);
for k = 1:n
    p1 = apply_T(R.transforms{k}, probe);
    p0 = apply_T(T_mean,          probe);
    disp_mm(k) = max(vecnorm(p1 - p0, 2, 2));
end
Stats.max_vertex_disp_mm = max(disp_mm);
Stats.vertex_disp_mm     = disp_mm;

% REPORT

fprintf('\n=== Coregistration repeatability (n = %d) ===\n\n', n);
fprintf('Fiducial selection:\n');
fprintf('  RMS displacement from mean selection : %6.2f mm\n', Stats.fid_rms_mm);
fprintf('  SD per fiducial (LS, RS, chin)       : %6.2f %6.2f %6.2f mm\n', ...
    Stats.fid_spread_mm);
fprintf('\nResulting transform spread (SD across repeats):\n');
fprintf('  Translation X/Y/Z : %6.2f %6.2f %6.2f mm\n', Stats.translation_sd);
fprintf('  Rotation Z/Y/X    : %6.2f %6.2f %6.2f deg\n', Stats.euler_sd_deg);
fprintf('  Scale X/Y/Z       : %6.4f %6.4f %6.4f\n',     Stats.scale_sd);
fprintf('\nGeometry displacement:\n');
fprintf('  Worst-case corner displacement vs mean transform : %6.2f mm\n', ...
    Stats.max_vertex_disp_mm);
fprintf('  Median across repeats                            : %6.2f mm\n', ...
    median(disp_mm));
fprintf('\n');

if do_plot
    figure('Color', 'w', 'Name', 'Coregistration repeatability', ...
           'Position', [100 100 1200 400]);

    subplot(1,3,1);
    bar(Stats.translation); grid on;
    xlabel('Repeat'); ylabel('Translation (mm)');
    legend({'X','Y','Z'}, 'Location','best'); title('Translation');
    set(gca,'FontSize',11,'TickDir','out');

    subplot(1,3,2);
    bar(Stats.euler_deg); grid on;
    xlabel('Repeat'); ylabel('Rotation (deg)');
    legend({'Z','Y','X'}, 'Location','best'); title('Rotation');
    set(gca,'FontSize',11,'TickDir','out');

    subplot(1,3,3);
    bar(disp_mm, 'FaceColor', [0.2 0.4 0.8]); grid on;
    xlabel('Repeat'); ylabel('Max corner displacement (mm)');
    title('Geometry displacement vs mean');
    set(gca,'FontSize',11,'TickDir','out');
end

end


% LOCAL FUNCTIONS

function p = apply_T(T, pts)
% Apply a 4x4 transform to [N x 3] points.
    p = (T * [pts, ones(size(pts,1),1)]')';
    p = p(:, 1:3);
end

function e = rotm_to_euler_zyx(Rot)
% ZYX Euler angles (radians) from a 3x3 rotation matrix.
% Returns [rz; ry; rx]. Handles the gimbal-lock case.
    sy = sqrt(Rot(1,1)^2 + Rot(2,1)^2);
    if sy > 1e-6
        rx = atan2(Rot(3,2), Rot(3,3));
        ry = atan2(-Rot(3,1), sy);
        rz = atan2(Rot(2,1), Rot(1,1));
    else
        rx = atan2(-Rot(2,3), Rot(2,2));
        ry = atan2(-Rot(3,1), sy);
        rz = 0;
    end
    e = [rz; ry; rx];
end
