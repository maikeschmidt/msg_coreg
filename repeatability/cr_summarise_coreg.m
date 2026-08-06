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

% Geometry-level displacement.
%
% Measured on the ACTUAL canonical torso vertices, not on an arbitrary
% probe cube. The transforms map canonical space (where the torso spans
% only about 8 units) to subject space, and they carry a scale factor of
% order 100. A probe sized in millimetres would therefore sit tens of
% metres from the origin once transformed, and a couple of degrees of
% rotation at that radius produces displacements of hundreds of
% millimetres that say nothing about the anatomy. Probing the real mesh
% keeps the number interpretable: it is how far the model itself moves.
%
% Computed PAIRWISE between repeats rather than against a mean transform.
% An elementwise mean of 4x4 matrices is not a valid mean of rigid
% transforms, and the pairwise form answers the question directly: how far
% apart are two independent coregistrations of the same subject?

probe = canonical_probe_points();

disp_pairwise = [];
for a = 1:n
    for b = a+1:n
        pa = apply_T(R.transforms{a}, probe);
        pb = apply_T(R.transforms{b}, probe);
        d  = vecnorm(pa - pb, 2, 2);
        disp_pairwise(end+1, :) = [max(d), median(d)]; %#ok<AGROW>
    end
end

Stats.pairwise_max_disp_mm    = disp_pairwise(:, 1);
Stats.pairwise_median_disp_mm = disp_pairwise(:, 2);
Stats.max_vertex_disp_mm      = max(disp_pairwise(:, 1));
Stats.median_vertex_disp_mm   = median(disp_pairwise(:, 2));

% Per-repeat displacement from the most central repeat, for the bar plot.
% The most central repeat is the one with the smallest total pairwise
% distance to all others — a genuine member of the set, unlike a mean.
tot = zeros(n, 1);
for a = 1:n
    for b = 1:n
        if a == b, continue; end
        pa = apply_T(R.transforms{a}, probe);
        pb = apply_T(R.transforms{b}, probe);
        tot(a) = tot(a) + median(vecnorm(pa - pb, 2, 2));
    end
end
[~, ref_rep] = min(tot);
Stats.reference_repeat = ref_rep;

disp_mm = zeros(n, 1);
p_ref   = apply_T(R.transforms{ref_rep}, probe);
for k = 1:n
    pk = apply_T(R.transforms{k}, probe);
    disp_mm(k) = median(vecnorm(pk - p_ref, 2, 2));
end
Stats.vertex_disp_mm = disp_mm;

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
fprintf('\nGeometry displacement (canonical torso vertices, pairwise):\n');
fprintf('  Median displacement between two repeats : %6.2f mm\n', ...
    Stats.median_vertex_disp_mm);
fprintf('  Worst-case vertex, worst pair           : %6.2f mm\n', ...
    Stats.max_vertex_disp_mm);
fprintf('  Most central repeat                     : %d\n', Stats.reference_repeat);
fprintf(['\n  This is the headline number: two independent manual\n' ...
         '  coregistrations of the same subject place the model about\n' ...
         '  %.0f mm apart. Report it alongside the fiducial RMS.\n'], ...
    Stats.median_vertex_disp_mm);
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
    xlabel('Repeat'); ylabel('Median vertex displacement (mm)');
    title(sprintf('Displacement vs repeat %d (most central)', ref_rep));
    set(gca,'FontSize',11,'TickDir','out');
end

end


% LOCAL FUNCTIONS

function p = apply_T(T, pts)
% Apply a 4x4 transform to [N x 3] points.
    p = (T * [pts, ones(size(pts,1),1)]')';
    p = p(:, 1:3);
end

function pts = canonical_probe_points()
% Vertices of the canonical torso, in canonical space — the same space the
% saved transforms take as input. Subsampled for speed; displacement
% statistics do not need every vertex.
%
% Falls back to the canonical bounding box if the STL cannot be read, so a
% missing mesh degrades the metric rather than erroring out of the whole
% summary.
    try
        stl_data = stlread(fullfile(coreg_path, 'meshes', 'canonical_torso.stl'));
        if isa(stl_data, 'triangulation')
            V = stl_data.Points;
        else
            V = stl_data.vertices;
        end
        step = max(1, round(size(V,1) / 500));
        pts  = V(1:step:end, :);
    catch
        warning('cr_summarise_coreg:noCanonicalMesh', ...
            ['Could not read canonical_torso.stl — falling back to its ' ...
             'bounding box. Displacement figures remain valid but coarser.']);
        lo = [-1.953, -4.874, -1.070];
        hi = [ 2.288,  3.139,  1.224];
        [x,y,z] = ndgrid([lo(1) hi(1)], [lo(2) hi(2)], [lo(3) hi(3)]);
        pts = [x(:), y(:), z(:)];
    end
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
