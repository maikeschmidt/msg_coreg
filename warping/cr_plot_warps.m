% cr_plot_warps - Visual sanity check on the generated anatomical warps
%
% Draws the original torso outline against a selection of warped torso
% outlines, plus a scatter of the scale factors, so you can confirm the
% warps are plausible body shapes and not degenerate before spending hours
% of FEM compute on them.
%
% ALWAYS RUN THIS before cr_build_warp_geometries. A warp that turns the
% torso inside out or collapses an axis will still produce a geometry file
% and will still run through the BEM, silently poisoning the group
% statistics.
%
% USAGE:
%   cr_plot_warps(geom_file)
%   cr_plot_warps(geom_file, warp_file)
%   cr_plot_warps(geom_file, W, n_show)
%
% INPUT:
%   geom_file - Path to the base geometry .mat (or a struct with
%               mesh_torso and mesh_wm)
%   warp_file - Path to the warp .mat from cr_generate_warps, or the W
%               struct directly (default: 'anatomical_warps.mat')
%   n_show    - How many warps to overlay (default: 8, evenly spaced
%               through the set)
%
% WHAT TO LOOK FOR:
%   - Warped outlines should look like plausible bodies: taller/narrower
%     or shorter/wider, all still recognisably a torso
%   - Each cord should sit inside the torso drawn in the SAME COLOUR
%   - Scale factors should cluster around 1 with no axis collapsed
%   - determinant should be ~1 for every warp when volume preservation is on
%
% WHY THE CORD CANNOT ESCAPE THE TORSO
%   The warp is a single invertible affine map applied to every mesh, the
%   source positions and the sensors alike. Affine maps preserve
%   containment exactly: if the cord was inside the torso before warping it
%   is inside afterwards, for any warp, with no tolerance argument needed.
%   The containment test below is therefore a guard against a coding error
%   (a mesh accidentally left untransformed), not against the warping
%   itself. If it ever fails, look for an un-warped field in
%   cr_build_warp_geometries rather than for a bad warp.
%
% SEE ALSO:
%   cr_generate_warps, cr_build_warp_geometries
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

function cr_plot_warps(geom_file, warp_src, n_show)

if nargin < 2 || isempty(warp_src), warp_src = 'anatomical_warps.mat'; end
if nargin < 3 || isempty(n_show),   n_show   = 8;                     end

% Load geometry
if ischar(geom_file) || isstring(geom_file)
    geom = load(char(geom_file));
else
    geom = geom_file;
end

% Load warps
if ischar(warp_src) || isstring(warp_src)
    tmp = load(char(warp_src), 'W');
    W   = tmp.W;
else
    W = warp_src;
end

n = numel(W.matrices);

% Recentre about the torso centroid, matching cr_build_warp_geometries
centre = mean(geom.mesh_torso.vertices, 1);
P          = W.params;
P.centre   = centre;
P.n_warps  = n;
P.outfile  = fullfile(tempdir, 'anatomical_warps_plotcheck.mat');
evalc('W = cr_generate_warps(P);');   % suppress the printed table

% Pick evenly spaced warps to display
show_idx = unique(round(linspace(1, n, min(n_show, n))));

V0 = geom.mesh_torso.vertices;
C0 = geom.mesh_wm.vertices;

figure('Color','w','Name','Anatomical warp check','Position',[80 80 1400 620]);

% The cord is warped by the SAME matrix as the torso, because
% cr_build_warp_geometries applies M to every mesh, to sources_cent.pos and
% to the sensors. Both must therefore be drawn warped, in matching colours.
%
% Drawing the UNWARPED cord against a WARPED torso — as an earlier version
% of this figure did — makes the cord appear to escape the body on any warp
% that shortens the torso, which is a drawing error, not a geometry error.

cmap = parula(numel(show_idx));

% PANEL 1: sagittal outlines (Y vs Z)
subplot(1,3,1); hold on;
plot_outline(V0(:,3), V0(:,2), [0 0 0], 2.5);
plot(C0(:,3), C0(:,2), '.', 'Color', [0 0 0], 'MarkerSize', 2);
for i = 1:numel(show_idx)
    M = W.matrices{show_idx(i)};
    V = apply_T(M, V0);
    C = apply_T(M, C0);
    plot_outline(V(:,3), V(:,2), cmap(i,:), 1.0);
    plot(C(:,3), C(:,2), '.', 'Color', cmap(i,:), 'MarkerSize', 2);
end
axis equal; grid on;
xlabel('Ventral-Dorsal (mm)'); ylabel('Rostral-Caudal (mm)');
title(sprintf('Sagittal — %d of %d warps\n(cord warped with its own torso)', ...
    numel(show_idx), n));
set(gca,'FontSize',11,'TickDir','out');

% PANEL 2: coronal outlines (X vs Y)
subplot(1,3,2); hold on;
plot_outline(V0(:,1), V0(:,2), [0 0 0], 2.5);
plot(C0(:,1), C0(:,2), '.', 'Color', [0 0 0], 'MarkerSize', 2);
for i = 1:numel(show_idx)
    M = W.matrices{show_idx(i)};
    V = apply_T(M, V0);
    C = apply_T(M, C0);
    plot_outline(V(:,1), V(:,2), cmap(i,:), 1.0);
    plot(C(:,1), C(:,2), '.', 'Color', cmap(i,:), 'MarkerSize', 2);
end
axis equal; grid on;
xlabel('Left-Right (mm)'); ylabel('Rostral-Caudal (mm)');
title('Coronal (black = original torso and cord)');
set(gca,'FontSize',11,'TickDir','out');

% PANEL 3: scale factors and determinants
subplot(1,3,3); hold on;
plot(1:n, W.scales(:,1), 'o-', 'DisplayName', 'X (Left-Right)');
plot(1:n, W.scales(:,2), 's-', 'DisplayName', 'Y (Rostral-Caudal)');
plot(1:n, W.scales(:,3), '^-', 'DisplayName', 'Z (Ventral-Dorsal)');
yline(1, '--k', 'Alpha', 0.5, 'HandleVisibility','off');
grid on; xlabel('Warp index'); ylabel('Scale factor');
title('Per-axis scale factors'); legend('Location','best');
set(gca,'FontSize',11,'TickDir','out');

% NUMERIC CHECKS

fprintf('\n=== Warp sanity check ===\n');
dets = cellfun(@(M) det(M(1:3,1:3)), W.matrices);
fprintf('Determinant : min %.4f  max %.4f  (1.0 = volume preserved)\n', ...
    min(dets), max(dets));

if any(dets <= 0)
    fprintf('*** WARNING: %d warp(s) have non-positive determinant ***\n', sum(dets<=0));
    fprintf('    These flip handedness and MUST NOT be used.\n');
end

min_scale = min(W.scales(:));
fprintf('Scale range : %.3f to %.3f\n', min_scale, max(W.scales(:)));
if min_scale < 0.5
    fprintf('*** WARNING: a scale factor below 0.5 will badly distort anatomy ***\n');
end

% Cord containment check across all warps
fprintf('Cord containment: ');
bad = 0;
for k = 1:n
    Vk = apply_T(W.matrices{k}, V0);
    Ck = apply_T(W.matrices{k}, C0);
    inside = all(Ck >= min(Vk) - 1e-6, 2) & all(Ck <= max(Vk) + 1e-6, 2);
    if ~all(inside), bad = bad + 1; end
end
if bad == 0
    fprintf('OK — cord within torso bounding box in all %d warps\n\n', n);
else
    fprintf('*** %d warp(s) put cord outside the torso bounding box ***\n\n', bad);
end

end


% LOCAL FUNCTIONS

function p = apply_T(T, pts)
    p = (T * [pts, ones(size(pts,1),1)]')';
    p = p(:, 1:3);
end

function plot_outline(x, y, col, lw)
% Convex-hull outline of a projected point cloud — a cheap silhouette.
    k = convhull(x, y);
    plot(x(k), y(k), '-', 'Color', col, 'LineWidth', lw);
end
