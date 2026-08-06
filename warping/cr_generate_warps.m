% cr_generate_warps - Generate N random anatomical warp matrices
%
% Produces a reproducible set of 4x4 affine warps that make the anatomical
% model taller and thinner, shorter and wider, and everything in between,
% with optional mild shear so the variants are not purely axis-aligned
% rescalings.
%
% WHY THIS EXISTS
%   Reviewer 1 asked for 3-5 participants; Reviewer 3 asked for the n = 1
%   anatomy to be listed as a limitation. Only one participant was
%   scanned. Warping the single anatomical model into N plausible body
%   shapes gives a population of DISTINCT GEOMETRIES on which BEM and FEM
%   can each be run and compared.
%
%   State clearly in the manuscript what this does and does not show.
%     It DOES show: the BEM-FEM agreement is not an artefact of one
%       particular set of vertebral proportions, and it quantifies how the
%       reported 35-72% LR amplitude effect varies with body shape.
%     It does NOT show: true inter-subject anatomical variability. The
%       warps are affine, so vertebra count, spinal curvature, relative
%       organ placement and all non-affine anatomical detail are preserved
%       from the single scanned participant. These are synthetic
%       geometries, not participants, and must never be described as n=30
%       subjects.
%
% VOLUME PRESERVATION
%   By default each warp is volume-preserving (scale factors are rescaled
%   so their product is 1). "Taller and thinner" therefore genuinely trades
%   height against width rather than just making the whole body bigger,
%   which is the anatomically meaningful axis of variation and avoids
%   confounding shape change with a pure size change.
%
% USAGE:
%   W = cr_generate_warps();
%   W = cr_generate_warps(S);
%
% INPUT:
%   S - (optional) struct with fields:
%     .n_warps        Number of warps (default: 30)
%     .seed           RNG seed for reproducibility (default: 20260806)
%     .scale_range    [lo hi] per-axis scale factor before volume
%                     normalisation (default: [0.85 1.15], i.e. +/-15%)
%     .shear_max      Maximum shear coefficient (default: 0.03).
%                     Set 0 for pure axis-aligned scaling.
%     .preserve_volume  (default: true)
%     .centre         [1 x 3] centre about which to warp, in the SAME
%                     space the warp will be applied (subject space, mm).
%                     Warping about the origin would translate the model
%                     as well as reshape it. Default [0 0 0]; set this to
%                     the model centroid for best behaviour.
%     .outfile        Path to save the warp set (default: 'anatomical_warps.mat')
%
% OUTPUT:
%   W - struct with fields:
%     .matrices    {1 x N} cell of 4x4 warp matrices
%     .scales      [N x 3] per-axis scale factors actually applied
%     .shears      [N x 3] shear coefficients applied
%     .labels      {1 x N} short descriptive labels, e.g. 'taller/thinner'
%     .params      the settings struct used
%   Also saved to S.outfile.
%
% AXIS CONVENTION (matches the rest of the toolbox):
%   X = Left-Right      (width)
%   Y = Rostral-Caudal  (height / along the cord)
%   Z = Ventral-Dorsal  (depth)
%   So scale = [0.9 1.15 0.9] is taller and narrower.
%
% NEXT STEPS:
%   cr_build_warp_geometries — apply each warp and write geometry files
%
% SEE ALSO:
%   cr_build_warp_geometries, cr_plot_warps, cr_check_registration
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

function W = cr_generate_warps(S)

if nargin < 1, S = struct(); end

if ~isfield(S, 'n_warps'),         S.n_warps         = 30;          end
if ~isfield(S, 'seed'),            S.seed            = 20260806;    end
if ~isfield(S, 'scale_range'),     S.scale_range     = [0.85 1.15]; end
if ~isfield(S, 'shear_max'),       S.shear_max       = 0.03;        end
if ~isfield(S, 'preserve_volume'), S.preserve_volume = true;        end
if ~isfield(S, 'centre'),          S.centre          = [0 0 0];     end
if ~isfield(S, 'outfile'),         S.outfile         = 'anatomical_warps.mat'; end

rng(S.seed, 'twister');

n  = S.n_warps;
lo = S.scale_range(1);
hi = S.scale_range(2);

W.matrices = cell(1, n);
W.scales   = zeros(n, 3);
W.shears   = zeros(n, 3);
W.labels   = cell(1, n);
W.params   = S;

c = S.centre(:);

% Translate to centre, warp, translate back
T_to   = eye(4); T_to(1:3, 4)   = -c;
T_back = eye(4); T_back(1:3, 4) =  c;

for k = 1:n

    % Per-axis scale factors
    sc = lo + (hi - lo) * rand(1, 3);

    if S.preserve_volume
        % Rescale so det = 1: pure shape change, no size change
        sc = sc / (prod(sc)^(1/3));
    end

    % Mild shear so warps are not purely axis-aligned
    sh = (2 * rand(1, 3) - 1) * S.shear_max;

    A = [sc(1),  sh(1),  sh(2);
         0,      sc(2),  sh(3);
         0,      0,      sc(3)];

    M            = eye(4);
    M(1:3, 1:3)  = A;

    W.matrices{k} = T_back * M * T_to;
    W.scales(k,:) = sc;
    W.shears(k,:) = sh;
    W.labels{k}   = describe_warp(sc);
end

save(S.outfile, 'W');

% REPORT

fprintf('\n=== Generated %d anatomical warps (seed %d) ===\n\n', n, S.seed);
fprintf('Scale range     : [%.2f, %.2f] per axis\n', lo, hi);
fprintf('Shear max       : %.3f\n', S.shear_max);
fprintf('Volume preserved: %d\n', S.preserve_volume);
fprintf('Warp centre     : [%.1f %.1f %.1f] mm\n\n', S.centre);
fprintf('%4s  %6s %6s %6s   %s\n', 'idx', 'sX(LR)', 'sY(RC)', 'sZ(VD)', 'shape');
fprintf('%s\n', repmat('-', 1, 52));
for k = 1:n
    fprintf('%4d  %6.3f %6.3f %6.3f   %s\n', k, W.scales(k,:), W.labels{k});
end
fprintf('\nSaved to: %s\n', S.outfile);
fprintf('Next: cr_build_warp_geometries\n\n');

end


% LOCAL FUNCTIONS

function lbl = describe_warp(sc)
% Human-readable shape label from the per-axis scale factors.
% X = width (LR), Y = height (RC), Z = depth (VD).
    height = sc(2);
    girth  = mean(sc([1 3]));

    if height > 1.03 && girth < 0.99
        lbl = 'taller / thinner';
    elseif height < 0.97 && girth > 1.01
        lbl = 'shorter / wider';
    elseif height > 1.03
        lbl = 'taller';
    elseif height < 0.97
        lbl = 'shorter';
    elseif girth > 1.03
        lbl = 'wider';
    elseif girth < 0.97
        lbl = 'thinner';
    else
        lbl = 'near-original';
    end
end
