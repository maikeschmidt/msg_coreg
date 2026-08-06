% cr_diagnose_coreg - Diagnose implausible spread in repeated coregistrations
%
% Recomputes, per repeat, the quantities that drive cr_register_torso, so
% that a large spread reported by cr_summarise_coreg can be attributed to a
% specific cause rather than assumed to be manual imprecision.
%
% RUN THIS WHENEVER cr_summarise_coreg REPORTS SOMETHING IMPLAUSIBLE
%   The giveaway is an internal inconsistency: small translation and
%   rotation spread but very large geometry displacement. Translation SD of
%   ~10 mm and rotation SD of ~2 deg can only move a 100 mm probe by of
%   order 15 mm. Displacements of hundreds of millimetres must come from
%   the SCALE term, and scale should be near-constant across repeats of the
%   same subject.
%
% THE MOST LIKELY CAUSE
%   cr_register_torso derives its unit-normalisation scale factor from the
%   AREA OF THE FIDUCIAL TRIANGLE:
%
%       pow = round(log10(sqrt(body_area / thorax_area)));
%       sf  = 10^pow;
%
%   This is a step function of the fiducial geometry. The shoulders are
%   broad, poorly defined landmarks and are typically the least repeatable
%   picks, so the triangle area varies between repeats. If log10(...) sits
%   near a .5 boundary, a few millimetres of shoulder wobble flips sf by a
%   FACTOR OF TEN, and that repeat's model ends up an order of magnitude
%   too large or too small.
%
%   This report prints the pre-rounding value and its distance from the
%   nearest boundary, so the failure mode is visible rather than inferred.
%
% USAGE:
%   cr_diagnose_coreg('coreg_repeats.mat')
%   D = cr_diagnose_coreg(R)
%
% INPUT:
%   src - path to the .mat saved by cr_repeat_coreg (containing R), or the
%         R struct itself
%
% OUTPUT:
%   D - struct with per-repeat fields:
%       .sf              scale factor actually applied
%       .log10_val       the value passed to round() before rounding
%       .margin          distance from the nearest rounding boundary
%                        (small = fragile; 0.5 = maximally safe)
%       .body_area       fiducial triangle cross-product magnitude, subject
%       .thorax_area     same, canonical torso
%       .total_scale     column norms of the composite transform
%       .translation     translation component
%       .suspect         logical, this repeat looks unsafe
%
% WHAT TO DO WITH THE RESULT
%   - If sf is NOT the same for every repeat, that is the bug. The affected
%     repeats are unusable and must be recollected, or the scale factor
%     pinned (see the printed guidance).
%   - If sf is constant but total_scale still varies, the problem is
%     downstream in the ICP step, not in unit normalisation.
%   - If margin is small (< 0.15) even where sf agrees, the pipeline is one
%     bad pick away from a 10x error and the scale should be pinned anyway.
%
% SEE ALSO:
%   cr_repeat_coreg, cr_summarise_coreg, cr_register_torso
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of the MSG Coregistration Toolbox.

function D = cr_diagnose_coreg(src)

if ischar(src) || isstring(src)
    tmp = load(char(src), 'R');
    R   = tmp.R;
else
    R = src;
end

n = R.n_repeats;

% Canonical torso fiducials, resolved exactly as cr_register_torso does
torso_file = fullfile(coreg_path, 'meshes', 'canonical_torso.stl');
stl_data   = stlread(torso_file);
if isa(stl_data, 'triangulation')
    V = stl_data.Points;
elseif isfield(stl_data, 'vertices')
    V = stl_data.vertices;
else
    error('Unsupported STL format returned for %s', torso_file);
end

torso_fid_ref = [ ...
     2.0898,   0.2269,  -0.4374;   % Left shoulder
    -1.8253,   0.0529,  -0.2542;   % Right shoulder
     0.1836,   1.1304,   1.2111];  % Chin

fid_idx = zeros(3,1);
for f = 1:3
    [~, fid_idx(f)] = min(vecnorm(V - torso_fid_ref(f,:), 2, 2));
end
torso_fids = V(fid_idx, :);

thorax_vec  = torso_fids([1 2],:) - torso_fids(3,:);
thorax_area = norm(cross(thorax_vec(1,:), thorax_vec(2,:)));

D = struct('sf', nan(n,1), 'log10_val', nan(n,1), 'margin', nan(n,1), ...
           'body_area', nan(n,1), 'thorax_area', thorax_area, ...
           'total_scale', nan(n,3), 'translation', nan(n,3), ...
           'suspect', false(n,1));

for k = 1:n
    fids = R.fiducials{k};

    body_vec  = fids([1 2],:) - fids(3,:);
    body_area = norm(cross(body_vec(1,:), body_vec(2,:)));

    val = log10(sqrt(body_area / thorax_area));
    pow = round(val);

    D.body_area(k)   = body_area;
    D.log10_val(k)   = val;
    D.sf(k)          = 10^pow;
    D.margin(k)      = abs(val - pow);   % 0 = on a boundary is 0.5 away
    D.margin(k)      = 0.5 - D.margin(k);

    T = R.transforms{k};
    D.total_scale(k,:) = vecnorm(T(1:3,1:3), 2, 1);
    D.translation(k,:) = T(1:3,4)';
end

% A repeat is suspect if its sf disagrees with the modal sf, or if it sits
% close to a rounding boundary
modal_sf     = mode(D.sf);
D.suspect    = (D.sf ~= modal_sf) | (D.margin < 0.15);


% REPORT

fprintf('\n=== Coregistration diagnostic (n = %d) ===\n\n', n);
fprintf('Canonical fiducial triangle magnitude : %.4f\n', thorax_area);
fprintf('Scale factor rule : sf = 10^round(log10(sqrt(body/thorax)))\n');
fprintf('  -> a STEP function; log10 near x.5 makes sf flip by 10x\n\n');

fprintf('%4s %12s %11s %9s %8s %10s %10s %10s  %s\n', ...
    'rep', 'body_area', 'log10_val', 'margin', 'sf', ...
    'scaleX', 'scaleY', 'scaleZ', 'status');
fprintf('%s\n', repmat('-', 1, 98));
for k = 1:n
    fprintf('%4d %12.1f %11.4f %9.4f %8g %10.4f %10.4f %10.4f  %s\n', ...
        k, D.body_area(k), D.log10_val(k), D.margin(k), D.sf(k), ...
        D.total_scale(k,1), D.total_scale(k,2), D.total_scale(k,3), ...
        ternary_str(D.suspect(k), '<<< SUSPECT', 'ok'));
end

fprintf('\n');

n_distinct = numel(unique(D.sf));
if n_distinct > 1
    fprintf('*** SCALE FACTOR IS NOT CONSTANT ACROSS REPEATS ***\n');
    fprintf('    Distinct values: %s\n', mat2str(unique(D.sf)'));
    fprintf('    This is the cause of the large geometry displacement.\n');
    fprintf('    Repeats whose sf differs from the modal value (%g) are\n', modal_sf);
    fprintf('    scaled by the wrong order of magnitude and MUST NOT be used.\n\n');
    fprintf('    FIX: pin the scale factor so it cannot flip. In\n');
    fprintf('    cr_register_torso, replace the call to\n');
    fprintf('    determine_body_scan_units with a fixed value, e.g.\n\n');
    fprintf('        sf = %g;   %% pinned; was derived from fiducial area\n\n', modal_sf);
    fprintf('    then recollect the affected repeats. The correct value is\n');
    fprintf('    the one that matches your scan units, not necessarily the\n');
    fprintf('    modal one - check that the registration figure overlays.\n');
elseif any(D.margin < 0.15)
    fprintf('Scale factor is consistent (sf = %g) but FRAGILE:\n', modal_sf);
    fprintf('  %d repeat(s) sit within 0.15 of a rounding boundary.\n', ...
        sum(D.margin < 0.15));
    fprintf('  One more bad shoulder pick would flip sf by 10x.\n');
    fprintf('  Consider pinning sf = %g in cr_register_torso.\n', modal_sf);
else
    fprintf('Scale factor is consistent (sf = %g) and robust.\n', modal_sf);
    fprintf('Unit normalisation is NOT the cause of the spread.\n');
    fprintf('Check the ICP step (M2) in cr_register_torso instead:\n');
    fprintf('  - S.dist too large, pulling in far-away vertices\n');
    fprintf('  - subject mesh containing background or floor geometry\n');
end

fprintf('\nTranslation per repeat (mm):\n');
for k = 1:n
    fprintf('  %2d: [%9.2f %9.2f %9.2f]\n', k, D.translation(k,:));
end
fprintf('\n');

end


function s = ternary_str(c, a, b)
    if c, s = a; else, s = b; end
end
