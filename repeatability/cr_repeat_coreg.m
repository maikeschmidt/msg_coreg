% cr_repeat_coreg - Repeat manual canonical coregistration N times and save
%                   every transform
%
% Prompts you to select the three torso fiducials (left shoulder, right
% shoulder, chin) N separate times on the SAME subject surface, running
% the full cr_register_torso pipeline after each selection, and saves
% every resulting transform matrix together with the fiducials that
% produced it.
%
% WHY THIS EXISTS
%   Reviewer 1 asked for 3-5 participants to assess variability. A single
%   participant was scanned, so inter-subject variability cannot be
%   measured directly. What CAN be measured, and is arguably the more
%   relevant quantity for a canonical-model pipeline, is how much the
%   forward solution moves when the SAME canonical model is coregistered
%   by hand several times. Each repeat is a legitimate independent
%   realisation of the coregistration procedure, so BEM and FEM can be
%   compared across them and the BEM-FEM agreement reported with a
%   distribution rather than a single number.
%
%   Be explicit in the manuscript that these are coregistration repeats,
%   NOT independent participants. They bound coregistration error, not
%   anatomical variability. Anatomical variability is addressed separately
%   by the warped models (see ../warping/cr_generate_warps.m).
%
% USAGE:
%   S.subject   = subject_surface_mesh;
%   S.n_repeats = 5;
%   S.outfile   = 'coreg_repeats.mat';
%   R = cr_repeat_coreg(S);
%
% INPUT:
%   S - struct with fields:
%
%   Required:
%     S.subject    - Subject surface mesh struct (.vertices, .faces)
%
%   Optional:
%     S.n_repeats  - Number of manual repeats (default: 5)
%     S.outfile    - Path to save results (default: 'coreg_repeats.mat')
%     S.plot       - Show the registration figure each repeat
%                    (default: false — set true to eyeball each fit, but
%                    it produces N figures)
%     S.dist       - ICP distance threshold passed to cr_register_torso
%                    (default: 0.02)
%     S.resume     - If true and S.outfile exists, load it and ADD the
%                    requested repeats to what is already there
%                    (default: true). Lets you do 2 today and 3 tomorrow.
%
% OUTPUT:
%   R - struct with fields:
%     .transforms  - {1 x N} cell array of 4x4 transform matrices
%     .fiducials   - {1 x N} cell array of 3x3 selected fiducial coords
%     .timestamps  - {1 x N} cell array of selection timestamps
%     .n_repeats   - number of repeats stored
%   Also saved to S.outfile.
%
% NEXT STEPS:
%   1. cr_summarise_coreg     — quantify the spread across repeats
%   2. cr_build_coreg_geometries — build one geometry .mat per repeat
%   3. run BEM and FEM in msg_fwd on each geometry
%   4. msg_fwd/stats/st_collect_replicates — group-level statistics
%
% NOTES:
%   - Use torso_mode 'canonical' semantics: this registers the canonical
%     torso, which is the model the reviewer's point concerns.
%   - Select fiducials as consistently as you can, exactly as you would
%     in real use. Deliberately sloppy selection would inflate the spread
%     and misrepresent the pipeline.
%   - Each repeat is independent: do NOT copy a previous selection.
%
% SEE ALSO:
%   cr_register_torso, cr_check_registration, cr_summarise_coreg,
%   cr_build_coreg_geometries
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

function R = cr_repeat_coreg(S)

if ~isfield(S, 'subject'),   error('Please provide S.subject (surface mesh).'); end
if ~isfield(S, 'n_repeats'), S.n_repeats = 5;                     end
if ~isfield(S, 'outfile'),   S.outfile   = 'coreg_repeats.mat';   end
if ~isfield(S, 'plot'),      S.plot      = false;                 end
if ~isfield(S, 'dist'),      S.dist      = 0.02;                  end
if ~isfield(S, 'resume'),    S.resume    = true;                  end

% Resume from an existing file if asked
R = struct('transforms', {{}}, 'fiducials', {{}}, 'timestamps', {{}}, ...
           'n_repeats', 0);

if S.resume && isfile(S.outfile)
    prev = load(S.outfile, 'R');
    if isfield(prev, 'R')
        R = prev.R;
        fprintf('Resuming: %d repeat(s) already stored in %s\n', ...
            R.n_repeats, S.outfile);
    end
end

start_idx = R.n_repeats + 1;
end_idx   = R.n_repeats + S.n_repeats;

fprintf('\n=== Manual coregistration repeatability ===\n');
fprintf('Collecting repeats %d to %d.\n', start_idx, end_idx);
fprintf('For each repeat select THREE fiducials in this order:\n');
fprintf('   1. Left shoulder\n   2. Right shoulder\n   3. Chin\n\n');

for k = start_idx:end_idx

    fprintf('--- Repeat %d of %d ---\n', k, end_idx);
    fprintf('Select the three fiducials now...\n');

    fids_sel = spm_mesh_select(S.subject);
    fids     = fids_sel';

    if ~isequal(size(fids), [3 3])
        warning(['Expected 3 fiducials (3x3), got %dx%d. ' ...
                 'Repeat %d skipped — rerun to collect it.'], ...
                 size(fids,1), size(fids,2), k);
        continue;
    end

    regS           = struct();
    regS.subject   = S.subject;
    regS.fiducials = fids;
    regS.dist      = S.dist;
    regS.plot      = S.plot;

    T = cr_register_torso(regS);

    R.transforms{end+1} = T;
    R.fiducials{end+1}  = fids;
    R.timestamps{end+1} = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    R.n_repeats         = numel(R.transforms);

    % Save after EVERY repeat so an interrupted session loses nothing
    save(S.outfile, 'R');
    fprintf('Repeat %d stored. Saved to %s\n\n', k, S.outfile);
end

fprintf('=== Complete: %d repeat(s) stored ===\n', R.n_repeats);
fprintf('Next: cr_summarise_coreg(''%s'')\n', S.outfile);

end
