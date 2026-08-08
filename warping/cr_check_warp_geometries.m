function T = cr_check_warp_geometries(S)
% cr_check_warp_geometries - Run cr_check_geometry over a whole warp set
%
% Answers the question you actually need answered before starting a
% multi-day FEM run: WHICH of the 30 warps will fail to tetrahedralise, and
% why. Cheap — seconds per geometry — so run it before every FEM batch.
%
% If only a handful fail, drop them and report n honestly; the warping
% analysis does not need all 30. If most fail, the warp amplitudes are too
% aggressive for the base mesh quality and cr_repair_geometry (or a smaller
% scale_range in cr_generate_warps) is the fix.
%
% USAGE:
%   S.dir = 'D:\Simulations\Replicates\geometries';
%   T = cr_check_warp_geometries(S);
%
% INPUT:
%   S - struct:
%     .dir       (required) folder holding the geometry .mat files
%     .reference a geometry KNOWN to tetrahedralise (filename in .dir, or a
%                full path). Everything else is judged against it. Without
%                this you get absolute thresholds, which flag geometries
%                that mesh perfectly well.
%     .pattern   file mask (default 'geometries_*.mat')
%     .opts      options passed through to cr_check_geometry
%     .csv       path to write a summary CSV (default: <dir>/geometry_check.csv)
%
% OUTPUT:
%   T - struct array, one per file: name, ok, n_fatal, n_warn, min_angle,
%       worst_mesh, first_problem
%
% USAGE WITH A REFERENCE:
%   S.dir       = 'D:\...\geometries';
%   S.reference = 'geometries_coreg_rep01_cont.mat';   % this one meshed
%   T = cr_check_warp_geometries(S);
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of msg_coreg.

if ~isfield(S,'pattern'), S.pattern = 'geometries_*.mat'; end
if ~isfield(S,'opts'),    S.opts    = struct();           end
if ~isfield(S,'csv'),     S.csv     = fullfile(S.dir, 'geometry_check.csv'); end
S.opts.verbose = false;

d = dir(fullfile(S.dir, S.pattern));
if isempty(d)
    error('No files matching %s in %s', S.pattern, S.dir);
end

% REFERENCE GEOMETRY
% Absolute quality thresholds are guesswork — TetGen accepts far thinner
% triangles than intuition suggests. Point S.reference at a geometry you
% have ALREADY tetrahedralised successfully and every other geometry is
% judged against it, which is the only calibration that means anything.
if isfield(S,'reference') && ~isempty(S.reference)
    ref = S.reference;
    if ~isfile(ref), ref = fullfile(S.dir, ref); end
    fprintf('Reference (known to mesh): %s\n', ref);
    ro = S.opts; ro.verbose = false; ro.ref_mesh = []; ro.tolerated = {};
    Rref = cr_check_geometry(ref, ro);
    S.opts.ref_mesh = Rref.mesh;

    % Anything the reference ALSO fails cannot be what stops TetGen — this
    % geometry demonstrably meshes. Those categories are demoted to warnings
    % for every file, which is what keeps an invented threshold from
    % condemning geometries that work.
    S.opts.tolerated = unique(Rref.fatal_cat);
    if isempty(Rref.mesh)
        fprintf('  WARNING: no meshes found in the reference.\n');
    else
        fprintf('  reference worst angle : %.4f deg\n', min([Rref.mesh.min_angle]));
    end
    if isempty(S.opts.tolerated)
        fprintf('  reference is clean on every check.\n\n');
    else
        fprintf('  reference also fails  : %s\n', strjoin(S.opts.tolerated, ', '));
        fprintf('  -> those are demoted to warnings everywhere.\n\n');
    end
end

fprintf('=== Checking %d geometries ===\n\n', numel(d));
fprintf('%-42s %6s %8s %8s %10s %-16s\n', ...
    'file', 'ok', 'fatals', 'warns', 'min angle', 'worst mesh');
fprintf('%s\n', repmat('-', 1, 96));

T = struct('name',{},'ok',{},'n_fatal',{},'n_warn',{},'min_angle',{}, ...
           'worst_mesh',{},'first_problem',{});

for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);
    try
        R = cr_check_geometry(f, S.opts);
    catch err
        T(end+1) = struct('name', d(k).name, 'ok', false, 'n_fatal', NaN, ...
            'n_warn', NaN, 'min_angle', NaN, 'worst_mesh', '', ...
            'first_problem', ['check errored: ' err.message]); %#ok<AGROW>
        fprintf('%-42s %6s %8s %8s %10s %-16s\n', d(k).name, 'ERR', '-', '-', '-', '-');
        fprintf('%52s %s\n', '', err.message);
        if k == 1
            % A failure on the very first file is almost always a bug in the
            % checker or a wrong path, not 100 bad geometries. Show where.
            fprintf('\n  First file errored — stack:\n');
            for s = 1:min(3, numel(err.stack))
                fprintf('    %s line %d\n', err.stack(s).name, err.stack(s).line);
            end
            fprintf('\n');
        end
        continue;
    end

    if isempty(R.mesh)
        ma = NaN; wm = '';
    else
        [ma, wi] = min([R.mesh.min_angle]);
        wm = R.mesh(wi).name;
    end

    if isempty(R.fatal), fp = ''; else, fp = R.fatal{1}; end

    T(end+1) = struct('name', d(k).name, 'ok', R.ok, ...
        'n_fatal', numel(R.fatal), 'n_warn', numel(R.warnings), ...
        'min_angle', ma, 'worst_mesh', wm, 'first_problem', fp); %#ok<AGROW>

    fprintf('%-42s %6s %8d %8d %10.4f %-16s\n', d(k).name, ...
        yesno(R.ok), numel(R.fatal), numel(R.warnings), ma, wm);
    if ~R.ok
        fprintf('%52s %s\n', '', fp);
    elseif ~isempty(R.warnings)
        fprintf('%52s %s\n', '', R.warnings{1});
    end
end

n_ok = sum([T.ok]);
fprintf('\n%d of %d geometries pass (%d flagged).\n', ...
    n_ok, numel(T), numel(T) - n_ok);

if ~isfield(S,'reference') || isempty(S.reference)
    fprintf(['\nNo reference given, so these verdicts rest on absolute\n' ...
             'thresholds. Pass S.reference = a geometry you have already\n' ...
             'tetrahedralised and re-run — it is the only calibration that\n' ...
             'means anything.\n']);
elseif n_ok == numel(T)
    fprintf(['\nNothing here distinguishes these geometries from the one that\n' ...
             'meshed. If the FEM still fails, the cause is not surface\n' ...
             'quality — look at the compartment seed points in\n' ...
             'run_fem_leadfields, which are sampled per run and can land in\n' ...
             'the wrong region on a deformed geometry.\n']);
end

fid = fopen(S.csv, 'w');
fprintf(fid, 'file,ok,n_fatal,n_warn,min_angle_deg,worst_mesh,first_problem\n');
for k = 1:numel(T)
    fprintf(fid, '%s,%d,%g,%g,%.6f,%s,"%s"\n', T(k).name, T(k).ok, ...
        T(k).n_fatal, T(k).n_warn, T(k).min_angle, T(k).worst_mesh, ...
        T(k).first_problem);
end
fclose(fid);
fprintf('Summary written to %s\n', S.csv);

end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end
