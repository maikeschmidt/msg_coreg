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
%     .pattern   file mask (default 'geometries_*.mat')
%     .opts      options passed through to cr_check_geometry
%     .csv       path to write a summary CSV (default: <dir>/geometry_check.csv)
%
% OUTPUT:
%   T - struct array, one per file: name, ok, n_fatal, min_angle, worst_mesh,
%       first_problem
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

fprintf('=== Checking %d geometries ===\n\n', numel(d));
fprintf('%-42s %6s %8s %10s %-16s\n', 'file', 'ok', 'fatals', 'min angle', 'worst mesh');
fprintf('%s\n', repmat('-', 1, 88));

T = struct('name',{},'ok',{},'n_fatal',{},'min_angle',{},'worst_mesh',{}, ...
           'first_problem',{});

for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);
    try
        R = cr_check_geometry(f, S.opts);
    catch err
        T(end+1) = struct('name', d(k).name, 'ok', false, 'n_fatal', NaN, ...
            'min_angle', NaN, 'worst_mesh', '', ...
            'first_problem', ['check errored: ' err.message]); %#ok<AGROW>
        fprintf('%-42s %6s %8s %10s %-16s\n', d(k).name, 'ERR', '-', '-', '-');
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
        'n_fatal', numel(R.fatal), 'min_angle', ma, 'worst_mesh', wm, ...
        'first_problem', fp); %#ok<AGROW>

    fprintf('%-42s %6s %8d %10.4f %-16s\n', d(k).name, ...
        yesno(R.ok), numel(R.fatal), ma, wm);
    if ~R.ok
        fprintf('%50s %s\n', '', fp);
    end
end

n_ok = sum([T.ok]);
fprintf('\n%d of %d geometries pass (%d would fail TetGen).\n', ...
    n_ok, numel(T), numel(T) - n_ok);

fid = fopen(S.csv, 'w');
fprintf(fid, 'file,ok,n_fatal,min_angle_deg,worst_mesh,first_problem\n');
for k = 1:numel(T)
    fprintf(fid, '%s,%d,%g,%.6f,%s,"%s"\n', T(k).name, T(k).ok, ...
        T(k).n_fatal, T(k).min_angle, T(k).worst_mesh, T(k).first_problem);
end
fclose(fid);
fprintf('Summary written to %s\n', S.csv);

end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end
