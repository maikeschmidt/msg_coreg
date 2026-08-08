function Scan = cr_scan_intersections(S)
% cr_scan_intersections - Run cr_find_intersections over a folder of geometries
%
% Answers, in one pass, which geometries will fail TetGen with a
% self-intersection and which meshes are responsible — before committing
% days of FEM time to finding out one geometry at a time.
%
% Each geometry is judged by TetGen itself (-d), so a PASS here means TetGen
% found no collision, not that some proxy check was satisfied.
%
% USAGE:
%   S.dir = 'D:\Simulations\...\geometries';
%   Scan  = cr_scan_intersections(S);
%
% INPUT:
%   S - struct:
%     .dir        (required) folder holding the geometry .mat files
%     .pattern    file mask (default 'geometries_*.mat')
%     .outfile    .mat to write (default <dir>/intersection_scan.mat)
%     .csv        .csv to write (default <dir>/intersection_scan.csv)
%     .resume     skip files already in an existing .outfile (default true)
%     .max_report collisions listed per geometry (default 40)
%     .split      split meshes into components, so a 23-vertebra bone mesh
%                 reports WHICH vertebra (default true)
%
% OUTPUT / SAVED
%   Scan - struct array, one entry per geometry:
%     .name        filename
%     .ok          true if TetGen found no self-intersection
%     .n_int       number of collisions reported
%     .top_pair    the most frequent colliding mesh pair
%     .pairs       every colliding pair, with counts
%     .I           the full cr_find_intersections output
%     .err         error message, if the scan itself failed
%
%   The .mat holds Scan. The .csv holds one row per geometry for a quick
%   look, and one row per pair in a second file (<csv stem>_pairs.csv) so
%   the recurring culprits are obvious.
%
% READING THE RESULT
%   Rigid and affine transforms cannot create a collision that was not
%   already present — they move every mesh together. So if a geometry family
%   fails uniformly, suspect the BASE model, not the coregistration or the
%   warp. The variant breakdown printed at the end is there to make that
%   visible: if every 'inhomo' fails and every 'realistic' passes, the
%   problem is the toroidal bone model, not the replicate generation.
%
% -------------------------------------------------------------------------
% Copyright (c) 2026 University College London
% Department of Imaging Neuroscience
%
% Author: Maike Schmidt
% Email:  maike.schmidt.23@ucl.ac.uk
%
% This file is part of msg_coreg.

if ~isfield(S,'pattern'),    S.pattern    = 'geometries_*.mat'; end
if ~isfield(S,'outfile'),    S.outfile    = fullfile(S.dir, 'intersection_scan.mat'); end
if ~isfield(S,'csv'),        S.csv        = fullfile(S.dir, 'intersection_scan.csv'); end
if ~isfield(S,'resume'),     S.resume     = true;  end
if ~isfield(S,'max_report'), S.max_report = 40;    end
if ~isfield(S,'split'),      S.split      = true;  end

d = dir(fullfile(S.dir, S.pattern));
if isempty(d)
    error('No files matching %s in %s', S.pattern, S.dir);
end

Scan = struct('name',{},'ok',{},'n_int',{},'top_pair',{},'pairs',{}, ...
              'I',{},'err',{});
done = {};

% Resume: TetGen -d is not free, and a 100-geometry scan is worth not
% repeating because MATLAB was closed halfway through.
if S.resume && isfile(S.outfile)
    old = load(S.outfile, 'Scan');
    if isfield(old, 'Scan') && ~isempty(old.Scan)
        Scan = old.Scan;
        done = {Scan.name};
        fprintf('Resuming: %d geometries already scanned.\n', numel(done));
    end
end

fprintf('\n=== Intersection scan: %d geometries ===\n\n', numel(d));
fprintf('%-44s %6s %8s  %s\n', 'file', 'ok', 'collis', 'worst pair');
fprintf('%s\n', repmat('-', 1, 100));

opts = struct('verbose', false, 'split', S.split, 'max_report', S.max_report);

for k = 1:numel(d)
    if any(strcmp(done, d(k).name))
        continue;
    end

    f = fullfile(d(k).folder, d(k).name);
    e = struct('name', d(k).name, 'ok', false, 'n_int', NaN, ...
               'top_pair', '', 'pairs', {{}}, 'I', [], 'err', '');

    try
        I = cr_find_intersections(f, opts);

        if isempty(I)
            e.ok = true; e.n_int = 0;
        else
            if isfield(I, 'n_total') && ~isempty(I(1).n_total)
                e.n_int = I(1).n_total;
            else
                e.n_int = numel(I);
            end
            e.I = I;

            lab = arrayfun(@(x) sprintf('%s x %s', ...
                pair_label(x.mesh_a, x.comp_a), ...
                pair_label(x.mesh_b, x.comp_b)), I, 'UniformOutput', false);
            [u, ~, ix] = unique(lab);
            cnt = accumarray(ix, 1);
            [cnt, o] = sort(cnt, 'descend');
            e.pairs    = arrayfun(@(j) sprintf('%s (%d)', u{o(j)}, cnt(j)), ...
                                  1:numel(u), 'UniformOutput', false);
            e.top_pair = u{o(1)};
        end
    catch err
        e.err = err.message;
    end

    Scan(end+1) = e; %#ok<AGROW>

    if ~isempty(e.err)
        fprintf('%-44s %6s %8s  %s\n', e.name, 'ERR', '-', e.err);
    else
        fprintf('%-44s %6s %8d  %s\n', e.name, yesno(e.ok), e.n_int, e.top_pair);
    end

    % Saved as we go, so an interrupted scan is not a wasted one
    save(S.outfile, 'Scan', '-v7.3');
end


% SUMMARY

ok_mask = arrayfun(@(x) isequal(x.ok, true), Scan);
fprintf('\n%d of %d geometries have NO self-intersection.\n', ...
    sum(ok_mask), numel(Scan));

failed = Scan(~ok_mask & cellfun(@isempty, {Scan.err}));
if ~isempty(failed)
    fprintf('\nRecurring culprits across all failing geometries:\n');
    allp = {};
    for k = 1:numel(failed)
        if isempty(failed(k).I), continue; end
        allp = [allp, arrayfun(@(x) sprintf('%s x %s', ...
            pair_label(x.mesh_a, x.comp_a), ...
            pair_label(x.mesh_b, x.comp_b)), failed(k).I, ...
            'UniformOutput', false)]; %#ok<AGROW>
    end
    if ~isempty(allp)
        [u, ~, ix] = unique(allp);
        cnt = accumarray(ix, 1);
        [cnt, o] = sort(cnt, 'descend');
        for k = 1:min(15, numel(u))
            fprintf('  %5d  %s\n', cnt(k), u{o(k)});
        end
    end
end

% Breakdown by bone variant — this is what separates "the replicates are
% broken" from "this bone model was never FEM-viable".
fprintf('\nBy bone variant:\n');
for v = {'cont','inhomo','homo','realistic'}
    sel = contains({Scan.name}, ['_' v{1}]);
    if ~any(sel), continue; end
    fprintf('  %-10s %3d of %3d clean\n', v{1}, ...
        sum(ok_mask(sel)), sum(sel));
end

fprintf('\nBy replicate family:\n');
for fam = {'coreg','warp'}
    sel = contains({Scan.name}, fam{1});
    if ~any(sel), continue; end
    fprintf('  %-10s %3d of %3d clean\n', fam{1}, ...
        sum(ok_mask(sel)), sum(sel));
end


% CSV

fid = fopen(S.csv, 'w');
fprintf(fid, 'file,ok,n_intersections,top_pair,error\n');
for k = 1:numel(Scan)
    fprintf(fid, '%s,%d,%g,"%s","%s"\n', Scan(k).name, Scan(k).ok, ...
        Scan(k).n_int, Scan(k).top_pair, Scan(k).err);
end
fclose(fid);

[p, n, ~] = fileparts(S.csv);
pcsv = fullfile(p, [n '_pairs.csv']);
fid = fopen(pcsv, 'w');
fprintf(fid, 'file,pair_and_count\n');
for k = 1:numel(Scan)
    for j = 1:numel(Scan(k).pairs)
        fprintf(fid, '%s,"%s"\n', Scan(k).name, Scan(k).pairs{j});
    end
end
fclose(fid);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n', S.outfile, S.csv, pcsv);
fprintf(['\nTo see every collision for one geometry:\n' ...
         '  cr_find_intersections(fullfile(S.dir, ''<name>.mat''));\n\n']);

end


% LOCAL FUNCTIONS

function s = pair_label(name, comp)
    if isempty(name), s = '?'; return; end
    if isnan(comp) || comp <= 1
        s = name;
    else
        s = sprintf('%s #%d', name, comp);
    end
end

function s = yesno(b)
    if b, s = 'yes'; else, s = 'NO'; end
end
