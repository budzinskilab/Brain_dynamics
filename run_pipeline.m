function run_pipeline(steps)
% run_pipeline.m
%
% PURPOSE: Runs the whole pipeline (steps 1-6 from the README) and then
%          all the figure scripts, in order, with a log and a summary.
% ROLE:    Entry point for a complete run, e.g. overnight on the server.
%          Needs all the data folders (see README). Takes several hours.
%
% Usage (from the repository root):
%   run_pipeline                         % everything
%   run_pipeline({'figure_2_our_replication', 'figure_EDRLR_components'})
%
% Each script runs in its own workspace so variables don't leak between
% them. If a pipeline step fails, the run stops, since later steps would
% otherwise pick up old outputs. A failing figure script is logged and the
% run carries on. Output goes to run_pipeline_log.txt.

pipeline = {
    'connectome_harmonics'                          % 1: EDR+LR connectome + eigenmodes
    'step4_5_full_resolution_reconstruction'        % 2: task reconstruction, four bases
    'step5_static_reconstruction'                   % 3: parcel-level version + parcel distances
    'step6b_fc_longrange_reconstruction_fullres'    % 4: Fig. 2 long-range FC (~1 h)
    'generate_connectivity_matrices'                % 5: structural matrices
    'generate_functional_connectivity_matrices'     % 6: functional matrices (~30 min)
    };
figures = {
    'figure_2_our_replication'
    'generate_our_figure3A'
    'generate_our_figure3B'
    'figure_EDRLR_components'
    'figure_fc_longrange_mask'
    'figure_2_fMRI_reconstruction'
    'figure_3_tfMRI_reconstruction'
    'figure_3D_surface_reconstruction'
    };

if nargin < 1
    steps = [pipeline; figures];
    % step 6 only redraws if its saved results exist, so set them aside
    % to make it recompute (the old file is kept as *_previous.mat)
    fc_file = fullfile('results', 'functional_connectivity_matrices.mat');
    if exist(fc_file, 'file')
        movefile(fc_file, strrep(fc_file, '.mat', '_previous.mat'));
    end
end

log_file = 'run_pipeline_log.txt';
diary off; diary(log_file);
fprintf('run_pipeline started %s (MATLAB %s)\n', char(datetime('now')), version);

status = cell(numel(steps), 1);
minutes = nan(numel(steps), 1);
for i = 1:numel(steps)
    name = steps{i};
    fprintf('\n========== [%d/%d] %s  (%s) ==========\n', i, numel(steps), name, char(datetime('now')));
    t = tic;
    try
        run_one(name);
        status{i} = 'OK';
    catch err
        status{i} = ['FAILED: ' err.message];
        fprintf(2, '%s failed: %s\n', name, err.message);
    end
    minutes(i) = toc(t) / 60;
    fprintf('---------- %s: %s (%.1f min)\n', name, status{i}, minutes(i));

    if ~strcmp(status{i}, 'OK') && ismember(name, pipeline)
        fprintf(2, '\nStopping: later steps depend on %s.\n', name);
        status(i+1:end) = {'not run'};
        break
    end
end

fprintf('\n========== Summary (%s) ==========\n', char(datetime('now')));
for i = 1:numel(steps)
    if isnan(minutes(i))
        fprintf('%-45s %11s  %s\n', steps{i}, '', status{i});
    else
        fprintf('%-45s %7.1f min  %s\n', steps{i}, minutes(i), status{i});
    end
end
diary off;
end


function run_one(name)
% runs one script in this function's own (fresh) workspace
close all
run(name);
close all
end
