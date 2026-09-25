% generate_our_figure3B.m
%
% PURPOSE: Our Fig. 3B: per-task and task-averaged difference in
%          reconstruction error vs. Geometry over the first 20 modes.
% ROLE:    Figure (after step 2). Makes
%          results/our_figure3B_replication.png.
%
% Same statistic and plot as the authors' Fig. 3B code, applied to our
% step4_5 results so the two figures can be compared.

%% Our results and the redblue colormap

addpath(genpath(fullfile('Vohryzek2024', 'vohryzek2024_EDRLR', ...
    'Data', 'template_surfaces', 'colormaps_add'))); % redblue

load(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
    'mse_all', 'graph_names', 'task_names');

idx_geom  = find(strcmp(graph_names, 'Geometry'));
idx_bin   = find(strcmp(graph_names, 'EDR binary'));
idx_edr   = find(strcmp(graph_names, 'EDR continuous'));
idx_edrlr = find(strcmp(graph_names, 'EDR+LR'));

%% First 20 modes only, as in the paper (mse_all{g} is tasks x modes)
th_mode = 20;
mse_geom  = mse_all{idx_geom}(:, 1:th_mode);
mse_bin   = mse_all{idx_bin}(:, 1:th_mode);
mse_edr   = mse_all{idx_edr}(:, 1:th_mode);
mse_edrlr = mse_all{idx_edrlr}(:, 1:th_mode);

%% MSE minus Geometry's MSE (negative = better than Geometry)
diff_bin   = mse_bin   - mse_geom;
diff_edr   = mse_edr   - mse_geom;
diff_edrlr = mse_edrlr - mse_geom;

%% Task domains for the heatmap rows
% The paper labels rows by domain rather than by task. The domain comes
% from each task name's prefix, and the rows are already grouped by domain.
domain_prefix = {'social', 'motor', 'gambling', 'wm', 'language', 'emotion', 'relational'};
domain_label  = {'Social', 'Motor', 'Gambling', 'Working Memory', 'Language', 'Emotion', 'Relational'};
domain_color  = [0.75 0.15 0.20; 0.95 0.70 0.10; 0.95 0.55 0.15; 0.45 0.28 0.15; ...
                 0.15 0.45 0.75; 0.45 0.25 0.55; 0.40 0.70 0.25]; % roughly the paper's colours
task_domain = zeros(numel(task_names), 1);
for d = 1:numel(domain_prefix)
    task_domain(startsWith(task_names, [domain_prefix{d} '_'])) = d;
end
assert(all(task_domain > 0) && issorted(task_domain), 'Unexpected task order/domains');

%% Plot: heatmap per task and mode, bar chart averaged over tasks
col_titles = {'EDR binary', 'EDR weighted', 'EDR+LR'}; % paper's name for EDR continuous
diffs = {diff_bin, diff_edr, diff_edrlr};

figB = figure('Name', 'Our replication - tfMRI reconstruction accuracy difference', ...
    'Position', [100 100 1150 800]);
for k = 1:3
    % top: per-task heatmap
    ax = subplot(2,3,k);
    imagesc(diffs{k}, [-0.25 0.25]); colormap(ax, redblue);
    title(col_titles{k});
    set(ax, 'XTick', 5:5:th_mode, 'XTickLabel', [], 'YTick', []);
    if k == 1
        % domain colour bar and names down the left side
        for d = 1:numel(domain_prefix)
            rows = find(task_domain == d);
            rectangle('Position', [-0.6, rows(1)-0.5, 0.9, numel(rows)], ...
                'FaceColor', domain_color(d,:), 'EdgeColor', 'none', 'Clipping', 'off');
            text(-1.2, mean(rows), domain_label{d}, 'HorizontalAlignment', 'right', ...
                'FontSize', 8, 'Clipping', 'off');
        end
    end
    if k == 3
        cb = colorbar(ax);
        cb.Label.String = 'MSE difference (EDR - Geometry)';
    end

    % bottom: averaged over tasks
    ax2 = subplot(2,3,k+3);
    bar(mean(diffs{k},1), 'k', 'LineWidth', 2); title(col_titles{k});
    ylim([-0.25 0.25]); xlim([0 th_mode+0.5]);
    set(ax2, 'YTick', [-0.25 0 0.25]);
    xlabel('Modes');
    if k == 1
        ylabel({'mean mse. diff.', 'over tasks'});
    end
end
set(gcf, 'Color', [1 1 1]);
sgtitle('Our replication: differences in TA reconstruction for 47 HCP tasks (EDR - Geometry, N\leq20 modes)');

% the colorbar squeezes the top-right heatmap; line the heatmaps up with the bars
drawnow;
for k = 1:3
    top_ax = subplot(2,3,k); bot_ax = subplot(2,3,k+3);
    top_ax.Position([1 3]) = bot_ax.Position([1 3]);
end

exportgraphics(figB, fullfile('results', 'our_figure3B_replication.png'), 'Resolution', 200);
fprintf('Saved our_figure3B_replication.png\n');

% quick summary in the console
fprintf('\nMean MSE relative to Geometry, averaged across N<=20 modes and all 47 tasks:\n');
fprintf('  EDR binary:      %+.4f\n', mean(diff_bin(:)));
fprintf('  EDR continuous:  %+.4f\n', mean(diff_edr(:)));
fprintf('  EDR+LR:          %+.4f\n', mean(diff_edrlr(:)));
fprintf('(negative = better than Geometry)\n');
