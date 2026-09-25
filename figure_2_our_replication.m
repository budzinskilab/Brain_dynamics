% figure_2_our_replication.m
%
% PURPOSE: Plots our Fig. 2 (long-range FC reconstruction error vs.
%          number of modes) in the paper's exact plot style.
% ROLE:    Figure. Runs from the included results (made by step 4). Makes
%          results/our_figure2_replication.png.
%
% Uses the same axes, order, labels and colours as the authors' plot
% (run_paper_figure2.m -> results/paper_figure2_authors_own_plot.png) so
% the two can sit side by side. step6b only ran 10 values of N, so our
% lines join those points (marked) rather than every N from 1 to 200.

load(fullfile('results', 'step6b_fc_longrange_reconstruction_fullres.mat'), ...
    'graph_names', 'mse_LR', 'N_test_list');
% mse_LR: [num_sbj x 4 graphs x numel(N_test_list)]

% paper's order: EDR binary, EDR continuous, Geometry, EDR+LR
paper_order = {'EDR binary', 'EDR continuous', 'Geometry', 'EDR+LR'};
[~, idx] = ismember(paper_order, graph_names);
mean_mse = squeeze(mean(mse_LR(:, idx, :), 1)); % [4 x numel(N_test_list)]

fig = figure('Name', 'Our replication - Fig 2C style - LR-FC reconstruction MSE vs modes');
plot(N_test_list, mean_mse', '.-', 'Linewidth', 2, 'MarkerSize', 12)
axis square; grid on; ylim([0, 0.05]); xlim([0, 200])
legend({'EDR binary', 'EDR weighted', 'Geometry', 'EDR+LR'}, 'Location', 'northeast')
ylabel('MSE'); xlabel('Modes');
title('Our replication: long-range FC reconstruction MSE vs. modes (full resolution, 255 subjects)');
set(gcf, 'Color', [1 1 1]);

exportgraphics(fig, fullfile('results', 'our_figure2_replication.png'), 'Resolution', 200);
fprintf('Saved results/our_figure2_replication.png\n');
