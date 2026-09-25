% generate_our_figure3A.m
%
% PURPOSE: Our Fig. 3A: task reconstruction error vs. number of modes,
%          and each mode's contribution, for the 7 representative tasks.
% ROLE:    Figure (after step 2). Makes
%          results/our_figure3A_replication.png.
%
% The plotting follows the authors' Fig. 3A code, applied to our step4_5
% results so the two figures can be compared panel by panel.

%% Our results and the 7 example tasks

load(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
    'mse_all', 'graph_names', 'task_names');

idx_geom  = find(strcmp(graph_names, 'Geometry'));
idx_bin   = find(strcmp(graph_names, 'EDR binary'));
idx_edr   = find(strcmp(graph_names, 'EDR continuous'));
idx_edrlr = find(strcmp(graph_names, 'EDR+LR'));

% the authors' indices; task_names comes from the same data file, so they
% pick out the same tasks
representTask = [3, 11, 19, 30, 41, 44, 47];
fprintf('Representative tasks (matching the authors'' indices):\n');
for i = 1:numel(representTask)
    fprintf('  %d: %s\n', representTask(i), task_names{representTask(i)});
end

nMod = 200;

%% MSE for the 7 tasks, as [nMod x 7] like the authors' variables
mse_7tk_bin   = mse_all{idx_bin}(representTask, 1:nMod)';
mse_7tk_edr   = mse_all{idx_edr}(representTask, 1:nMod)';
mse_7tk_geom  = mse_all{idx_geom}(representTask, 1:nMod)';
mse_7tk_edrlr = mse_all{idx_edrlr}(representTask, 1:nMod)';

%% Figure

figA = figure('Name', 'Our replication - Fig 3A style - 7 representative tasks');
subplot(2,4,2)
plot(1:nMod, mse_7tk_bin./max(mse_7tk_bin), 'b-', 'linewidth', 2);grid on; title('EDR binary');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,3);axis square
plot(1:nMod, mse_7tk_edr./max(mse_7tk_edr), 'g-', 'linewidth', 2');grid on; title('EDR weighted');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,1);axis square
plot(1:nMod, mse_7tk_geom./max(mse_7tk_geom), 'm-', 'linewidth', 2);grid on;ylim([0,1]); title('Geometry');hold on;xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,4);
plot(1:nMod, mse_7tk_edrlr./max(mse_7tk_edrlr), 'k-', 'linewidth', 2);grid on; title('EDR+LR');xlim([1,nMod])
axis square;ylim([0, 1])
subplot(2,4,6);
plot((2:nMod), abs(diff([zeros(1,7); mse_7tk_bin(2:nMod,:)])), 'b-', 'linewidth', 2);grid on; title('EDR binary');
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,7);
plot((2:nMod), abs(diff([zeros(1,7); mse_7tk_edr(2:nMod,:)])), 'g-', 'linewidth', 2);grid on; title('EDR weighted');
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,5);
plot((2:nMod), abs(diff([zeros(1,7); mse_7tk_geom(2:nMod,:)])), 'm-', 'linewidth', 2);grid on; title('Geometry');
xlim([4,nMod]);ylim([0, 0.3])
subplot(2,4,8);
plot((2:nMod), abs(diff([zeros(1,7); mse_7tk_edrlr(2:nMod,:)])), 'k-', 'linewidth', 2);grid on; title('EDR+LR');
xlim([4,nMod]);ylim([0, 0.3])
set(gcf,'Color', [1 1 1]);
sgtitle('Our replication: normalized MSE (top) and per-mode contribution (bottom), same 7 representative tasks');

exportgraphics(figA, fullfile('results', 'our_figure3A_replication.png'), 'Resolution', 200);
fprintf('Saved our_figure3A_replication.png\n');
