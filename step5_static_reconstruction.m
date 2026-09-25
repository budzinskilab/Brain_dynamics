% step5_static_reconstruction.m
%
% PURPOSE: Quick parcel-level (180-region) version of the task
%          reconstruction; also saves the region-to-region distances used
%          by later steps.
% ROLE:    Pipeline step 3. Needs the Pang + Vohryzek data folders. Makes
%          results/step5_static_reconstruction.mat.
%
% A fast version of the Fig. 3 task reconstruction at Glasser360 (180
% parcels) instead of full vertex resolution. It differs from step4_5 in
% two ways:
%   1. everything is done at parcel level, and
%   2. it reconstructs the group-average map for each task, rather than
%      each subject separately. step4_5 follows the paper (per subject);
%      this version is kept to show how much the averaging order matters.
%
% It rebuilds the EDR and EDR+LR connectomes itself (same steps as
% Stage 1 of connectome_harmonics.m), so it doesn't need step 1.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

if ~exist('results', 'dir'); mkdir('results'); end

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
parc_name = 'Glasser360';

%% Surface, cortex mask and parcellation

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);

parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);
fprintf('Loaded %s parcellation: %d parcels\n', parc_name, num_parcels);

%% EDR and EDR+LR connectomes at parcel resolution

vertices_cortex = surface_midthickness.vertices(cortex_ind, :);
centroids = zeros(num_parcels, 3);
for p = 1:num_parcels
    centroids(p,:) = mean(vertices_cortex(parc_cortex==parcels(p), :), 1);
end
rr_parc = squareform(pdist(centroids));

load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
connectome_parc = calc_parcellate_matrix(parc_cortex, avgSC_L);
C_parc = connectome_parc / max(connectome_parc(:));
clear avgSC_L

NR = 60; NSTD = 3; DistRange = 40;
range_dist = max(rr_parc(:));
delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);
index_parc = floor(rr_parc/delta) + 1;
index_parc(index_parc > NR) = NR;

sc_density = cell(1, NR); sc_density_i = cell(1, NR); sc_density_j = cell(1, NR);
ycoor2 = nan(1, NR);
for n = 1:NR
    [idx_i, idx_j] = find(index_parc == n);
    idx = find(index_parc == n);
    sc_density{n} = C_parc(idx);
    sc_density_i{n} = idx_i;
    sc_density_j{n} = idx_j;
    if ~isempty(idx)
        ycoor2(n) = mean(C_parc(idx));
    end
end

fit_start = find(xcoor >= 10, 1);
fit_ind = fit_start:NR;
fit_ind = fit_ind(~isnan(ycoor2(fit_ind)));
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - ycoor2(fit_ind)).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
Afit = fminsearch(sse, [0.15, 0.18], options);
lambda = Afit(2);
fprintf('EDR fit (Glasser360): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

Clong = zeros(num_parcels, num_parcels);
for i = fit_start:NR
    if isempty(sc_density{i}); continue; end
    mv = mean(sc_density{i}); st = std(sc_density{i});
    ind_exc = find(sc_density{i} > mv + NSTD*st);
    for n = 1:numel(ind_exc)
        ii = sc_density_i{i}(ind_exc(n)); jj = sc_density_j{i}(ind_exc(n));
        if rr_parc(ii,jj) > DistRange
            Clong(ii,jj) = sc_density{i}(ind_exc(n));
        end
    end
end

EDR_conn = Afit(1)*exp(-Afit(2)*rr_parc);
EDR_LRE_parc = EDR_conn;
EDR_LRE_parc(Clong>0) = Clong(Clong>0);

fprintf('Long-range exceptions found: %d parcel pairs (out of %d)\n', ...
    nnz(triu(Clong,1)), num_parcels*(num_parcels-1)/2);

clear Clong index_parc sc_density sc_density_i sc_density_j C_parc connectome_parc

%% EDR continuous and EDR binary
%
% Built as in the authors' code, from Pang et al.'s fixed lambda =
% 0.12/mm rather than the fit above (the fit is only used as the EDR+LR
% baseline).
lambda_pang = 0.12;
EDR_conn_pang = single(exp(-lambda_pang * rr_parc));
EDR_conn_pang(1:(num_parcels+1):end) = 0;
EDR_conn_pang = EDR_conn_pang / max(EDR_conn_pang(:));

rng(1); % same seed as the authors' code
rand_prob = rand(size(EDR_conn_pang));
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
EDR_binary_parc = double(rand_prob < EDR_conn_pang);
EDR_binary_parc(1:(num_parcels+1):end) = 0;

%% Geometric eigenmodes (precomputed by Pang et al.)
%
% These come from the surface mesh rather than a graph. After
% parcellation they're no longer exactly orthonormal, which is fine for
% the least-squares fit below.

geom_evec_vertex = dlmread(fullfile(eigenmode_toolbox_dir, 'data', 'results', ...
    sprintf('basis_geometric_%s-%s_evec_200.txt', mesh_interest, hemisphere)));
eig_vec_geom_full = calc_parcellate(parc, geom_evec_vertex);
clear geom_evec_vertex

%% Eigenmodes for all four graphs
%
% Dense eig() rather than eigs(): at 180 parcels it's cheap, and EDR binary
% has 21 isolated nodes, i.e. many repeated zero eigenvalues, which eigs()
% handled badly (it returned complex eigenvectors).
num_modes = 100; % well past the 1-20 modes the paper focuses on

fprintf('\nComputing eigenmodes (Geometry, EDR binary, EDR continuous, EDR+LR) at %d-parcel resolution...\n', num_parcels);
eig_vec_geom = eig_vec_geom_full(:, 1:num_modes);
[eig_vec_EDRbin, eig_val_EDRbin] = calc_network_eigenmode_dense(EDR_binary_parc, num_modes);
[eig_vec_EDR, eig_val_EDR] = calc_network_eigenmode_dense(EDR_conn_pang, num_modes);
[eig_vec_EDRLR, eig_val_EDRLR] = calc_network_eigenmode_dense(EDR_LRE_parc, num_modes);

graph_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
graph_colors = {'m', 'b', 'g', 'k'};
graph_eigvecs = {eig_vec_geom, eig_vec_EDRbin, eig_vec_EDR, eig_vec_EDRLR};
num_graphs = numel(graph_names);

%% Task fMRI data: group-average map per task, parcellated

fprintf('\nLoading task fMRI data (47 contrasts, 255-subject group average)...\n');
data = load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
task_names = fieldnames(data.zstat);
num_tasks = numel(task_names);

task_maps_parc = zeros(num_parcels, num_tasks);
for t = 1:num_tasks
    group_avg = nanmean(data.zstat.(task_names{t}), 2); % [32492 x 1]
    task_maps_parc(:, t) = calc_parcellate(parc, group_avg);
end
clear data

%% Reconstruction accuracy vs. number of modes

fprintf('\nRunning reconstruction (accuracy vs. N modes, %d tasks x %d graphs)...\n', num_tasks, num_graphs);

corr_all = cell(1, num_graphs);
mse_all  = cell(1, num_graphs);
for g = 1:num_graphs
    corr_all{g} = nan(num_tasks, num_modes);
    mse_all{g}  = nan(num_tasks, num_modes);
end

for t = 1:num_tasks
    y = task_maps_parc(:, t);

    for g = 1:num_graphs
        basis = graph_eigvecs{g};
        beta = calc_eigendecomposition(y, basis, 'matrix');
        for N = 1:num_modes
            y_hat = double(basis(:,1:N) * beta(1:N));
            corr_all{g}(t,N) = corr(y, y_hat);
            mse_all{g}(t,N)  = immse(y, y_hat);
        end
    end
end
fprintf('Done.\n');

mean_corr = cellfun(@(x) mean(x,1), corr_all, 'UniformOutput', false);
mean_mse  = cellfun(@(x) mean(x,1), mse_all, 'UniformOutput', false);

%% Paired comparisons across tasks at N = 20

N_test = min(20, num_modes);
idx_EDR   = find(strcmp(graph_names, 'EDR continuous'));
idx_EDRLR = find(strcmp(graph_names, 'EDR+LR'));
idx_geom  = find(strcmp(graph_names, 'Geometry'));

fprintf('\n=== Paired comparisons at N=%d modes, across %d tasks ===\n', N_test, num_tasks);
for g = 1:num_graphs
    fprintf('%-16s mean corr = %.4f, mean MSE = %.4f\n', graph_names{g}, ...
        mean_corr{g}(N_test), mean_mse{g}(N_test));
end

[~, p_corr_edrlr_edr, ~, s_corr_edrlr_edr] = ttest(corr_all{idx_EDRLR}(:,N_test), corr_all{idx_EDR}(:,N_test));
[~, p_mse_edrlr_edr,  ~, s_mse_edrlr_edr]  = ttest(mse_all{idx_EDRLR}(:,N_test),  mse_all{idx_EDR}(:,N_test));
fprintf('\nPrimary comparison, EDR+LR vs. EDR-only:\n');
fprintf('  Correlation: paired t(%d) = %.3f, p = %.4g\n', s_corr_edrlr_edr.df, s_corr_edrlr_edr.tstat, p_corr_edrlr_edr);
fprintf('  MSE:         paired t(%d) = %.3f, p = %.4g\n', s_mse_edrlr_edr.df, s_mse_edrlr_edr.tstat, p_mse_edrlr_edr);

[~, p_corr_edrlr_geom, ~, s_corr_edrlr_geom] = ttest(corr_all{idx_EDRLR}(:,N_test), corr_all{idx_geom}(:,N_test));
[~, p_mse_edrlr_geom,  ~, s_mse_edrlr_geom]  = ttest(mse_all{idx_EDRLR}(:,N_test),  mse_all{idx_geom}(:,N_test));
fprintf('\nReplication check, EDR+LR vs. Geometry (Vohryzek''s own headline comparison):\n');
fprintf('  Correlation: paired t(%d) = %.3f, p = %.4g\n', s_corr_edrlr_geom.df, s_corr_edrlr_geom.tstat, p_corr_edrlr_geom);
fprintf('  MSE:         paired t(%d) = %.3f, p = %.4g\n', s_mse_edrlr_geom.df, s_mse_edrlr_geom.tstat, p_mse_edrlr_geom);

%% Save

save(fullfile('results', 'step5_static_reconstruction.mat'), ...
    'task_names', 'graph_names', 'corr_all', 'mse_all', ...
    'eig_vec_geom', 'eig_vec_EDRbin', 'eig_val_EDRbin', ...
    'eig_vec_EDR', 'eig_val_EDR', 'eig_vec_EDRLR', 'eig_val_EDRLR', ...
    'EDR_conn', 'EDR_conn_pang', 'EDR_binary_parc', 'EDR_LRE_parc', 'rr_parc', 'Afit', 'lambda', '-v7.3');
fprintf('\nSaved results to results/step5_static_reconstruction.mat\n');

%% Figure

figure('Name', 'Stage 5 - static reconstruction, all four graphs');
subplot(1,2,1);
hold on
for g = 1:num_graphs
    plot(1:num_modes, mean_corr{g}, [graph_colors{g} '-'], 'LineWidth', 2);
end
xlabel('Number of modes (N)'); ylabel('Correlation with empirical activity');
legend(graph_names, 'Location', 'southeast'); grid on
title('Reconstruction accuracy (correlation)'); xlim([1 num_modes]);

subplot(1,2,2);
hold on
for g = 1:num_graphs
    plot(1:num_modes, mean_mse{g}, [graph_colors{g} '-'], 'LineWidth', 2);
end
xlabel('Number of modes (N)'); ylabel('MSE vs. empirical activity');
legend(graph_names, 'Location', 'northeast'); grid on
title('Reconstruction error (MSE)'); xlim([1 num_modes]);

sgtitle(sprintf('Glasser360, %d-task average, group-level HCP activation maps', num_tasks));
saveas(gcf, fullfile('results', 'step5_static_reconstruction.png'));
fprintf('Saved figure to results/step5_static_reconstruction.png\n');
