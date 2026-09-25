% step4_5_full_resolution_reconstruction.m
%
% PURPOSE: Reconstructs all 47 HCP task maps for each of the 255 subjects
%          from the four eigenmode bases (Geometry, EDR binary, EDR
%          continuous, EDR+LR) at full vertex resolution.
% ROLE:    Pipeline step 2 (after step 1). Makes
%          results/step4_5_full_resolution_reconstruction.mat (the four
%          bases + task reconstruction results).
%
% Replicates Fig. 3A/B of Vohryzek et al. (2025): reconstruction accuracy
% (correlation and MSE) vs. number of modes. step5 runs the same analysis
% at parcel resolution as a quick check.
%
% Notes on the method (following the paper):
%   - EDR+LR comes from step 1. EDR continuous and EDR binary use Pang et
%     al.'s fixed lambda = 0.12/mm instead of our fitted value, otherwise
%     EDR continuous and EDR+LR would share the same baseline.
%   - Eigenmodes are computed per vertex, but the error is measured after
%     parcellating both the data and the reconstruction to Glasser360.
%   - Each subject is reconstructed separately and the error is averaged
%     over subjects afterwards. Averaging the maps first smooths away the
%     subject-level detail and makes all four bases look alike.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

if ~exist('results', 'dir'); mkdir('results'); end

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
num_modes = 200; % same as the paper

%% Surface and cortex mask

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);

fprintf('Loaded surface: %d total vertices, %d cortex vertices\n', num_vertices, numel(cortex_ind));

% Glasser360 is only used to score the reconstructions, not to build the modes
parc_name = 'Glasser360';
parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
fprintf('Loaded %s parcellation for evaluation: %d parcels\n', ...
    parc_name, numel(unique(parc_cortex(parc_cortex>0))));

%% EDR+LR connectome from step 1

checkpoint_file = fullfile('results', 'EDR_LR_connectome_full_resolution.mat');
if ~exist(checkpoint_file, 'file')
    error(['Checkpoint file not found: %s\nRun Stage 2 of connectome_harmonics.m ' ...
        'first (the connectome-derivation part only -- the eigendecomposition ' ...
        'that used to crash there is superseded by this script).'], checkpoint_file);
end
fprintf('Loading checkpointed EDR+LR connectome and EDR fit...\n');
load(checkpoint_file, 'EDR_LRE', 'Afit', 'lambda');
fprintf('Loaded EDR fit: A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% EDR continuous and EDR binary (Pang et al.'s lambda = 0.12/mm)

lambda_pang = 0.12;

tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed %dx%d vertex distance matrix in %.1f s\n', size(rr,1), size(rr,2), toc);

Pspace = single(exp(-lambda_pang * double(rr)));
N = size(Pspace, 1);
Pspace(1:(N+1):end) = 0;
Pspace = Pspace / max(Pspace(:));

EDR_conn_full = Pspace;

% EDR binary: keep each vertex pair with probability Pspace. Same seed as
% the authors' code, and drawn in double precision like theirs (single
% gives a different random sequence for the same seed).
rng(1);
rand_prob = rand(size(Pspace));
rand_prob = triu(rand_prob,1) + triu(rand_prob,1)';
EDR_binary_full = single(rand_prob < Pspace);
clear rand_prob Pspace rr

EDR_conn_full(1:(N+1):end) = 0;
EDR_binary_full(1:(N+1):end) = 0;

%% Eigenmodes (Geometry modes are precomputed by Pang et al.)

fprintf('\nEigendecomposing EDR continuous (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDR, eig_val_EDR] = calc_network_eigenmode_lowmem(EDR_conn_full, num_modes); toc
clear EDR_conn_full

fprintf('\nEigendecomposing EDR binary (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDRbin, eig_val_EDRbin] = calc_network_eigenmode_lowmem(EDR_binary_full, num_modes); toc
clear EDR_binary_full

fprintf('\nEigendecomposing EDR+LR (%dx%d, %d modes)...\n', N, N, num_modes);
tic; [eig_vec_EDRLR, eig_val_EDRLR] = calc_network_eigenmode_lowmem(EDR_LRE, num_modes); toc
clear EDR_LRE

fprintf('\nLoading precomputed Geometry eigenmodes...\n');
geom_evec_full = dlmread(fullfile(eigenmode_toolbox_dir, 'data', 'results', ...
    sprintf('basis_geometric_%s-%s_evec_200.txt', mesh_interest, hemisphere)));
eig_vec_geom = single(geom_evec_full(cortex_ind, 1:num_modes));
clear geom_evec_full

graph_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
graph_colors = {'m', 'b', 'g', 'k'};
num_graphs = numel(graph_names);

%% Orthonormalize each basis (economy QR)
%
% The loop below fits all 200 modes once and truncates to N modes, which
% is only valid for an orthonormal basis. The EDR bases already are, but
% the Geometry modes are orthonormal under a surface-area-weighted inner
% product, not the plain dot product. QR without pivoting keeps the span
% of the first N columns for every N, so the reconstructions are the same;
% it just makes the truncation valid. step6b reuses these saved bases.
graph_eigvecs_raw = {eig_vec_geom, eig_vec_EDRbin, eig_vec_EDR, eig_vec_EDRLR};
graph_eigvecs = cell(1, num_graphs);
for g = 1:num_graphs
    [Qg, ~] = qr(double(graph_eigvecs_raw{g}), 0);
    graph_eigvecs{g} = single(Qg);
end
eig_vec_geom  = graph_eigvecs{1};
eig_vec_EDRbin = graph_eigvecs{2};
eig_vec_EDR    = graph_eigvecs{3};
eig_vec_EDRLR  = graph_eigvecs{4};
clear graph_eigvecs_raw Qg

%% Task fMRI data (all 255 subjects, not group-averaged)

fprintf('\nLoading task fMRI data (47 contrasts, 255 subjects each)...\n');
data = load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'));
task_names = fieldnames(data.zstat);
num_tasks = numel(task_names);
num_sbj = size(data.zstat.(task_names{1}), 2);
fprintf('  %d tasks x %d subjects\n', num_tasks, num_sbj);

%% Reconstruction accuracy vs. number of modes
%
% The coefficients are fit at vertex level for all subjects at once, then
% scored per subject at parcel level and averaged over subjects. Since
% parcellating is just averaging, we parcellate each basis once instead of
% every reconstruction.

fprintf('\nParcellating eigenmode bases to %s for reconstruction-error evaluation...\n', parc_name);
graph_eigvecs_parc = cellfun(@(v) calc_parcellate(parc_cortex, double(v)), graph_eigvecs, 'UniformOutput', false);

fprintf('\nRunning reconstruction (accuracy vs. N modes, %d tasks x %d subjects x %d graphs, %d modes)...\n', ...
    num_tasks, num_sbj, num_graphs, num_modes);

corr_all = cell(1, num_graphs);
mse_all  = cell(1, num_graphs);
for g = 1:num_graphs
    corr_all{g} = nan(num_tasks, num_modes);
    mse_all{g}  = nan(num_tasks, num_modes);
end

tic
for t = 1:num_tasks
    y = double(data.zstat.(task_names{t})(cortex_ind, :)); % [Ncortex x num_sbj]
    y_parc = calc_parcellate(parc_cortex, y);               % [num_parcels x num_sbj]

    for g = 1:num_graphs
        basis = double(graph_eigvecs{g});
        basis_parc = graph_eigvecs_parc{g};
        beta = calc_eigendecomposition(y, basis, 'matrix'); % [num_modes x num_sbj]
        for N_modes = 1:num_modes
            y_hat_parc = basis_parc(:,1:N_modes) * beta(1:N_modes,:); % [num_parcels x num_sbj]

            % correlation and MSE per subject, then averaged
            An = y_parc - mean(y_parc,1);
            Bn = y_hat_parc - mean(y_hat_parc,1);
            r_sbj   = sum(An.*Bn,1) ./ sqrt(sum(An.^2,1) .* sum(Bn.^2,1));
            mse_sbj = mean((y_parc - y_hat_parc).^2, 1);

            corr_all{g}(t,N_modes) = mean(r_sbj, 'omitnan');
            mse_all{g}(t,N_modes)  = mean(mse_sbj, 'omitnan');
        end
    end
    fprintf('  ...task %d/%d (%s) done (%.1f s elapsed)\n', t, num_tasks, task_names{t}, toc);
end
fprintf('Done in %.1f s.\n', toc);
clear data

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

save(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
    'task_names', 'graph_names', 'corr_all', 'mse_all', ...
    'eig_vec_geom', 'eig_vec_EDRbin', 'eig_val_EDRbin', ...
    'eig_vec_EDR', 'eig_val_EDR', 'eig_vec_EDRLR', 'eig_val_EDRLR', ...
    'Afit', 'lambda', '-v7.3');
fprintf('\nSaved results to results/step4_5_full_resolution_reconstruction.mat\n');

%% Figure

figure('Name', 'Full-resolution static reconstruction, all four graphs');
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

sgtitle(sprintf('Full vertex resolution (%d cortex vertices), %d-task average', numel(cortex_ind), num_tasks));
saveas(gcf, fullfile('results', 'step4_5_full_resolution_reconstruction.png'));
fprintf('Saved figure to results/step4_5_full_resolution_reconstruction.png\n');
