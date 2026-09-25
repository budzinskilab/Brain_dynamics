% step6b_fc_longrange_reconstruction_fullres.m
%
% PURPOSE: Main Fig. 2C replication: reconstructs each subject's resting-
%          state functional connectivity (FC) from the four bases at full
%          resolution, scores it on the long-range FC connections, and
%          runs the paper's statistics.
% ROLE:    Pipeline step 4 (after steps 2 and 3; needs the HCP resting
%          scans; ~1 h). Makes
%          results/step6b_fc_longrange_reconstruction_fullres.mat.
%
% Same analysis as step6 (see its header for the paper's numbers), but
% done the way the authors do it (PNAS_Figure_2_main.m): each subject's
% resting scan (neurodynamics_data/empirical/rfMRI_raw/) is fit and
% reconstructed at vertex level, then both the data and the
% reconstruction are parcellated to Glasser360, turned into FC, and
% scored on the long-range FC connections.
%

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');
addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

%% The four eigenmode bases (from step4_5)

fprintf('Loading eigenmode bases from step4_5 results...\n');
basesfile = load(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
    'eig_vec_geom', 'eig_vec_EDRbin', 'eig_vec_EDR', 'eig_vec_EDRLR');

graph_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
graph_eigvecs = {basesfile.eig_vec_geom, basesfile.eig_vec_EDRbin, ...
    basesfile.eig_vec_EDR, basesfile.eig_vec_EDRLR};
num_graphs = numel(graph_names);
num_modes = size(graph_eigvecs{1}, 2); % 200

%% Surface, cortex mask and parcellation (for scoring)

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
parc_name = 'Glasser360';

cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);

parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
num_parcels = numel(unique(parc_cortex(parc_cortex>0)));
fprintf('Loaded surface/parcellation: %d cortex vertices, %d parcels\n', numel(cortex_ind), num_parcels);

% step4_5 saves the bases for cortex vertices only, in cortex_ind order
graph_eigvecs_cortex = cellfun(@double, graph_eigvecs, 'UniformOutput', false);

%% Parcel distances (step5) and subject list

distfile = load(fullfile('results', 'step5_static_reconstruction.mat'), 'rr_parc');
rr_parc = distfile.rr_parc; % [180 x 180]

raw_dir = fullfile('neurodynamics_data', 'empirical', 'rfMRI_raw');
subject_list = importdata(fullfile('neurodynamics_data', 'empirical', 'subject_list_HCP.txt'));
num_sbj = numel(subject_list);
fprintf('Subject list: %d subjects\n', num_sbj);

N_test_list = [5 10 15 20 30 50 80 120 179 200]; % the paper uses 200
n_N = numel(N_test_list);

%% One pass over subjects: fit, reconstruct, parcellate, store FC
% Scoring waits until after the loop because the long-range mask needs the
% group-mean FC of all subjects. The stored 180x180 FC matrices come to
% about 1.9 GB.

FC_emp_all = zeros(num_parcels, num_parcels, num_sbj, 'single');
FC_recon_all = zeros(num_parcels, num_parcels, num_sbj, num_graphs, n_N, 'single');

fprintf('\nReading %d subjects'' full-vertex resting scans and reconstructing (this is the slow part)...\n', num_sbj);
tic
for s = 1:num_sbj
    cifti_path = fullfile(raw_dir, sprintf('%d_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii', subject_list(s)));
    d = ft_read_cifti(cifti_path);
    y_full = double(d.dtseries(d.brainstructure==1, :)); % [32492 x 1200], left hemisphere
    y_full = y_full(cortex_ind, :); % [29696 x 1200], cortex only
    clear d

    y_parc_emp = calc_parcellate(parc_cortex, y_full); % [180 x 1200]
    FC_emp_all(:,:,s) = single(corr(y_parc_emp'));

    for g = 1:num_graphs
        basis = graph_eigvecs_cortex{g}; % [29696 x 200]
        % fit all 200 modes once and truncate for smaller N; fine here
        % because step4_5 orthonormalized the bases (this doesn't work at
        % parcel level, which is why step6 refits for each N)
        beta_full = calc_eigendecomposition(y_full, basis, 'matrix'); % [200 x 1200]
        for ni = 1:n_N
            N = N_test_list(ni);
            y_hat_full = basis(:,1:N) * beta_full(1:N,:); % [29696 x 1200]
            y_hat_parc = calc_parcellate(parc_cortex, y_hat_full); % [180 x 1200]
            FC_recon_all(:,:,s,g,ni) = single(corr(y_hat_parc'));
        end
    end

    if mod(s,10) == 0
        fprintf('  ...subject %d/%d done (%.1f min elapsed)\n', s, num_sbj, toc/60);
    end
end
fprintf('Done in %.1f min.\n', toc/60);

%% Long-range FC mask from the group-mean empirical FC

th_dist = 40; th_corr = 0.5;
FC_mean = mean(double(FC_emp_all), 3);
LR_mask = (rr_parc > th_dist) & (FC_mean > th_corr);
LR_mask_triu = LR_mask & triu(true(num_parcels), 1);
n_LR_edges = nnz(LR_mask_triu);
fprintf('\nLR-FC mask: %d edges (>%gmm, empirical group FC >%g)\n', n_LR_edges, th_dist, th_corr);

%% MSE on the long-range connections, per subject, basis and N

mse_LR = nan(num_sbj, num_graphs, n_N);
for s = 1:num_sbj
    FC_emp_s = double(FC_emp_all(:,:,s));
    for g = 1:num_graphs
        for ni = 1:n_N
            FC_recon_s = double(FC_recon_all(:,:,s,g,ni));
            mse_LR(s,g,ni) = immse(FC_emp_s(LR_mask_triu), FC_recon_s(LR_mask_triu));
        end
    end
end

%% Paired t-tests across subjects at N = 20 and N = 200

idx_geom  = find(strcmp(graph_names, 'Geometry'));
idx_bin   = find(strcmp(graph_names, 'EDR binary'));
idx_edr   = find(strcmp(graph_names, 'EDR continuous'));
idx_edrlr = find(strcmp(graph_names, 'EDR+LR'));

fprintf('\n=== Mean LR-FC reconstruction MSE by N (lower = better), full vertex resolution ===\n');
fprintf('%-16s', 'N ->');
for ni = 1:n_N; fprintf('%10d', N_test_list(ni)); end
fprintf('\n');
for g = 1:num_graphs
    fprintf('%-16s', graph_names{g});
    for ni = 1:n_N
        fprintf('%10.5f', mean(mse_LR(:,g,ni)));
    end
    fprintf('\n');
end

for N_report = [20, 200]
    ni = find(N_test_list == N_report, 1);
    fprintf('\n=== Paired t-tests across %d subjects at N=%d (full vertex resolution) ===\n', num_sbj, N_report);
    [~, p1, ~, s1] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_geom,ni));
    fprintf('EDR+LR vs. geometry:         t(%d) = %+.3f, p = %.4g\n', s1.df, s1.tstat, p1);
    [~, p2, ~, s2] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_edr,ni));
    fprintf('EDR+LR vs. EDR continuous:   t(%d) = %+.3f, p = %.4g\n', s2.df, s2.tstat, p2);
    [~, p3, ~, s3] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_bin,ni));
    fprintf('EDR+LR vs. EDR binary:       t(%d) = %+.3f, p = %.4g\n', s3.df, s3.tstat, p3);
    [~, p4, ~, s4] = ttest(mse_LR(:,idx_edr,ni), mse_LR(:,idx_geom,ni));
    fprintf('EDR continuous vs. geometry: t(%d) = %+.3f, p = %.4g\n', s4.df, s4.tstat, p4);
end
fprintf('\n[paper''s Fig 2C, N=200: EDR+LR vs geometry P<0.0005, EDR+LR vs EDR continuous P<1e-4, ');
fprintf('EDR+LR vs EDR binary P<1e-4, EDR continuous vs geometry n.s.]\n');
fprintf('(negative t / lower MSE for the first-named graph means it reconstructs LR-FC better)\n');

%% Save

save(fullfile('results', 'step6b_fc_longrange_reconstruction_fullres.mat'), ...
    'graph_names', 'mse_LR', 'N_test_list', 'LR_mask_triu', 'n_LR_edges', 'num_sbj', '-v7.3');
fprintf('\nSaved results to results/step6b_fc_longrange_reconstruction_fullres.mat\n');
