% step6_fc_longrange_reconstruction.m
%
% PURPOSE: Earlier parcel-resolution approximation of the long-range FC
%          reconstruction (the paper's Fig. 2C).
% ROLE:    Reference only; replaced by
%          step6b_fc_longrange_reconstruction_fullres.m.
%
% Fig. 2C: reconstruct resting-state FC from each basis and compare the
% error on the long-range FC connections, paired across the 255 subjects.
% The paper reports (paired t-tests, N = 200 modes):
%   EDR+LR vs. geometry         P < 0.0005
%   EDR+LR vs. EDR continuous   P < 1e-4
%   EDR+LR vs. EDR binary       P < 1e-4
%   EDR continuous vs. geometry n.s.
%
% This version uses Pang's resting time series, which are already
% parcellated (S255_resting_empirical_Glasser360.mat), so the fit is done
% on 180 parcels instead of at vertex level like the paper. That caps N
% at 179, and the fit is redone for each N because the parcellated bases
% aren't orthonormal. step6b does the proper vertex-level version.
%
% In Pang's bilateral file the left hemisphere is labels 181-360 (checked
% against the lh/rh label files). If every basis gives near-zero
% correlations, check this mapping first.

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
num_modes = size(graph_eigvecs{1}, 2);
fprintf('  %d graphs x %d modes\n', num_graphs, num_modes);

%% Surface, cortex mask and parcellation (to parcellate the bases)

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
fprintf('Parcellating eigenmode bases to %s (%d parcels)...\n', parc_name, num_parcels);

graph_eigvecs_parc = cellfun(@(v) calc_parcellate(parc_cortex, double(v)), graph_eigvecs, 'UniformOutput', false);

%% Parcel distances (step5) and resting time series (Pang)

fprintf('Loading parcel distances and resting-state time series...\n');
distfile = load(fullfile('results', 'step5_static_reconstruction.mat'), 'rr_parc');
rr_parc = distfile.rr_parc; % [180 x 180] centroid distances, mm

restfile = load(fullfile(eigenmode_toolbox_dir, 'data', 'results', 'S255_resting_empirical_Glasser360.mat'));
ts_bilateral = restfile.resting_emp.time_series; % [360 x 1200 x 255]
ts_lh = ts_bilateral(181:360, :, :);              % left hemisphere = labels 181-360
clear ts_bilateral restfile

num_sbj = size(ts_lh, 3);
T = size(ts_lh, 2);
fprintf('  %d subjects, %d timepoints, %d LH parcels\n', num_sbj, T, size(ts_lh,1));

%% Empirical FC per subject; the group mean defines the long-range mask

fprintf('Computing empirical FC per subject...\n');
FC_emp_all = zeros(num_parcels, num_parcels, num_sbj, 'single');
for s = 1:num_sbj
    FC_emp_all(:,:,s) = single(corr(ts_lh(:,:,s)'));
end
FC_mean = mean(FC_emp_all, 3);

th_dist = 40; th_corr = 0.5;
triu_ind = find(triu(ones(num_parcels), 1));
LR_mask = (rr_parc > th_dist) & (FC_mean > th_corr);
LR_mask_triu = LR_mask & triu(true(num_parcels), 1);
n_LR_edges = nnz(LR_mask_triu);
fprintf('LR-FC mask: %d edges (>%gmm, empirical group FC >%g)\n', n_LR_edges, th_dist, th_corr);

if n_LR_edges < 10
    error(['Only %d LR-FC edges found -- too few to draw a conclusion. ' ...
        'This usually means the hemisphere/parcel-order mapping for the ' ...
        'bilateral resting-state data is wrong (see header note).'], n_LR_edges);
end

%% Reconstruct each subject's FC and score it on the long-range connections
%
% N stops at 179: 200 modes can't be fit in a 180-parcel space (the paper
% can use 200 because it fits at vertex level).

N_max = num_parcels - 1; % 179: highest N at parcel resolution
N_test_list = [5 10 15 20 30 50 80 120 N_max];
fprintf('\nReconstructing FC for %d subjects x %d graphs x N in {%s} (capped at %d, see header note)...\n', ...
    num_sbj, num_graphs, num2str(N_test_list), N_max);

mse_LR = nan(num_sbj, num_graphs, numel(N_test_list));

tic
for s = 1:num_sbj
    y = double(ts_lh(:,:,s)); % [180 x 1200]
    FC_emp_s = double(FC_emp_all(:,:,s));

    for g = 1:num_graphs
        basis_parc = graph_eigvecs_parc{g}; % [180 x 200]
        for ni = 1:numel(N_test_list)
            N = N_test_list(ni);
            basis_N = basis_parc(:,1:N);
            beta = calc_eigendecomposition(y, basis_N, 'matrix'); % [N x 1200], refit for each N
            y_hat = basis_N * beta;
            FC_recon = corr(y_hat');
            mse_LR(s,g,ni) = immse(FC_emp_s(LR_mask_triu), FC_recon(LR_mask_triu));
        end
    end
    if mod(s,50) == 0
        fprintf('  ...subject %d/%d done (%.1f s elapsed)\n', s, num_sbj, toc);
    end
end
fprintf('Done in %.1f s.\n', toc);

%% Paired t-tests across subjects at N = 20 and N = 179 (closest to the paper's 200)

idx_geom  = find(strcmp(graph_names, 'Geometry'));
idx_bin   = find(strcmp(graph_names, 'EDR binary'));
idx_edr   = find(strcmp(graph_names, 'EDR continuous'));
idx_edrlr = find(strcmp(graph_names, 'EDR+LR'));

fprintf('\n=== Mean LR-FC reconstruction MSE by N (lower = better) ===\n');
fprintf('%-16s', 'N ->');
for ni = 1:numel(N_test_list); fprintf('%10d', N_test_list(ni)); end
fprintf('\n');
for g = 1:num_graphs
    fprintf('%-16s', graph_names{g});
    for ni = 1:numel(N_test_list)
        fprintf('%10.5f', mean(mse_LR(:,g,ni)));
    end
    fprintf('\n');
end

for N_report = [20, N_max]
    ni = find(N_test_list == N_report, 1);
    fprintf('\n=== Paired t-tests across %d subjects at N=%d ===\n', num_sbj, N_report);
    [~, p1, ~, s1] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_geom,ni));
    fprintf('EDR+LR vs. geometry:         t(%d) = %+.3f, p = %.4g\n', s1.df, s1.tstat, p1);
    [~, p2, ~, s2] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_edr,ni));
    fprintf('EDR+LR vs. EDR continuous:   t(%d) = %+.3f, p = %.4g\n', s2.df, s2.tstat, p2);
    [~, p3, ~, s3] = ttest(mse_LR(:,idx_edrlr,ni), mse_LR(:,idx_bin,ni));
    fprintf('EDR+LR vs. EDR binary:       t(%d) = %+.3f, p = %.4g\n', s3.df, s3.tstat, p3);
    [~, p4, ~, s4] = ttest(mse_LR(:,idx_edr,ni), mse_LR(:,idx_geom,ni));
    fprintf('EDR continuous vs. geometry: t(%d) = %+.3f, p = %.4g\n', s4.df, s4.tstat, p4);
end
fprintf('\n[paper''s Fig 2C, at their N=200 in full 29,696-vertex space: EDR+LR vs geometry P<0.0005, ');
fprintf('EDR+LR vs EDR continuous P<1e-4, EDR+LR vs EDR binary P<1e-4, EDR continuous vs geometry n.s.]\n');
fprintf('(negative t / lower MSE for the first-named graph means it reconstructs LR-FC better)\n');

%% Save
save(fullfile('results', 'step6_fc_longrange_reconstruction.mat'), ...
    'graph_names', 'mse_LR', 'N_test_list', 'LR_mask_triu', 'n_LR_edges', 'num_sbj', '-v7.3');
fprintf('\nSaved results to results/step6_fc_longrange_reconstruction.mat\n');
