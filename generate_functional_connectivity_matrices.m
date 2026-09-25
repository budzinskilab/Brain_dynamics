% generate_functional_connectivity_matrices.m
%
% PURPOSE: Functional counterpart of generate_connectivity_matrices.m:
%          empirical and reconstructed resting-state FC matrices plus the
%          functional long-range connections.
% ROLE:    Pipeline step 6. If
%          results/functional_connectivity_matrices.mat exists it only
%          redraws the figure (no data needed); otherwise needs steps 2,
%          3, 5 and the HCP resting scans (~30 min). Makes
%          results/functional_connectivity_matrices.png.
%
% Panels (group mean over 255 subjects, Glasser360 left hemisphere):
%   1. empirical FC
%   2. FC reconstructed from the EDR (continuous) modes
%   3. FC reconstructed from the EDR+LR modes
%   4. the functional long-range connections: pairs > 40 mm apart with
%      group FC > 0.5 (the paper's "LR connections", which Fig. 2 is
%      scored on)
%
% The reconstruction is the same as in step6b (per subject, vertex level,
% then parcellated and averaged). step6b only saved the MSE, so this
% redoes the pass and saves the group-mean FC at N = 20 and 200 in
% results/functional_connectivity_matrices.mat. If that file exists, only
% the figure is redrawn.

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');
addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath(genpath(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', 'colormaps_add'))); % redblue
addpath('functions');

out_file = fullfile('results', 'functional_connectivity_matrices.mat');
N_save = [20 200];
N_plot = 200; % as in the paper's Fig. 2

%% Group-mean empirical and reconstructed FC (slow, skipped if already saved)

if ~exist(out_file, 'file')
    basesfile = load(fullfile('results', 'step4_5_full_resolution_reconstruction.mat'), ...
        'eig_vec_geom', 'eig_vec_EDRbin', 'eig_vec_EDR', 'eig_vec_EDRLR');
    graph_names = {'Geometry', 'EDR binary', 'EDR continuous', 'EDR+LR'};
    graph_eigvecs = cellfun(@double, {basesfile.eig_vec_geom, basesfile.eig_vec_EDRbin, ...
        basesfile.eig_vec_EDR, basesfile.eig_vec_EDRLR}, 'UniformOutput', false);
    clear basesfile
    num_graphs = numel(graph_names);

    cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', 'fsLR_32k_cortex-lh_mask.txt'));
    cortex_ind = find(cortex);
    parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', 'fsLR_32k_Glasser360-lh.txt'));
    parc_cortex = parc(cortex_ind);
    num_parcels = numel(unique(parc_cortex(parc_cortex>0)));

    distfile = load(fullfile('results', 'step5_static_reconstruction.mat'), 'rr_parc');
    rr_parc = distfile.rr_parc;

    raw_dir = fullfile('neurodynamics_data', 'empirical', 'rfMRI_raw');
    subject_list = importdata(fullfile('neurodynamics_data', 'empirical', 'subject_list_HCP.txt'));
    num_sbj = numel(subject_list);

    FC_emp_sum = zeros(num_parcels);
    FC_recon_sum = zeros(num_parcels, num_parcels, num_graphs, numel(N_save));

    fprintf('Reconstructing resting FC for %d subjects (full vertex resolution)...\n', num_sbj);
    tic
    for s = 1:num_sbj
        cifti_path = fullfile(raw_dir, sprintf('%d_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii', subject_list(s)));
        d = ft_read_cifti(cifti_path);
        y_full = double(d.dtseries(d.brainstructure==1, :));
        y_full = y_full(cortex_ind, :); % [29696 x 1200]
        clear d

        FC_emp_sum = FC_emp_sum + corr(calc_parcellate(parc_cortex, y_full)');
        for g = 1:num_graphs
            basis = graph_eigvecs{g};
            beta_full = calc_eigendecomposition(y_full, basis, 'matrix'); % orthonormal bases, so fit once
            for ni = 1:numel(N_save)
                N = N_save(ni);
                y_hat_parc = calc_parcellate(parc_cortex, basis(:,1:N) * beta_full(1:N,:));
                FC_recon_sum(:,:,g,ni) = FC_recon_sum(:,:,g,ni) + corr(y_hat_parc');
            end
        end
        if mod(s,10) == 0
            fprintf('  ...subject %d/%d done (%.1f min elapsed)\n', s, num_sbj, toc/60);
        end
    end
    fprintf('Done in %.1f min.\n', toc/60);

    FC_emp_mean = FC_emp_sum / num_sbj;
    FC_recon_mean = FC_recon_sum / num_sbj; % [180 x 180 x graph x N_save]
    LR_mask = (rr_parc > 40) & (FC_emp_mean > 0.5); % as in step6b

    save(out_file, 'FC_emp_mean', 'FC_recon_mean', 'graph_names', 'N_save', 'LR_mask', 'num_sbj', '-v7.3');
    fprintf('Saved %s\n', out_file);
end

%% Plot

load(out_file, 'FC_emp_mean', 'FC_recon_mean', 'graph_names', 'N_save', 'LR_mask', 'num_sbj');
num_parcels = size(FC_emp_mean, 1);
ni = find(N_save == N_plot);
FC_EDR   = FC_recon_mean(:,:, strcmp(graph_names, 'EDR continuous'), ni);
FC_EDRLR = FC_recon_mean(:,:, strcmp(graph_names, 'EDR+LR'), ni);
FC_LR = FC_emp_mean; FC_LR(~LR_mask) = NaN;

mats = {FC_emp_mean, FC_EDR, FC_EDRLR, FC_LR};
titles = {sprintf('Empirical FC (group mean, %d subjects)', num_sbj), ...
    sprintf('FC from EDR modes (N = %d)', N_plot), ...
    sprintf('FC from EDR+LR modes (N = %d)', N_plot), ...
    sprintf('Functional LR connections (r > 0.5, > 40 mm; %d pairs)', nnz(triu(LR_mask,1)))};

% Same colours as the structural figure: grey for panels 1-2, blue/red
% for EDR+LR, red for LR. FC itself can't be split into EDR and LR parts,
% so the EDR+LR panel is coloured by where each pair's wiring comes from
% (the red pairs are the same ones as in the structural figure). One 0-1
% scale for all panels; the few slightly negative values (under 1% of
% pairs, min about -0.14) show as the lightest colour.
sc = load(fullfile('results', 'connectivity_matrices_parc.mat'), 'EDRLR_parc', 'exception_strength_parc');
offdiag = ~eye(num_parcels);
is_splice = (sc.exception_strength_parc ./ sc.EDRLR_parc > 0.5) & offdiag;

n_col = 256;
light = [0.94 0.95 0.97];
cmap_gray = flipud(bone(n_col));
cmap_edr = interp1([0 1], [light; 0.03 0.15 0.40], linspace(0, 1, n_col));
cmap_lr  = interp1([0 1], [light; 0.55 0.02 0.05], linspace(0, 1, n_col));
cmaps = {cmap_gray, cmap_gray, cmap_edr, cmap_lr};
clim_fc = [0 1];

fig = figure('Name', 'Functional connectivity matrices', 'Position', [100 100 1100 950], 'Color', [1 1 1]);
for k = 1:4
    ax = subplot(2,2,k);
    m = mats{k}; m(~offdiag) = NaN; % blank the diagonal (r = 1 by definition)
    if k == 3
        % place each value on the shared scale, colour by structural origin
        t = min(max((m - clim_fc(1)) / diff(clim_fc), 0), 1);
        idx = round(t * (n_col - 1)) + 1; idx(isnan(idx)) = 1;
        rgb = ones(num_parcels, num_parcels, 3);
        for c = 1:3
            layer = reshape(cmap_edr(idx, c), num_parcels, num_parcels);
            layer(is_splice) = cmap_lr(idx(is_splice), c);
            layer(isnan(m)) = 1;
            rgb(:, :, c) = layer;
        end
        image(rgb);
    else
        imagesc(m, 'AlphaData', ~isnan(m));
    end
    axis square; colormap(ax, cmaps{k}); caxis(ax, clim_fc);
    set(ax, 'Color', [1 1 1]);
    cb = colorbar(ax);
    title(titles{k}); xlabel('Parcel'); ylabel('Parcel');
    if k == 3
        % second colorbar (red) next to the blue one
        drawnow;
        cb_pos = cb.Position;
        cb_lr = axes('Position', [cb_pos(1) + 2.2*cb_pos(3), cb_pos(2), cb_pos(3), cb_pos(4)]);
        image(cb_lr, permute(cmap_lr, [1 3 2]));
        ticks = 0:0.2:1;
        set(cb_lr, 'YDir', 'normal', 'XTick', [], 'YTick', ticks * (n_col - 1) + 1, ...
            'YTickLabel', arrayfun(@(v) sprintf('%.1f', v), ticks, 'UniformOutput', false), ...
            'YAxisLocation', 'right', 'FontSize', cb.FontSize);
        cb_lr.YLabel.String = 'Correlation (r)';
        cb.YTickLabel = [];
        title(cb, 'EDR', 'FontSize', 8, 'FontWeight', 'bold');
        title(cb_lr, 'LR', 'FontSize', 8, 'FontWeight', 'bold');
        cb_edr = cb;
    else
        cb.Label.String = 'Correlation (r)';
    end
end
sgtitle('Functional Connectivity Matrices for EDR+LR');

% line the red colorbar up with the blue one after layout
drawnow;
cb_pos = cb_edr.Position;
cb_lr.Position = [cb_pos(1) + 2.2*cb_pos(3), cb_pos(2), cb_pos(3), cb_pos(4)];

exportgraphics(fig, fullfile('results', 'functional_connectivity_matrices.png'), 'Resolution', 200);
fprintf('Saved results/functional_connectivity_matrices.png\n');
