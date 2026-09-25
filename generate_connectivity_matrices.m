% generate_connectivity_matrices.m
%
% PURPOSE: Plots the four structural connectivity matrices (Connectome,
%          EDR baseline, EDR+LR, structural long-range exceptions) at
%          Glasser360 resolution.
% ROLE:    Pipeline step 5 (after step 1). Makes
%          results/connectivity_matrices.png and
%          results/connectivity_matrices_parc.mat.
%
% Everything is built at vertex level exactly as in Stage 2 of
% connectome_harmonics.m (reusing its fitted A and lambda), then averaged
% down to Glasser360 so the matrices are readable. A 29,696 x 29,696
% heatmap just looks like noise.

%% Surface, cortex mask and parcellation

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');
addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

hemisphere = 'lh'; surface_interest = 'fsLR_32k'; mesh_interest = 'midthickness'; parc_name = 'Glasser360';

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
vertices_cortex = surface_midthickness.vertices(cortex_ind, :);

parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
num_parcels = numel(unique(parc_cortex(parc_cortex>0)));
fprintf('Loaded surface/parcellation: %d cortex vertices, %d parcels\n', numel(cortex_ind), num_parcels);

%% EDR fit from step 1

load(fullfile('results', 'EDR_LR_connectome_full_resolution.mat'), 'EDR_LRE', 'Afit', 'lambda');
fprintf('Loaded EDR fit: A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% Empirical connectome and vertex distances

fprintf('\nLoading empirical connectome (Pang''s intact copy)...\n');
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
C = avgSC_L / max(avgSC_L(:));
clear avgSC_L

tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed %dx%d vertex distance matrix in %.1f s\n', size(rr,1), size(rr,2), toc);

%% EDR baseline and long-range exception mask (same settings as Stage 2)

NR = 400; NRini = 20; NRfin = 380; NSTD = 3; DistRange = 40;

EDR_conn_full = single(Afit(1) * exp(-Afit(2) * double(rr)));

range_dist = max(rr(:));
delta = range_dist / NR;
index = uint16(min(floor(double(rr)/delta) + 1, NR));
bin_mean = accumarray(index(:), C(:), [NR,1], @mean, single(NaN));
bin_std  = accumarray(index(:), C(:), [NR,1], @std, single(NaN));

local_mean = bin_mean(index);
local_std = bin_std(index);
is_exception = (index >= NRini) & (index <= NRfin) & (rr > DistRange) & (C > local_mean + NSTD*local_std);
clear local_mean local_std index bin_mean bin_std rr

n_exceptions = nnz(triu(is_exception,1));
fprintf('Long-range exceptions: %d vertex pairs\n', n_exceptions);

%% Parcellate to Glasser360
%
% Each parcel pair covers roughly 165 x 165 vertex pairs, so almost every
% pair contains at least one exception. Instead of a yes/no mask we save
% the fraction of exception vertex pairs in each block (density).
exception_density_parc = calc_parcellate_matrix(parc_cortex, single(is_exception));

% For the plot we use the strength the splice adds (the empirical value
% at each exception, zero elsewhere), so all four panels have the same units.
C_exc = C; C_exc(~is_exception) = 0;
exception_strength_parc = calc_parcellate_matrix(parc_cortex, C_exc); clear C_exc

fprintf('\nParcellating to %s (%d parcels)...\n', parc_name, num_parcels);
C_parc = calc_parcellate_matrix(parc_cortex, C); clear C
EDR_parc = calc_parcellate_matrix(parc_cortex, EDR_conn_full); clear EDR_conn_full
EDRLR_parc = calc_parcellate_matrix(parc_cortex, EDR_LRE); clear EDR_LRE
clear is_exception

save(fullfile('results', 'connectivity_matrices_parc.mat'), ...
    'C_parc', 'EDR_parc', 'EDRLR_parc', 'exception_strength_parc', 'exception_density_parc', ...
    'Afit', 'lambda', '-v7.3');

%% Plot: Connectome, EDR, EDR+LR, long-range exceptions
%
% One shared log colour scale for all panels. On a linear scale EDR+LR
% looks the same as EDR, because the EDR values near the diagonal (~1e-2)
% drown out the splice (~1e-5). The diagonal is blanked and zeros show as
% white.
%
% EDR+LR is drawn in two colours (as in figure_EDRLR_components.m): blue
% where a pair's strength comes mostly from the EDR baseline, red where it
% comes mostly from the splice.

mats = {C_parc, EDR_parc, EDRLR_parc, exception_strength_parc};
titles = {'Connectome (empirical, group-averaged)', ...
    sprintf('EDR baseline (\\lambda = %.3f /mm)', lambda), ...
    'EDR+LR (baseline + long-range splice)', ...
    sprintf('Structural long-range exceptions (SC; %d vertex pairs)', n_exceptions)};

offdiag = ~eye(num_parcels);
is_splice = (exception_strength_parc ./ EDRLR_parc > 0.5) & offdiag;
for k = 1:numel(mats)
    mats{k}(~offdiag) = NaN;
    mats{k}(mats{k} <= 0) = NaN;
end
all_vals = cell2mat(cellfun(@(m) m(~isnan(m)), mats, 'UniformOutput', false)');
clim_shared = [prctile(all_vals, 1), max(all_vals)];

fig = figure('Name', 'Connectivity matrices (Glasser360, from full-resolution derivation)', ...
    'Position', [100 100 1100 950]);
cmap = flipud(bone);
% blue for EDR, red for LR (same as figure_EDRLR_components.m)
n_col = 256;
light = [0.94 0.95 0.97];
cmap_edr = interp1([0 1], [light; 0.03 0.15 0.40], linspace(0, 1, n_col));
cmap_lr  = interp1([0 1], [light; 0.55 0.02 0.05], linspace(0, 1, n_col));
cmaps = {cmap, cmap, cmap_edr, cmap_lr};

for k = 1:numel(mats)
    ax = subplot(2,2,k);
    if k == 3
        % EDR+LR: place each value on the shared scale, then colour it
        % blue or red depending on where its strength comes from
        t = (log10(mats{k}) - log10(clim_shared(1))) / diff(log10(clim_shared));
        idx = round(min(max(t, 0), 1) * (n_col - 1)) + 1;
        idx(isnan(idx)) = 1;
        rgb = ones(num_parcels, num_parcels, 3);
        for c = 1:3
            layer = reshape(cmap_edr(idx, c), num_parcels, num_parcels);
            layer(is_splice) = cmap_lr(idx(is_splice), c);
            layer(isnan(mats{k})) = 1;
            rgb(:, :, c) = layer;
        end
        image(rgb);
    else
        imagesc(mats{k}, 'AlphaData', ~isnan(mats{k}));
    end
    axis square; colormap(ax, cmaps{k});
    set(ax, 'ColorScale', 'log', 'CLim', clim_shared, 'Color', [1 1 1]);
    cb = colorbar(ax);
    title(titles{k}); xlabel('Parcel'); ylabel('Parcel');
    if k == 3
        % second colorbar (red) next to the blue one
        drawnow;
        cb_pos = cb.Position;
        cb_lr = axes('Position', [cb_pos(1) + 2.2*cb_pos(3), cb_pos(2), cb_pos(3), cb_pos(4)]);
        image(cb_lr, permute(cmap_lr, [1 3 2]));
        tick_exp = 2*ceil(log10(clim_shared(1))/2):2:floor(log10(clim_shared(2))); % same ticks as the other colorbars
        tick_pos = (tick_exp - log10(clim_shared(1))) / diff(log10(clim_shared)) * (n_col - 1) + 1;
        set(cb_lr, 'YDir', 'normal', 'XTick', [], 'YTick', tick_pos, 'YAxisLocation', 'right', ...
            'YTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), tick_exp, 'UniformOutput', false), ...
            'TickLabelInterpreter', 'tex', 'FontSize', cb.FontSize);
        cb_lr.YLabel.String = 'Connection strength (normalized, log scale)';
        cb.YTickLabel = [];
        title(cb, 'EDR', 'FontSize', 8, 'FontWeight', 'bold');
        title(cb_lr, 'LR', 'FontSize', 8, 'FontWeight', 'bold');
        cb_edr = cb;
    else
        cb.Label.String = 'Connection strength (normalized, log scale)';
    end
end

set(gcf, 'Color', [1 1 1]);
sgtitle('Connectivity Matrices for EDR+LR');

% line the red colorbar up with the blue one after layout
drawnow;
cb_pos = cb_edr.Position;
cb_lr.Position = [cb_pos(1) + 2.2*cb_pos(3), cb_pos(2), cb_pos(3), cb_pos(4)];

exportgraphics(fig, fullfile('results', 'connectivity_matrices.png'), 'Resolution', 200);
fprintf('\nSaved results/connectivity_matrices.png\n');
