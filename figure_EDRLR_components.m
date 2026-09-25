% figure_EDRLR_components.m
%
% PURPOSE: Two-colour EDR+LR matrix: blue where a region pair's wiring
%          comes mainly from the EDR baseline, red where it comes mainly
%          from the long-range splice.
% ROLE:    Figure. Runs from the included results (made by step 5). Makes
%          results/connectivity_EDRLR_components.png.
%
% Darker means stronger, on one log scale for both colours.
%
% Each parcel pair averages many vertex pairs, some from the EDR baseline
% and some from the splice, so the splice's share of a pair is
% exception_strength_parc ./ EDRLR_parc. In practice pairs are almost all
% one or the other (median share ~0.03, upper quartile ~0.95), so each
% pair is coloured by whichever part is larger.

load(fullfile('results', 'connectivity_matrices_parc.mat'), ...
    'EDRLR_parc', 'exception_strength_parc', 'lambda');

num_parcels = size(EDRLR_parc, 1);
offdiag = ~eye(num_parcels);

splice_share = exception_strength_parc ./ EDRLR_parc;
is_splice = (splice_share > 0.5) & offdiag;
fprintf('Parcel pairs dominated by the long-range splice: %.1f%%\n', ...
    100 * nnz(is_splice) / nnz(offdiag));

%% Map strength onto a shared log scale, then onto blue or red

logv = log10(double(EDRLR_parc));
clim_log = [prctile(logv(offdiag), 1), max(logv(offdiag))];
t = (logv - clim_log(1)) / diff(clim_log);
t = min(max(t, 0), 1);

n_col = 256;
light = [0.94 0.95 0.97];
cmap_edr = interp1([0 1], [light; 0.03 0.15 0.40], linspace(0, 1, n_col)); % light -> navy
cmap_lr  = interp1([0 1], [light; 0.55 0.02 0.05], linspace(0, 1, n_col)); % light -> dark red

idx = round(t * (n_col - 1)) + 1;
rgb = zeros(num_parcels, num_parcels, 3);
for c = 1:3
    ch_edr = cmap_edr(:, c);
    ch_lr = cmap_lr(:, c);
    layer = ch_edr(idx);
    layer(is_splice) = ch_lr(idx(is_splice));
    layer(~offdiag) = 1; % blank the diagonal
    rgb(:, :, c) = layer;
end

%% Plot

fig = figure('Name', 'EDR+LR components', 'Position', [100 100 900 720], 'Color', [1 1 1]);

ax = axes('Position', [0.08 0.10 0.66 0.80]);
image(ax, rgb); axis(ax, 'square');
xlabel(ax, 'Parcel'); ylabel(ax, 'Parcel');
title(ax, sprintf('EDR+LR connectome: EDR baseline vs. long-range splice (\\lambda = %.3f /mm)', lambda));

% two colorbars drawn by hand, same log scale
tick_exp = ceil(clim_log(1)):floor(clim_log(2));
tick_pos = (tick_exp - clim_log(1)) / diff(clim_log) * (n_col - 1) + 1;
tick_lbl = arrayfun(@(e) sprintf('10^{%d}', e), tick_exp, 'UniformOutput', false);

cb_specs = {cmap_edr, {'EDR', 'baseline'}, 0.78; cmap_lr, {'Long-range', 'splice'}, 0.88};
for k = 1:size(cb_specs, 1)
    cax = axes('Position', [cb_specs{k,3} 0.20 0.03 0.60]);
    image(cax, permute(cb_specs{k,1}, [1 3 2]));
    set(cax, 'YDir', 'normal', 'XTick', [], 'YTick', tick_pos, 'YTickLabel', tick_lbl, ...
        'YAxisLocation', 'right', 'TickLabelInterpreter', 'tex');
    title(cax, cb_specs{k,2}, 'FontSize', 9);
end
annotation(fig, 'textbox', [0.76 0.11 0.2 0.06], 'String', {'Connection strength', '(normalized, log scale)'}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);

exportgraphics(fig, fullfile('results', 'connectivity_EDRLR_components.png'), 'Resolution', 200);
fprintf('Saved results/connectivity_EDRLR_components.png\n');
