% figure_fc_longrange_mask.m
%
% PURPOSE: Plots the functional long-range connections (r > 0.5, > 40 mm
%          apart) that Fig. 2 is scored on (the paper's Fig. 1 "LR
%          connections").
% ROLE:    Figure. Runs from the included results (made by step 4). Makes
%          results/fc_longrange_mask.png.
%
% Drawn in the style of the paper's Fig. 1B-iii/1D. Note this is the
% functional mask (from resting-state fMRI), not the structural
% exceptions spliced into EDR+LR (panel 4 of
% generate_connectivity_matrices.m). step6b saves the mask as an upper
% triangle.

s = load(fullfile('results', 'step6b_fc_longrange_reconstruction_fullres.mat'), ...
    'LR_mask_triu', 'n_LR_edges');
LR_mask = s.LR_mask_triu | s.LR_mask_triu'; % make it symmetric
num_parcels = size(LR_mask, 1);
fprintf('Functional LR connections: %d region pairs (%.2f%% of all pairs)\n', ...
    s.n_LR_edges, 100 * s.n_LR_edges / (num_parcels*(num_parcels-1)/2));

fig = figure('Name', 'Functional long-range connection mask', ...
    'Position', [100 100 520 480], 'Color', [1 1 1]);
imagesc(double(LR_mask)); axis square; colormap(parula); % same colours as the paper
set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 14);
xlabel('Regions'); ylabel('Regions');
title(sprintf('Functional LR connections (r > 0.5, > 40 mm; %d pairs)', s.n_LR_edges), 'FontSize', 11);

exportgraphics(fig, fullfile('results', 'fc_longrange_mask.png'), 'Resolution', 200);
fprintf('Saved results/fc_longrange_mask.png\n');
