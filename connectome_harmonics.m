% connectome_harmonics.m
%
% PURPOSE: Builds the EDR+LR structural connectome (exponential-distance-
%          rule fit plus the spliced-in long-range exceptions) and
%          computes its eigenmodes.
% ROLE:    Pipeline step 1. Needs the Pang + Vohryzek data folders. Makes
%          results/EDR_LR_connectome_full_resolution.mat.
%
% Method (Vohryzek et al., 2025):
%   1. Fit an exponential distance rule (EDR) to the structural
%      connectome: strength falls off exponentially with distance.
%   2. Find the long-range (LR) exceptions: connections > 40 mm long
%      that are more than 3 SD stronger than typical for their distance.
%   3. Put those exceptions back on top of the EDR baseline (EDR+LR).
%   4. Compute the eigenmodes of EDR+LR.
%
% Stage 1 runs the whole thing at parcel level (180 regions) as a quick
% check. Stage 2 is the real version at full vertex resolution (29,696
% vertices), and its output is what the rest of the pipeline uses.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

%% Stage 1: load data (Glasser360, left hemisphere)

hemisphere = 'lh';
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';
parc_name = 'Glasser360';

% surface geometry, used to get parcel centroid distances
[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

% cortex mask (which vertices are actually cortex vs. medial wall)
cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
num_vertices = length(cortex);

% parcellation labels (given at full-vertex resolution; keep only cortex vertices)
parc = dlmread(fullfile(edrlr_data_dir, 'Data', 'parcellations', ...
    sprintf('%s_%s-%s.txt', surface_interest, parc_name, hemisphere)));
parc_cortex = parc(cortex_ind);
parcels = unique(parc_cortex(parc_cortex>0));
num_parcels = length(parcels);
fprintf('Loaded %s parcellation: %d parcels\n', parc_name, num_parcels);

% group-average structural connectome (29,696 x 29,696). The copy in
% Vohryzek2024/ is truncated, so this loads Pang's copy.
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');

%% Stage 1: parcel centroid distances and parcellated connectome

vertices_cortex = surface_midthickness.vertices(cortex_ind, :);
centroids = zeros(num_parcels, 3);
for p = 1:num_parcels
    centroids(p,:) = mean(vertices_cortex(parc_cortex==parcels(p), :), 1);
end
rr_parc = squareform(pdist(centroids)); % [num_parcels x num_parcels] Euclidean distance

connectome_parc = calc_parcellate_matrix(parc_cortex, avgSC_L);
C_parc = connectome_parc / max(connectome_parc(:));
clear avgSC_L

%% Stage 1: fit the EDR, C(i,j) = A * exp(-lambda * d(i,j))
% Connections are binned by distance and the curve is fit to the mean
% strength in each bin.

NR = 60;
NSTD = 3;       % SDs above the bin mean to count as an exception (paper: 3)
DistRange = 40; % minimum length of an LR connection, mm (paper: 40)

range_dist = max(rr_parc(:));
delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);

index_parc = floor(rr_parc/delta) + 1;
index_parc(index_parc > NR) = NR;

sc_density = cell(1, NR);
sc_density_i = cell(1, NR);
sc_density_j = cell(1, NR);
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

% skip the shortest-distance bins and any empty ones
fit_start = find(xcoor >= 10, 1);
fit_ind = fit_start:NR;
fit_ind = fit_ind(~isnan(ycoor2(fit_ind)));

% the authors use lsqcurvefit; used fminsearch instead because I did not want 
% to download  the Optimization Toolbox
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - ycoor2(fit_ind)).^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
A0 = [0.15, 0.18];
Afit = fminsearch(sse, A0, options);
lambda = Afit(2);
yl = Afit(1)*exp(-Afit(2)*xcoor);

fprintf('EDR fit (Glasser360, validation only): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

%% Stage 1: long-range exceptions and the EDR+LR connectome
% EDR+LR is the EDR curve with the real strength of each exception put back.

Clong = zeros(num_parcels, num_parcels);     % LR exceptions (> 40 mm)
Clong_all = zeros(num_parcels, num_parcels); % all > 3 SD outliers, any distance

for i = fit_start:NR
    if isempty(sc_density{i})
        continue
    end
    mv = mean(sc_density{i});
    st = std(sc_density{i});
    ind_exc = find(sc_density{i} > mv + NSTD*st);
    for n = 1:numel(ind_exc)
        ii = sc_density_i{i}(ind_exc(n));
        jj = sc_density_j{i}(ind_exc(n));
        Clong_all(ii,jj) = sc_density{i}(ind_exc(n));
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

%% Stage 1: plots

figure('Name', 'Stage 1 validation - EDR fit');
errorbar(xcoor, cellfun(@(x) mean(x,'omitnan'), sc_density), cellfun(@(x) std(x,'omitnan'), sc_density), 'o');
hold on
plot(xcoor, yl, 'r-', 'linewidth', 2)
xlabel('Distance (mm)'); ylabel('Connection strength (normalized)')
legend('binned SC (mean +/- SD)', 'EDR fit')
title(sprintf('Glasser360 validation: lambda = %.4f /mm', lambda))
grid on

figure('Name', 'Stage 1 validation - Connectome comparison');
subplot(2,2,1); imagesc(C_parc); axis square; colorbar; title('Connectome (parcellated)')
subplot(2,2,2); imagesc(EDR_conn); axis square; colorbar; title('EDR fit')
subplot(2,2,3); imagesc(EDR_LRE_parc); axis square; colorbar; title('EDR+LR')
subplot(2,2,4); imagesc(Clong>0); axis square; colorbar; title(sprintf('Long-range exceptions (n=%d)', nnz(triu(Clong,1))))
colormap(flipud(bone))

clear Clong Clong_all index_parc sc_density sc_density_i sc_density_j

%% Stage 2: EDR+LR at full vertex resolution
%
% Same steps as Stage 1, on all 29,696 cortical vertices.

if ~exist('results', 'dir'); mkdir('results'); end

num_modes = 200;
NR = 400; NRini = 20; NRfin = 380; NSTD = 3; DistRange = 40; % as in the original paper

fprintf('\n=== STAGE 2: full vertex-resolution EDR+LR derivation ===\n');

% empirical connectome (Pang's copy, single precision, 29696x29696)
load(fullfile(eigenmode_toolbox_dir, 'data', 'empirical', ...
    'S255_high-resolution_group_average_connectome_cortex_nomedial-lh.mat'), 'avgSC_L');
C = avgSC_L / max(avgSC_L(:));
clear avgSC_L

% vertex-level Euclidean distances
tic
rr = squareform(pdist(single(vertices_cortex)));
fprintf('Computed %dx%d vertex distance matrix in %.1f s\n', size(rr,1), size(rr,2), toc);

% bin by distance with accumarray (a per-bin loop is far too slow here)
range_dist = max(rr(:));
delta = range_dist / NR;
xcoor = delta/2 + delta*(0:NR-1);

index = uint16(min(floor(double(rr)/delta) + 1, NR));

tic
bin_mean = accumarray(index(:), C(:), [NR,1], @mean, single(NaN));
bin_std  = accumarray(index(:), C(:), [NR,1], @std, single(NaN));
fprintf('Binned connectivity by distance (%d bins) in %.1f s\n', NR, toc);

fit_ind = 25:NR; % the authors leave out the 24 shortest-distance bins
fit_ind = fit_ind(~isnan(bin_mean(fit_ind)));
expfunc = @(A, x) (A(1)*exp(-A(2)*x));
sse = @(A) sum((expfunc(A, xcoor(fit_ind)) - double(bin_mean(fit_ind))').^2);
options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');
Afit = fminsearch(sse, [0.15, 0.18], options);
lambda = Afit(2);
fprintf('EDR fit (full resolution): A = %.4f, lambda = %.4f /mm\n', Afit(1), lambda);

% long-range exceptions: > 3 SD above their bin mean and > 40 mm apart
tic
local_mean = bin_mean(index);
local_std = bin_std(index);
is_exception = (index >= NRini) & (index <= NRfin) & (rr > DistRange) & (C > local_mean + NSTD*local_std);
clear local_mean local_std
fprintf('Identified %d long-range exceptions in %.1f s\n', nnz(triu(is_exception,1)), toc);

EDR_LRE = single(Afit(1) * exp(-Afit(2) * double(rr)));
EDR_LRE(is_exception) = C(is_exception);
clear is_exception index C rr

save(fullfile('results', 'EDR_LR_connectome_full_resolution.mat'), 'EDR_LRE', 'Afit', 'lambda', '-v7.3');
fprintf('Saved full-resolution EDR+LR connectome to results/EDR_LR_connectome_full_resolution.mat\n');

%% Stage 2: EDR+LR eigenmodes

fprintf('\nComputing EDR+LR eigenmodes (dense eig on %dx%d)...\n', size(EDR_LRE,1), size(EDR_LRE,2));
tic
[eig_vec_temp, eig_val] = calc_network_eigenmode_lowmem(EDR_LRE, num_modes);
fprintf('Eigendecomposition done in %.1f s\n', toc);
clear EDR_LRE

% put the medial wall back in as zeros (same layout as Pang's files)
eig_vec_EDRLR = zeros(num_vertices, num_modes, 'single');
eig_vec_EDRLR(cortex_ind, :) = eig_vec_temp;
clear eig_vec_temp

save(fullfile('results', 'synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat'), ...
    'eig_vec_EDRLR', 'eig_val', '-v7.3');
fprintf('Saved EDR+LR eigenmodes to results/synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat\n');

%% Stage 2: plot EDR+LR eigenmodes on the cortical surface

mode_interest = [2, 3, 4, 5, 7, 16]; % modes shown in the paper's Fig. 3
surface_to_plot = surface_midthickness;
data_to_plot = eig_vec_EDRLR(:, mode_interest);
medial_wall = find(cortex==0);
with_medial = 1;

fig = draw_surface_bluewhitered_gallery_dull(surface_to_plot, data_to_plot, hemisphere, medial_wall, with_medial);
fig.Name = 'Stage 2 - EDR+LR eigenmodes';
