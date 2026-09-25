% figure_3D_surface_reconstruction.m
%
% PURPOSE: Brain-surface renderings of a motor task map rebuilt from 5-20
%          Geometry vs. EDR+LR modes (the paper's Fig. 3C/D).
% ROLE:    Figure. Needs
%          results/synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat, which
%          step 1 has not produced yet.
%
% Adapted from the authors' PNAS_Figure_3_main_part_CD.m. Their version
% loads a reconstruction file that isn't in the release, so this one
% reconstructs the map directly. Also plots EDR+LR mode 16, which the
% paper highlights.

%% Setup

eigenmode_toolbox_dir = fullfile('pang2023_BrainEigenmodes', 'BrainEigenmodes-main');
edrlr_data_dir = fullfile('Vohryzek2024', 'vohryzek2024_EDRLR');

addpath(genpath(fullfile(eigenmode_toolbox_dir, 'functions_matlab')));
addpath('functions');

hemisphere = 'lh';
num_modes = 200;
surface_interest = 'fsLR_32k';
mesh_interest = 'midthickness';

[vertices, faces] = read_vtk(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_%s-%s.vtk', surface_interest, mesh_interest, hemisphere)));
surface_midthickness.vertices = vertices';
surface_midthickness.faces = faces';

cortex = dlmread(fullfile(edrlr_data_dir, 'Data', 'template_surfaces', ...
    sprintf('%s_cortex-%s_mask.txt', surface_interest, hemisphere)));
cortex_ind = find(cortex);
medial_wall = find(cortex==0);

%% Activation map (motor, right foot vs. average)

data = load(fullfile(edrlr_data_dir, 'Data', 'empirical', 'S255_tfMRI_ALLTASKS_raw_lh.mat'), 'zstat');
activation_map_motor_rf_avg = mean(data.zstat.motor_rf_avg, 2, 'omitnan');
clear data

%% Fig. 3C: activation map

fig = draw_surface_bluewhitered_gallery_dull(surface_midthickness, activation_map_motor_rf_avg, ...
    hemisphere, medial_wall, 1);
fig.Name = 'Figure 3C - activation map (motor rf_avg)';

%% Fig. 3C: EDR+LR mode 16

edrlr = load(fullfile('results', 'synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat'), 'eig_vec_EDRLR');
eig_vec_EDRLR = edrlr.eig_vec_EDRLR;
clear edrlr

fig = draw_surface_bluewhitered_gallery_dull(surface_midthickness, eig_vec_EDRLR(:, 16), ...
    hemisphere, medial_wall, 1);
fig.Name = 'Figure 3C - EDR+LR mode 16';

%% Geometry eigenmodes (Pang's template, 200 modes)

eig_vec_geometry = dlmread(fullfile(eigenmode_toolbox_dir, 'data', 'template_eigenmodes', ...
    sprintf('%s_%s-%s_emode_%i.txt', surface_interest, mesh_interest, hemisphere, num_modes)));

%% Reconstruct the map from each basis with increasing N

mode_interest = [5, 10, 15, 20];

recon_geometry = zeros(length(cortex), numel(mode_interest));
recon_edrlr = zeros(length(cortex), numel(mode_interest));
for m = 1:numel(mode_interest)
    N = mode_interest(m);

    basis_g = eig_vec_geometry(cortex_ind, 1:N);
    beta_g = calc_eigendecomposition(activation_map_motor_rf_avg(cortex_ind), basis_g, 'matrix');
    recon_geometry(:, m) = eig_vec_geometry(:, 1:N) * beta_g;

    basis_e = eig_vec_EDRLR(cortex_ind, 1:N);
    beta_e = calc_eigendecomposition(activation_map_motor_rf_avg(cortex_ind), basis_e, 'matrix');
    recon_edrlr(:, m) = eig_vec_EDRLR(:, 1:N) * beta_e;
end

%% Fig. 3D: Geometry, N = 20, 15, 10, 5

fig = draw_surface_bluewhitered_gallery_dull(surface_midthickness, recon_geometry, hemisphere, medial_wall, 1);
fig.Name = 'Figure 3D - Geometry reconstruction (modes = 5,10,15,20)';

%% Fig. 3D: EDR+LR, N = 20, 15, 10, 5

fig = draw_surface_bluewhitered_gallery_dull(surface_midthickness, recon_edrlr, hemisphere, medial_wall, 1);
fig.Name = 'Figure 3D - EDR+LR reconstruction (modes = 5,10,15,20)';

%% Correlation vs. number of modes

corr_geometry = zeros(1, num_modes);
corr_edrlr = zeros(1, num_modes);
for N = 1:num_modes
    basis_g = eig_vec_geometry(cortex_ind, 1:N);
    beta_g = calc_eigendecomposition(activation_map_motor_rf_avg(cortex_ind), basis_g, 'matrix');
    corr_geometry(N) = corr(activation_map_motor_rf_avg(cortex_ind), eig_vec_geometry(cortex_ind,1:N)*beta_g);

    basis_e = eig_vec_EDRLR(cortex_ind, 1:N);
    beta_e = calc_eigendecomposition(activation_map_motor_rf_avg(cortex_ind), basis_e, 'matrix');
    corr_edrlr(N) = corr(activation_map_motor_rf_avg(cortex_ind), eig_vec_EDRLR(cortex_ind,1:N)*beta_e);
end

figure('Name', 'Figure 3D - reconstruction accuracy vs modes');
plot(1:num_modes, corr_geometry, 'm-', 'linewidth', 2, 'DisplayName', 'Geometry'); hold on
plot(1:num_modes, corr_edrlr, 'k-', 'linewidth', 2, 'DisplayName', 'EDR+LR');
xlabel('number of modes'); ylabel('reconstruction accuracy (corr)')
legend('location', 'southeast'); grid on
title('motor rf-avg activation map reconstruction: Geometry vs EDR+LR')
