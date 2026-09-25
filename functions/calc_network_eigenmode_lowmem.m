function [eig_vec, eig_val] = calc_network_eigenmode_lowmem(network, num_modes)
% calc_network_eigenmode_lowmem.m
%
% PURPOSE: Same eigenmodes as calc_network_eigenmode_dense.m, computed
%          with much less memory; for large (vertex-level) networks.
% ROLE:    Helper function. Used by steps 1-3.
%
% Same eigenmodes as calc_network_eigenmode_dense.m (normalized Laplacian,
% L_norm = D^(-1/2)*(D-A)*D^(-1/2)), but without building the Laplacian.
% At ~30,000 nodes a dense eig() runs MATLAB out of memory, and we only
% need a couple of hundred modes anyway.
%
% eigs() is given a function handle (Bmult) that only multiplies by
% `network`. It asks for the largest eigenvalues of
% M = D^(-1/2) A D^(-1/2), which is cheaper than the smallest of the
% Laplacian; the Laplacian eigenvalues are then 1 minus those.
%
% For small (parcel-level) networks the dense version is simpler.
%
% Inputs: network   : symmetric connectivity matrix [N x N]
%         num_modes : number of modes to return (int)
%
% Outputs: eig_vec  : eigenvectors (eigenmodes) [N x num_modes]
%          eig_val  : eigenvalues [num_modes x 1]

if nargin < 2
    num_modes = size(network,1);
end

N = size(network,1);

% remove self-connections
network(1:(N+1):end) = 0;

d = double(sum(network, 2));
dhalf = 1 ./ sqrt(d);

% isolated nodes (d = 0) would give 0/0; handle them separately
iso = (d == 0);
dhalf(iso) = 0;
iso_mask = double(iso);

% network is kept in single to save memory; the products are done in double
Bmult = @(x) iso_mask.*x + dhalf .* double(network * single(dhalf .* x));

opts.issym = true;
opts.isreal = true;
opts.p = min(N, 4*num_modes);   % larger subspace converges more reliably
opts.maxit = 1000;

[eig_vec, eig_val_M] = eigs(Bmult, N, num_modes, 'largestreal', opts);

eig_val_all = 1 - diag(eig_val_M);
[eig_val, sort_idx] = sort(eig_val_all, 'ascend');
eig_vec = single(eig_vec(:, sort_idx));
