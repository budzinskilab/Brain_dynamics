function [eig_vec, eig_val] = calc_network_eigenmode_dense(network, num_modes)
% calc_network_eigenmode_dense.m
%
% PURPOSE: Eigenmodes of a network's normalized graph Laplacian using
%          full eig(); for small (parcel-level) networks.
% ROLE:    Helper function. Used by step 3.
%
% L_norm = D^(-1/2) * (D - A) * D^(-1/2), as in Vohryzek et al. (2025).
% Same result as calc_network_eigenmode_lowmem.m. Use this one for a few
% hundred nodes; at vertex level (~30,000 nodes) it runs out of memory.
%
% Inputs: network   : symmetric connectivity matrix [N x N]
%         num_modes : number of modes to return (int)
%
% Outputs: eig_vec  : eigenvectors (eigenmodes) [N x num_modes], ascending
%                     eigenvalue order
%          eig_val  : eigenvalues [num_modes x 1], ascending

if nargin < 2
    num_modes = size(network,1);
end

N = size(network,1);

network(1:(N+1):end) = 0; % remove self-connections

d = sum(double(network), 2);
dhalf = 1 ./ sqrt(d);

Lnorm = -(dhalf .* double(network) .* dhalf');
Lnorm(1:(N+1):end) = Lnorm(1:(N+1):end) + 1;
Lnorm(isnan(Lnorm)) = 0; % isolated nodes (d=0) produce 0/0 = NaN

[V, D] = eig(Lnorm);

% eig() doesn't always return them sorted
eig_val_all = diag(D);
[eig_val_all, sort_idx] = sort(eig_val_all, 'ascend');
V = V(:, sort_idx);

eig_vec = V(:, 1:num_modes);
eig_val = eig_val_all(1:num_modes);
