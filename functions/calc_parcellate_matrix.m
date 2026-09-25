function matrix_parcellated = calc_parcellate_matrix(parc, matrix_input)
% calc_parcellate_matrix.m
%
% PURPOSE: Averages a vertex-by-vertex matrix down to a parcel-by-parcel
%          matrix.
% ROLE:    Helper function. Used by steps 1, 3 and 5.
%
% Each output entry is the mean of the vertex entries between two parcels.
% Vertices with parc == 0 are ignored. The authors' code calls a function
% with this name but doesn't include it, so this is our own version. It
% uses sparse matrix products rather than loops so it's fast at vertex
% level.
%
% Inputs: parc         : parcellation labels aligned to matrix rows/cols [Nx1]
%         matrix_input : vertex-level matrix [NxN]
%
% Output: matrix_parcellated : parcellated matrix [num_parcels x num_parcels]

parcels = unique(parc(parc>0));
num_parcels = length(parcels);
num_vertices = size(matrix_input, 1);

parcel_of_vertex = zeros(num_vertices, 1);
for p = 1:num_parcels
    parcel_of_vertex(parc==parcels(p)) = p;
end
valid = parcel_of_vertex > 0;

% P(p,v) = 1 if vertex v belongs to parcel p, else 0
P = sparse(parcel_of_vertex(valid), find(valid), 1, num_parcels, num_vertices);

counts = full(sum(P, 2)); % number of vertices in each parcel
sum_matrix = P * double(matrix_input) * P';
count_matrix = counts * counts';

matrix_parcellated = full(sum_matrix) ./ count_matrix;
