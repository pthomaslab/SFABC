function [us, ns] = get_counts(data)
%GET_COUNTS Extract unique count vectors and their frequencies.
%
%   [us, ns] = get_counts(data) identifies unique columns (us) in the
%   input data matrix and counts how many times each appears (ns).
%   data is S x N_obs: each column is one observation of S species.
%   us is S x B and ns is 1 x B, where B is the number of unique observations.

% Find unique observations and their indices using transposed data.
[us, ~, idx] = unique(data', 'rows','sorted');

% Count the frequency of each unique observation.
ns = accumarray(idx, 1).';

% Return unique observations as columns.
us = us.';
end
