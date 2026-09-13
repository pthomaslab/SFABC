function [Q_on, Q_off, Q_syn, Q_deg, I] = get_matrices(N)
%GET_MATRICES Construct coefficient matrices for the telegraph-model CME.
%
%   [Q_on, Q_off, Q_syn, Q_deg, I] = get_matrices(N) builds the coefficient
%   matrices (size 2N x (2N+2)) for gene switching, transcription, and
%   degradation. I has size (N+1) x (2N+2) and gives marginal mRNA probabilities.
%
%   Notes:
%     State order: [OFF_0, ON_0, OFF_1, ON_1, ...].
%     Even/odd indices below refer to the zero-based coefficient index x.

M = 2*N;

%% Q_on: gene activation
% Transition from OFF to ON (index shift +1).
% Main diagonal: -1 at even indices (outflow from OFF).
diag_on_main = arrayfun(@(x) -1 * (~mod(x, 2)), 0:M+1);
% Subdiagonal (-1): +1 at even indices (inflow to ON).
diag_on_sub  = arrayfun(@(x)  1 * (~mod(x, 2)), 0:M);

Q_on = diag(diag_on_main, 0) + diag(diag_on_sub, -1);
Q_on = Q_on(1:end-2, :);

%% Q_off: gene inactivation
% Transition from ON to OFF (index shift -1).
% Main diagonal: -1 at odd indices (outflow from ON).
diag_off_main = arrayfun(@(x) -1 * mod(x, 2), 0:M+1);
% Superdiagonal (+1): +1 at even indices (inflow to OFF).
diag_off_sub  = arrayfun(@(x)  1 * (~mod(x, 2)), 0:M);

Q_off = diag(diag_off_main, 0) + diag(diag_off_sub, 1);
Q_off = Q_off(1:end-2, :);

%% Q_syn: transcription
% Production of mRNA (index shift +2).
% Main diagonal: -1 at odd indices (outflow from the ON state).
diag_syn_main = arrayfun(@(x) -1 * mod(x, 2), 0:M+1);
% Second subdiagonal (-2): +1 at odd indices (inflow to the higher count).
diag_syn_sub  = arrayfun(@(x)  1 * mod(x, 2), 0:M-1);

Q_syn = diag(diag_syn_main, 0) + diag(diag_syn_sub, -2);
Q_syn = Q_syn(1:end-2, :);

%% Q_deg: degradation
% Decay of mRNA (index shift -2).
% Main diagonal: -floor(x/2), proportional to the mRNA count.
diag_deg_main = arrayfun(@(x) -floor(x/2) * (x > 1), 0:M+1);
% Second superdiagonal (+2): floor(x/2)+1, the source state's mRNA count.
diag_deg_sub  = arrayfun(@(x)  floor(x/2) + 1, 0:M-1);

Q_deg = diag(diag_deg_main, 0) + diag(diag_deg_sub, 2);
Q_deg = Q_deg(1:end-2, :);

%% I: probability aggregation
% Sum each OFF/ON pair to obtain marginal mRNA probabilities.
I = zeros(M/2 + 1, M + 2);
for i = 1:(M/2 + 1)
    % Select columns (2i-1, 2i), corresponding to mRNA count i-1.
    I(i, 2*i - 1 : 2*i) = ones(1, 2);
end
end
