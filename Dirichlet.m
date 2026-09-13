function Ps = Dirichlet(ns)
%DIRICHLET Sample a Dirichlet probability vector from observed frequencies.
%
%   Ps = Dirichlet(ns) draws from Dirichlet(alpha), where alpha = ns + 1
%   and ns is a vector of occurrence frequencies.

% Calculate the Dirichlet parameters.
alpha = ns + 1;

% Draw independent gamma samples with shape alpha and scale 1.
ys = gamrnd(alpha, ones(size(alpha)));

% Normalize to obtain probabilities summing to 1.
Ps = ys ./ sum(ys);
end
