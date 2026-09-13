function V = V_estimation(ns, exact_M, se, tol, Min_c, max_N)
%V_ESTIMATION Estimate relative acceptance-region volumes by Monte Carlo.
%
%   V = V_estimation(ns, exact_M, se, tol, Min_c, max_N)
%
%   Inputs:
%     ns        - B observed state values followed by a zero placeholder
%                 for unobserved states; these are not occurrence frequencies.
%     exact_M   - 1 x n_samples x (d+1) array of moments (orders 0,...,d)
%                 from feasible CME solutions.
%     se        - Standard errors for moments of orders 0,...,d.
%     tol       - Tolerance multiplier for the moment intervals.
%     Min_c     - Target successes per particle before stopping.
%     max_N     - Draw budget, checked between batches of 1000 draws.
%
%   Output:
%     V         - 1 x n_samples vector of volume estimates, scaled to mean 1.

n_samples = size(exact_M, 2);
d         = size(exact_M, 3) - 1;
K = length(ns);
ns = ns(:);

%% Build acceptance intervals
intvs_m   = zeros([n_samples, d+1, 2]);
for i = 1:d+1
    m_i = reshape(exact_M(1,:,i), [n_samples, 1]);
    intvs_m(:,i,1) = m_i - se(i) * tol;
    intvs_m(:,i,2) = m_i + se(i) * tol;
end

%% Estimate volumes by Monte Carlo
n_succ     = zeros([1, n_samples]);
n_rep      = 0;
batch_size = 1e3;
while min(n_succ) < Min_c && n_rep < max_N
    % Draw uniformly on the simplex: Dirichlet(1,...,1).
    U = exprnd(1, batch_size, K);
    n_rep = n_rep + batch_size;
    P = U ./ sum(U, 2);

    % Calculate moments over the observed states.
    M = ns .^ (0:d);
    M(end, 1) = 0; % Remove contribution from last entry (unobserved prob)
    X = P * M;

    % Check the moment intervals for all particles.
    succ = true(batch_size, n_samples);
    for i = 1:d+1
        lower = (intvs_m(:,i,1))';
        upper = (intvs_m(:,i,2))';
        succ = succ & (X(:,i) >= lower) & (X(:,i) <= upper);
    end
    n_succ = n_succ + sum(succ, 1); % Accumulate successes
end

%% Normalize volume estimates
V = max(n_succ, 1) / n_rep; % Avoid 0 volume estimate
V = reshape(V, [1, n_samples]);
V = V / mean(V); % Normalize volumes
end
