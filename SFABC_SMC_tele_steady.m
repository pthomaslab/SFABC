function [n_opt,G_on,ks_t,Vs_t,ws_t,exact_M,exact_Ps] = SFABC_SMC_tele_steady(...
    data, d, tols, N, n_samples, rs, Min_c, max_N, wait_bar)
%SFABC_SMC_TELE_STEADY Run SFABC-SMC for the steady-state telegraph model.
%
%   [n_opt,G_on,ks_t,Vs_t,ws_t,exact_M,exact_Ps] = SFABC_SMC_tele_steady(...
%       data, d, tols, N, n_samples, rs, Min_c, max_N, wait_bar)
%
%   Inputs:
%     data      - 1 x N_obs row vector of observed mRNA counts for one gene.
%     d         - Highest raw moment order in the feasibility constraints.
%     tols      - Acceptance tolerances, length L.
%     N         - State-space truncation: G in {0,1}, X in {0,...,N}.
%     n_samples - Accepted particles per population.
%     rs        - Relative perturbation widths, length L-1.
%     Min_c     - Target successes per particle in volume estimation.
%     max_N     - Volume-estimation draw budgets, length L.
%     wait_bar  - Whether to display a progress bar.
%
%   Outputs:
%     n_opt     - L x 1 vector of feasibility optimization counts.
%     G_on      - L x n_samples array of gene ON probability masses.
%     ks_t      - L x n_samples x 3 array of accepted parameters
%                 [k_off, k_on, k_syn].
%     Vs_t      - L x n_samples array of relative volumes, mean 1 per row.
%     ws_t      - L x n_samples array of importance weights, scaled by 1/V.
%     exact_M   - L x n_samples x (d+1) array of moments over observed states.
%     exact_Ps  - L x n_samples cell array of feasible state probability vectors.
%
%   Notes:
%     Rates have independent log-uniform priors on [0.1, 10]; k_deg = 1.
%     Prior-specific blocks use (*); perturbation-kernel blocks use (**).

%% Initialize and preallocate
L = length(tols);
n_opt = zeros([L,1]);
G_on = zeros([L,n_samples]);
ks_t = zeros([L,n_samples,3]);
exact_M = zeros([L,n_samples,d+1]);
exact_Ps = cell(L,n_samples);
Vs_t = zeros([L,n_samples]);
ws_t = zeros([L,n_samples]);

% Get coefficient matrices for the telegraph model.
[Q_on,Q_off,Q_syn,Q_deg,I] = get_matrices(N);
% Extract unique observed states and their frequencies.
[counts,ns] = get_counts(data);

p = sdpvar(2*(N+1),1); % Set optimization variable
p_agg = I * p; % Aggregate joint state to marginal mRNA state
if wait_bar
    f = waitbar(0,'1');
end

%% Estimate moment uncertainty
% Estimate standard errors using 1000 Dirichlet samples.
m = zeros([1000, d+1]);
pt = zeros([1000,1]);
for i = 1:1000
    Ps = Dirichlet([ns 0]);
    for j = 0:d
        m(i,j+1) = sum(Ps(1:end-1).*counts.^j);
    end
    pt(i) = Ps(end);
end
se = std(m)';
% Apply the threefold standard-error multiplier to the zeroth moment.
se(1) = se(1)*3;

%% Run sequential Monte Carlo
l = 0;
while l < L
    n = 0;
    if l == 0
        while n < n_samples
            % Sample parameters from the prior (*).
            k_off = 10^unifrnd(-1,1);
            k_on = 10^unifrnd(-1,1);
            k_syn = 10^unifrnd(-1,1);
            k_deg = 1; % Fix degradation in steady-state

            % Draw candidate probabilities; the last bin covers unobserved states.
            Ps = Dirichlet([ns 0]);

            % Test feasibility under the steady-state CME constraints.
            feasibility = feas_test_tele_steady(counts, tols(l+1), se, p, Ps, ...
                                    k_off, k_on, k_syn, k_deg, ...
                                    Q_off, Q_on, Q_syn, Q_deg, I);

            n_opt(l+1) = n_opt(l+1) + 1;

            % Accept feasible parameters.
            if feasibility == 0
                % Store the parameters and feasible solutions.
                ks_t(l+1,n+1,:) = [k_off, k_on, k_syn];
                exact_M(l+1,n+1,:) = value(p_agg(counts+1).' * (counts(:).^ (0:d)));
                exact_Ps{l+1,n+1} = value(p);
                G_on(l+1,n+1) = value(sum(p(2:2:end))); % Probabilities of ON states (optional)

                n = n + 1;
                if wait_bar
                    waitbar(n/n_samples, f, "t="+sprintf('%d',l)+", n="+sprintf('%d',n))
                end
            end
        end

        % Estimate acceptance-region volumes for accepted particles.
        Vs_t(l+1,:) = V_estimation([counts 0], exact_M(l+1,:,:), ...
                                    se, tols(l+1), Min_c, max_N(l+1));

        % Set weights 1/V in the first population.
        ws_t(l+1,:) = 1 ./ Vs_t(l+1,:);

    else % l > 0
        r = rs(l); % Perturbation kernel parameter
        w_next = ones([1,n_samples]);
        while n < n_samples
            % Resample particles using the previous population's weights.
            ind = randsample(n_samples,1,true,ws_t(l,:));

            % Apply the relative perturbation kernel (**).
            k_off = ks_t(l,ind,1)*(1+unifrnd(-r,r));
            k_on = ks_t(l,ind,2)*(1+unifrnd(-r,r));
            k_syn = ks_t(l,ind,3)*(1+unifrnd(-r,r));

            % Reject parameters outside the prior support (*).
            if log10(min([k_off,k_on,k_syn])) < -1 || log10(max([k_off,k_on,k_syn])) > 1
                continue
            end
            k_deg = 1; % Fix degradation in steady-state

            % Draw candidate probabilities; the last bin covers unobserved states.
            Ps = Dirichlet([ns 0]);

            % Test feasibility under the steady-state CME constraints.
            feasibility = feas_test_tele_steady(counts, tols(l+1), se, p, Ps, ...
                                    k_off, k_on, k_syn, k_deg, ...
                                    Q_off, Q_on, Q_syn, Q_deg, I);

            n_opt(l+1) = n_opt(l+1) + 1;

            % Accept feasible parameters.
            if feasibility == 0
                % Store the parameters and feasible solutions.
                ks_t(l+1,n+1,:) = [k_off, k_on, k_syn];
                exact_M(l+1,n+1,:) = value(p_agg(counts+1).' * (counts(:).^ (0:d)));
                exact_Ps{l+1,n+1} = value(p);
                G_on(l+1,n+1) = value(sum(p(2:2:end))); % Probabilities of ON states (optional)

                % Calculate importance weights.
                c = 0;
                for i = 1:n_samples
                    if abs(k_off/ks_t(l,i,1)-1) <= r && abs(k_on/ks_t(l,i,2)-1) <= r && abs(k_syn/ks_t(l,i,3)-1) <= r
                        c = c + ws_t(l,i)/ks_t(l,i,1)/ks_t(l,i,2)/ks_t(l,i,3); % Specific to perturbation kernel**
                    end
                end
                w_next(n+1) = 1/k_off/k_on/k_syn/c; % Specific to prior*

                n = n + 1;
                if wait_bar
                    waitbar(n/n_samples, f, "t="+sprintf('%d',l)+", n="+sprintf('%d',n))
                end
            end
        end
        % Estimate acceptance-region volumes for accepted particles.
        Vs_t(l+1,:) = V_estimation([counts 0], exact_M(l+1,:,:), ...
                                    se, tols(l+1),Min_c, max_N(l+1));

        % Scale weights by 1/V.
        ws_t(l+1,:) = w_next ./ Vs_t(l+1,:);
    end
    l = l + 1; % Move to next population
end
if wait_bar
    delete(f)
end
end

