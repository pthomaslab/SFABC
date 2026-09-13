function [n_sim,ks_t,ws_t] = ABC_SMC_tele_steady(data, d, tols, n_samples, rs, wait_bar, gillespie)

%ABC_SMC_TELE_STEADY Run ABC-SMC inference for the telegraph model.
%
%   [n_sim,ks_t,ws_t] = ABC_SMC_tele_steady(data, d, tols, n_samples, rs, wait_bar, gillespie)
%
%   Inputs:
%     data      - Vector of observed mRNA counts for one gene.
%     d         - Highest raw moment order in the distance function.
%     tols      - Acceptance tolerances, length L.
%     n_samples - Accepted particles per population.
%     rs        - Relative perturbation widths, length L-1.
%     wait_bar  - Whether to display a progress bar.
%     gillespie - true: Gillespie SSA to t=20; false: steady-state sampling.
%
%   Outputs:
%     n_sim     - L x 1 vector of simulated dataset counts.
%     ks_t      - L x n_samples x 3 array: [k_off, k_on, k_tx].
%     ws_t      - L x n_samples array of unnormalized importance weights.
%
%   k_tx denotes k_syn elsewhere in the project; k_deg is fixed at 1.
%   Rates have independent log-uniform priors on [0.1, 10].

t_max = length(tols);
se = get_se(data,d); n_obs = length(data);
n_sim = zeros([length(tols),1]);
ks_t = zeros([length(tols),n_samples,3]);
ws_t = zeros([length(tols),n_samples]);
if gillespie
    pfun = @propensities_2state;
    stoich_matrix = [ -1  1  0;
        1  -1  0;
        0   0  1;
        0   0  -1];
    tspan = [0 20];
    x0 = [1, 0, 0];
end
t = 0;
if wait_bar
    f = waitbar(0,'1');
end
while t < t_max
    n = 0;
    if t == 0
        while n < n_samples

            % Sample parameters from the prior.
            k_off = 10^unifrnd(-1,1);
            k_on = 10^unifrnd(-1,1);
            k_tx = 10^unifrnd(-1,1);
            k_deg = 1;

            if gillespie
                x_sim = zeros([n_obs, 1]);
                for row = 1:n_obs
                    [~,x] = directMethod(stoich_matrix, pfun, tspan,...
                        x0, k_off, k_on, k_tx, k_deg);
                    x_sim(row) = x(end,3);
                end
            else
                p = betarnd(k_on,k_off,[n_obs,1]);
                x_sim = poissrnd(k_tx.*p);
            end
            n_sim(t+1) = n_sim(t+1) + 1;

            distance = moment_distance(data, x_sim, d, se);
            if distance <= tols(t+1)
                ks_t(t+1,n+1,:) = [k_off, k_on, k_tx];

                n = n + 1;
                if wait_bar
                    waitbar(n/n_samples,f,...
                        "t="+sprintf('%d',t)+", n="+sprintf('%d',n))
                end
            end
        end
        ws_t(t+1,:) = 1;
    else
        r = rs(t);
        w_next = ones([1,n_samples]);
        while n < n_samples

            % Resample and perturb parameters from the previous population.
            ind = randsample(n_samples,1,true,ws_t(t,:));
            k_off = ks_t(t,ind,1)*(1+unifrnd(-r,r));
            k_on = ks_t(t,ind,2)*(1+unifrnd(-r,r));
            k_tx = ks_t(t,ind,3)*(1+unifrnd(-r,r));
            if log10(min([k_off,k_on,k_tx])) < -1 || log10(max([k_off,k_on,k_tx])) > 1
                % Reject parameters outside the prior support.
                continue
            end
            k_deg = 1;

            if gillespie
                x_sim = zeros([n_obs, 1]);
                for row = 1:n_obs
                    [~,x] = directMethod(stoich_matrix, pfun, tspan,...
                        x0, k_off, k_on, k_tx, k_deg);
                    x_sim(row) = x(end,3);
                end
            else
                p = betarnd(k_on,k_off,[n_obs,1]);
                x_sim = poissrnd(k_tx.*p);
            end
            n_sim(t+1) = n_sim(t+1) + 1;

            distance = moment_distance(data, x_sim, d, se);
            if distance <= tols(t+1)
                ks_t(t+1,n+1,:) = [k_off, k_on, k_tx];
                c = 0;
                for i = 1:n_samples
                    if abs(k_off/ks_t(t,i,1)-1) <= r && abs(k_on/ks_t(t,i,2)-1) <= r && abs(k_tx/ks_t(t,i,3)-1) <= r
                        c = c + ws_t(t,i)/ks_t(t,i,1)/ks_t(t,i,2)/ks_t(t,i,3);
                    end
                end
                w_next(n+1) = 1/k_off/k_on/k_tx/c;
                n = n + 1;
                if wait_bar
                    waitbar(n/n_samples,f,...
                        "t="+sprintf('%d',t)+", n="+sprintf('%d',n))
                end
            end
        end
        ws_t(t+1,:) = w_next;
    end
    t = t + 1;
end
if wait_bar
    delete(f)
end
end

function se = get_se(data, d)

%GET_SE Estimate raw moment standard errors by bootstrap resampling.
% Preallocate one standard error per requested moment.
se = zeros(d, 1);

% Loop over moment orders j = 1,...,d.
for j = 1:d
    % Bootstrap the j-th raw moment mean(x.^j) using 1000 resamples.
    means = bootstrp(1000, @(x) mean(x.^j), data);

    % The standard error is the standard deviation of bootstrap replicates.
    se(j) = sqrt(var(means));
end
end

function [distance] = moment_distance(x_obs, x_sim, d, se)
%MOMENT_DISTANCE Compute the maximum standardized raw moment difference.
% Compare orders 1,...,d; x_obs and x_sim may have different lengths.
distance = 0;
for i = 1:d
    distance = max(abs(mean(x_obs.^i) - mean(x_sim.^i)) / se(i),distance);
end
end


function a = propensities_2state(x, k1, k2, k3, k4)
%PROPENSITIES_2STATE Return reaction propensities for the telegraph model.
% x = [ON, OFF, mRNA]; rates are [k_off, k_on, k_syn, k_deg].
p1  = x(1);
p2 =  x(2);
p3 =  x(3);
a = [k1*p1;
    k2*p2;
    k3*p1;
    k4*p3];
end

function [ t, x ] = directMethod( stoich_matrix, propensity_fcn, tspan, x0,...
    k1, k2, k3, k4, output_fcn, MAX_OUTPUT_LENGTH)
%DIRECTMETHOD Simulate a reaction network using the Gillespie direct method.
%   Usage:
%       [t, x] = directMethod(stoich_matrix, propensity_fcn, tspan, x0, ...
%           k1, k2, k3, k4, output_fcn, MAX_OUTPUT_LENGTH)
%
%   Returns:
%       t              - Times of saved states (N_saved x 1).
%       x              - Species counts (N_saved x N_species).
%
%   Required:
%       tspan          - Initial and final times, [t_init, t_final].
%
%       x0             - Initial species counts, [S1_0, S2_0, ...].
%
%       stoich_matrix  - Stoichiometry matrix (N_reactions x N_species).
%                       Each row gives the stoichiometry of a reaction.
%
%       propensity_fcn - Function a = f(x, k1, k2, k3, k4), returning a
%                       column vector in the order of stoich_matrix rows.
%       k1,...,k4      - Reaction rate parameters passed to propensity_fcn.
%
%   Optional:
%       output_fcn        - Callback status = f(t, x), with x a column vector.
%                           A nonzero status stops simulation; default [].
%       MAX_OUTPUT_LENGTH - Maximum number of saved states; default 1e6.
%
%   Reference:
%       Gillespie, D.T. (1977) Exact Stochastic Simulation of Coupled
%       Chemical Reactions. J Phys Chem, 81:25, 2340-2361.
%
%   Nezar Abdennur, 2012 <nabdennur@gmail.com>
%   Dynamical Systems Biology Laboratory, University of Ottawa
%   www.sysbiolab.uottawa.ca
%   Created: 2012-01-19
if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end
if ~exist('params', 'var')
    params = [];
end

%% Initialize
%num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;
%% Simulate reaction events
while T(rxn_count) < tspan(2)
    % Calculate reaction propensities.
    a = propensity_fcn(X(rxn_count,:),  k1, k2, k3, k4);

    % Sample the waiting time to the next reaction (tau).
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0; %(1/a0)*log(1/r(1));

    % Sample the next reaction channel (mu).
    [~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);

    % ...alternatively...
    %mu = find((cumsum(a) >= r(2)*a0), 1,'first');

    % ...or...
    %mu=1; s=a(1); r0=r(2)*a0;
    %while s < r0
    %   mu = mu + 1;
    %   s = s + a(mu);
    %end
    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count,:);
        warning('SSA:ExceededCapacity',...
            'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
        return;
    end

    % Update time and carry out reaction mu.
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);
    rxn_count = rxn_count + 1;

    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count,:)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:TerminalEvent',...
                'Simulation was terminated by OutputFcn.');
            return;
        end
    end
end
% Return the simulation time course.
t = T(1:rxn_count);
x = X(1:rxn_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(rxn_count-1,:);
end
end
