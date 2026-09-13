%DATA_SIMULATION Generate steady-state counts from the telegraph model.
% Parameter order: [k_off, k_on, k_syn]; degradation rate is fixed at 1.

rng(314)

%% Configure simulation
n_obs = 1000;              % Observations per gene.

n_gene = 1;
threshold = 0.5;           % Minimum sample mean for retaining a gene.

%% Preallocate counts and true parameters
% Rows correspond to genes; columns of data correspond to observations.
data = zeros(n_gene, n_obs);
ks = zeros(n_gene, 3);

%% Simulate until enough genes pass the threshold
g = 1;
n_sim = 1;
while g <= n_gene
    % Draw independent log-uniform rates on [0.1, 10].
    k_on = 10.^(unifrnd(-1,1));
    k_off = 10.^(unifrnd(-1,1));
    k_syn = 10.^(unifrnd(-1,1));
    ks(g,:) = [k_off, k_on, k_syn];

    % Sample the beta-Poisson steady-state distribution.
    for rep = 1:n_obs
        p = betarnd(k_on,k_off);
        data(g,rep) = poissrnd(k_syn.*p);
    end
    if mean(data(g,:)) > threshold
        g = g + 1;
    end
    n_sim = n_sim + 1;
end

%% Save observations and true parameters
save('Data.mat','data','n_obs','ks')  
