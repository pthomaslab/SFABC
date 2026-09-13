%ABC_INFERENCE Run ABC-SMC inference for the steady-state telegraph model.
% Load Data.mat from the current folder; data has genes in rows and observations in columns.
% Save posterior particles and weights to ABC_SMC.mat for use with plot_results.m.
% Parameter order: [k_off, k_on, k_syn]; degradation rate is fixed at 1.
% Requires the Statistics and Machine Learning Toolbox.

% Fix the random seed for reproducible inference.
rng(314)

%% Load observed data
load("Data.mat")
[n_gene, n_obs] = size(data);

%% Configure inference
d = 5;                     % Highest raw moment order in the distance function.
tols = [10,5,2];           % Acceptance tolerances for successive SMC populations.
n_samples = 1000;          % Accepted particles per population.
rs = [0.25,0.25];          % Relative perturbation widths; one per SMC transition.
gillespie = true;         % true: Gillespie SSA; false: direct steady-state sampling.
wait_bar = true;           % Display inference progress.

%% Preallocate results
% ks_t: genes-by-populations-by-particles-by-3 array of rate parameters.
% ws_t: genes-by-populations-by-particles array of importance weights.
% n_sim: number of simulated datasets per gene and population.
% run_time: elapsed inference time in seconds for each gene.
n_sim = zeros([n_gene,length(tols)]);
ks_t = zeros([n_gene,length(tols),n_samples,3]);
ws_t = zeros([n_gene,length(tols),n_samples]);
run_time = zeros([n_gene,1]);

%% Infer parameters for each gene
for g = 1:n_gene
    data_g = data(g,:);

    % Run ABC-SMC using the selected simulation method.
    tic
    [n_sim_g,ks_t_g,ws_t_g] = ABC_SMC_tele_steady(data_g, d, tols, n_samples, rs, wait_bar, gillespie);
    run_time(g) = toc;

    % Store all SMC populations for this gene.
    n_sim(g,:) = n_sim_g;
    ks_t(g,:,:,:) = ks_t_g;
    ws_t(g,:,:) = ws_t_g;
end

%% Save inference results
% run_time remains in the workspace; it is not included in this MAT file.
save('ABC_SMC.mat','n_sim','ks_t','ws_t','d','tols','rs')
