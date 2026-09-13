%SFABC_INFERENCE Run SFABC-SMC inference for the steady-state telegraph model.
% Load Data.mat from the current folder; data has genes in rows and observations in columns.
% Save posterior particles and weights to SFABC_SMC.mat for use with plot_results.m.
% Parameter order: [k_off, k_on, k_syn]; degradation rate is fixed at 1.
% Requires the Statistics and Machine Learning Toolbox, YALMIP, and MOSEK.

% Fix the random seed for reproducible inference.
rng(314)

%% Load observed data
load("Data.mat")
[n_gene, n_obs] = size(data);

%% Configure inference
d = 5;                     % Highest raw moment order in the feasibility constraints.
tols = [10,5,2];           % Acceptance tolerances for successive SMC populations.
n_samples = 1000;          % Accepted particles per population.
rs = [0.25,0.25];          % Relative perturbation widths; one per SMC transition.
Min_c = 100;               % Target successes per particle in volume estimation.
max_N = [1e6,1e6,1e6];     % Maximum volume-estimation draws for each population.
wait_bar = true;           % Display inference progress.

%% Preallocate results
% ks_t: genes-by-populations-by-particles-by-3 array of rate parameters.
% ws_t: genes-by-populations-by-particles array of importance weights.
% Vs_t: same shape as ws_t; relative acceptance-region volumes.
% n_opt: number of feasibility optimizations per gene and population.
% run_time: elapsed inference time in seconds for each gene.
n_opt = zeros([n_gene,length(tols)]);
ks_t = zeros([n_gene,length(tols),n_samples,3]);
ws_t = zeros([n_gene,length(tols),n_samples]);
Vs_t = zeros([n_gene,length(tols),n_samples]);
run_time = zeros([n_gene,1]);

%% Infer parameters for each gene
for g = 1:n_gene
    data_g = data(g,:);
    N = max(data_g)*2;     % Truncate mRNA counts at twice the largest observation.

    % Run SFABC-SMC with steady-state CME feasibility constraints.
    tic
    [n_opt_g,G_on,ks_t_g,Vs_t_g,ws_t_g,exact_M,exact_Ps] = SFABC_SMC_tele_steady(data_g, d, tols, N,...
                                       n_samples, rs, Min_c, max_N, wait_bar);
    run_time(g) = toc;

    % Store all SMC populations for this gene.
    n_opt(g,:) = n_opt_g;
    ks_t(g,:,:,:) = ks_t_g;
    ws_t(g,:,:) = ws_t_g;
    Vs_t(g,:,:) = Vs_t_g;
end

%% Save inference results
% Include relative volumes, optimization counts, timings, and sampling settings.
save('SFABC_SMC.mat','n_opt','ks_t','ws_t','Vs_t','run_time','rs','tols','Min_c','max_N')
