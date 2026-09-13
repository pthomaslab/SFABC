%PLOT_RESULTS Compare saved ABC-SMC and SFABC-SMC posterior distributions.
% Parameter order: [k_off, k_on, k_syn].

rng(314)

%% Load saved results and true parameters from this folder
results_dir = fileparts(mfilename('fullpath'));
ABC = load(fullfile(results_dir, 'ABC_SMC.mat'), 'ks_t', 'ws_t', 'tols');
SFABC = load(fullfile(results_dir, 'SFABC_SMC.mat'), 'ks_t', 'ws_t', 'tols');
D = load(fullfile(results_dir, 'Data.mat'), 'ks');

%% Select the gene and SMC iteration
gene = 1;
ind = numel(SFABC.tols);  % Final SMC iteration.

% Compare populations at the same tolerance.
assert(ind <= numel(ABC.tols) && ABC.tols(ind) == SFABC.tols(ind), ...
    'ABC and SFABC must have the same tolerance at the selected iteration.');

%% Plot posterior distributions in log10 space
% Diagonal: both marginals. Lower triangle: ABC. Upper triangle: SFABC.
% Resample 1e5 particles per method; orange markers show the true parameters.
[fig, t, axs] = plot_gene_posteriors( ...
    gene, ind, ABC, SFABC, D.ks, ...
    'TrueColor', [255 164 0]/255, ...
    'N', 1e5, ...
    'NumBins', 25, ...
    'NumBoxes', 25, ...
    'Limits', [-1 1]);
