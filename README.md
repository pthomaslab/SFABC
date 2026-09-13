This repository contains example MATLAB code for the paper **Simulation-free approximate Bayesian computation for stochastic reaction networks**.

The example performs parameter inference for the **steady-state telegraph model** using approximate Bayesian computation (ABC) and simulation-free approximate Bayesian computation (SFABC), both implemented with sequential Monte Carlo (SMC). The prior, perturbation kernel, tolerance schedule, and related inference settings are the same as those used in **Fig. 2 E-F and Supplementary Fig. S2** of the paper.

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox.
- YALMIP and MOSEK for SFABC inference, configured on the MATLAB path with a working MOSEK license.

## Run the example

Set MATLAB's current folder to this directory, then run:

```matlab
Data_simulation
ABC_inference
SFABC_inference
Plot_results
```

1. `Data_simulation.m` generates 1,000 steady-state observations for one gene and saves the counts and true parameters to `Data.mat`. Genes are retained when their sample mean exceeds 0.5.
2. `ABC_inference.m` loads the data and saves ABC-SMC particles and weights to `ABC_SMC.mat`.
3. `SFABC_inference.m` loads the same data and saves SFABC-SMC particles, weights, volume estimates, and diagnostics to `SFABC_SMC.mat`.
4. `Plot_results.m` loads these files and plots the final SMC population for gene 1. Change `gene` and `ind` in that script to select another gene or iteration.

The plot shows both marginal posteriors on the diagonal, ABC joint posteriors below the diagonal, and SFABC joint posteriors above it. Orange markers indicate the true parameters.

## Model and inference settings

The inferred parameters are ordered as `[k_off, k_on, k_syn]`: gene inactivation, gene activation, and transcription rates. The degradation rate is fixed at 1. `k_tx` in `ABC_SMC.m` denotes the same transcription rate as `k_syn`.

| Setting | Value |
| --- | --- |
| Prior | Independent `log10(k) ~ Uniform(-1, 1)` for each inferred rate |
| Perturbation kernel | `k_new = k_old * (1 + u)`, with independent `u ~ Uniform(-r, r)` for each rate |
| Perturbation widths | `rs = [0.25, 0.25]`; proposals outside the prior support are rejected |
| Tolerance schedule | `tols = [10, 5, 2]` |
| Particles per SMC population | `n_samples = 1000` |
| Highest raw moment order | `d = 5` |
| SFABC state-space truncation | `N = 2 * max(data_g)` for each gene |
| SFABC volume estimation | `Min_c = 100`, `max_N = [1e6, 1e6, 1e6]` |
| Random seed | `rng(314)` in each example script |

ABC uses standardized differences between raw moments of orders 1 to 5. SFABC uses moment constraints of orders 0 to 5 with the truncated steady-state chemical master equation.

By default, ABC uses Gillespie SSA (`gillespie = true`), with a simulation end time of 20. Set `gillespie = false` in `ABC_inference.m` to draw simulated observations directly from the beta-Poisson steady-state distribution.

## Functions

| File | Purpose |
| --- | --- |
| `ABC_SMC.m` | ABC-SMC sampling and its local simulation and distance functions |
| `SFABC_SMC_tele_steady.m` | SFABC-SMC sampling |
| `feas_test_tele_steady.m` | Steady-state CME feasibility test |
| `get_matrices.m` | Telegraph-model CME coefficient matrices |
| `get_counts.m` | Unique observed states and their frequencies |
| `Dirichlet.m` | Dirichlet probability sampling from observed frequencies |
| `V_estimation.m` | Relative acceptance-region volume estimation |
| `plot_gene_posteriors.m` | Posterior comparison plots |

The local Gillespie `directMethod` implementation in `ABC_SMC.m` retains its original attribution to Nezar Abdennur (2012).
