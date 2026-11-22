# Green Simulation: Reusing the Output of Repeated Experiments — code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24431305.svg)](https://doi.org/10.5281/zenodo.24431305)

Code that reproduces the numerical examples (Section 6: Figures 2–5) of

> Feng, M., & Staum, J. (2017). Green simulation: Reusing the output of repeated experiments.
> *ACM Transactions on Modeling and Computer Simulation*, 27(4), Article 23, 1–28.
> https://doi.org/10.1145/3129130

If you use this code, please cite the paper. The conference version is Feng and Staum (2015),
Green simulation designs for repeated experiments, *Proceedings of the 2015 Winter Simulation
Conference*, 403–413 (code: github.com/mbfeng/GreenSim-WSC2015).

## Quick start

```matlab
% MATLAB R2023b or later with the Statistics and Machine Learning and the
% Optimization Toolboxes (Parallel Computing Toolbox optional, for reruns).
% From this folder:
reproduce_paper
```

By default (`RERUN_SIMULATIONS = false`) this draws Figures 2–5 (`figures/paramhist.pdf`,
`sig2_lomihi.pdf`, `sig2n.pdf`, `sig2n_example2.pdf`, the file names of the paper's LaTeX source)
from the saved summaries of the runs behind the published figures (`results/published_run/`) and
checks every plotted curve against the paper, in a few seconds. Set `RERUN_SIMULATIONS = true` to
rerun the 2 × 10,000 macro replications first (see [Run times](#run-times)).

## The experiments

A sequence of experiments n = 1, 2, … is run with the same simulation model at changing inputs
(*states*) x_n. Each experiment has r replications. Green simulation reuses the outputs of all
earlier experiments to estimate μ(x) = E[F(Y) | x]:

| Estimator | Uses | Paper |
|---|---|---|
| SMC, standard Monte Carlo | the r outputs of the current experiment | Eq. 2 |
| ILR, individual likelihood ratio | all n r outputs, each experiment's outputs reweighted by h(y; x)/h(y; x_i), experiments averaged equally | Eq. 6, Alg. 1 |
| EWLR, weighted ILR | as ILR, experiments weighted by the inverse sample variance of their likelihood-ratio terms | Eq. 9, Alg. 1 |
| MLR, mixture likelihood ratio | all n r outputs, likelihood ratio h(y; x)/((1/n) Σ_i h(y; x_i)) | Eqs. 14–15, Alg. 2 |
| GSK, green stochastic kriging | stochastic kriging fitted to the SMC estimates of all experiments so far (from n = 5) | Section 3, Eq. 19 |

All estimators use the same simulation output (common random numbers). The figure of merit is the
MSE over 10,000 independent macro replications of the whole sequence.

| | CAT bond (Section 6.1) | Credit risk (Section 6.2) |
|---|---|---|
| Output F(Y) | payoff per dollar 1{Y ≤ K} + p 1{Y > K}, Y compound Poisson–exponential loss | 1{L₁ + L₂ > κ}, portfolio loss of two loans (Merton model) |
| State x_n | AR(1) in ℝ², parameters (λ, θ) = ψ(x) by a sigmoid | asset values of two companies, weekly GBM with t-copula innovations |
| Settings | K = 25, p = 0.5, φ = (0.6, 0.5), σ_ε = (0.8, 0.5), λ ∈ (2, 4), θ ∈ (4, 6) | S₀ = (100, 90), D = (85, 85), a = (5, 4), η = (15%, 10%), ς = (30%, 20%), r = 5%, t₃ copula with correlation 0.5, t = 0.5, T = 5, κ = 6 |
| Experiments × replications | 100 × 100 | 52 × 1,000 |
| Targets | current state x_n; fixed states x_mi (λ, θ) = (3, 5) and x_hi = (4, 6) (and x_lo = (2, 4)) | current state |
| True value | Poisson–Gamma series (41 terms) | closed form (t-copula probabilities) |
| Macro replications | 10,000 | 10,000 |

## Paper result → code

| Result | Produced by | Output |
|---|---|---|
| Figure 1 (setting of repeated experiments) | drawn in TikZ in the LaTeX source | — |
| Figure 2 (histograms of λ and θ under the stationary distribution) | `scripts/make_figures.m` (10^6 samples, seed 2) | `figures/paramhist.pdf` |
| Figure 3 (MSEs at the fixed states x_mi, x_hi; CAT bond) | `scripts/run_simulation.m('catbond')` → `scripts/make_figures.m` | `figures/sig2_lomihi.pdf` |
| Figure 4 (MSEs at the current state; CAT bond) | same | `figures/sig2n.pdf` |
| Figure 5 (MSEs at the current state; credit risk) | `scripts/run_simulation.m('credit')` → `scripts/make_figures.m` | `figures/sig2n_example2.pdf` |

## Folder structure

```
reproduce_paper.m            entry point: Figures 2-5 and the check against the paper
functions/
  catbond_settings.m, catbond_states.m, catbond_parameters.m, catbond_simulate.m,
  catbond_payoff.m, catbond_density.m (Eq. 27), catbond_price.m (Eq. 26)
  credit_settings.m, credit_states.m, credit_scenarios.m, credit_density.m,
  credit_loss.m, credit_true_probability.m, scenario_marginals.m, t_copula_uniforms.m
  lr_estimators.m            ILR, EWLR, MLR (Eqs. 6, 9, 15) for one target
  sk_fit.m, sk_predict.m     stochastic kriging (linear trend, Gaussian correlation)
  catbond_one_path.m, credit_one_path.m   one macro replication of each example
scripts/
  run_simulation.m           macro replications (parallel, resumable, seeded) -> summaries
  make_figures.m             Figures 2-5
results/
  published_run/             summaries of the runs behind the published figures
  rerun_seed2017/, rerun_seed2018/   full reruns with this code
verify/
  verify_against_paper.m     checks the curves of Figures 3-5
  paper_*_points.csv         curves plotted in the paper (read from its vector figures)
tests/                       run_tests.m: estimators vs direct formulas, kriging, snapshot
tools/
  extract_published_runs.m   provenance: saved 2016 workspaces -> results/published_run
  digitize_published_figures.py   provenance: paper figures -> verify/paper_*_points.csv
original_2016/               the original 2016 code (reference; see its README)
```

## Notation: code → paper

| Code | Paper | Meaning |
|---|---|---|
| `Parameters(n, :)` | (λ_n, θ_n) = ψ(x_n) | CAT bond state as loss-model parameters |
| `States(n, :)` | x_n = S_0 | credit risk state (asset values) |
| `Losses`, `Scenarios` | Y | random input of the simulation |
| `Payoffs`, `HighLoss` | F(Y) | simulation output |
| `catbond_density`, `credit_density` | h(y; x) | density of Y given the state |
| `MixtureDensity` | h*_n | (1/n) Σ_i h(·; x_i) |
| `NumReplications` | r | replications per experiment |
| `Model.FixedStates.mi`, `.hi` | x_mi, x_hi | fixed states (as parameters) |
| `Model.Evaluation` | T − t | time from the risk horizon to maturity |
| `Fit.Tau2`, `Fit.Theta`, `Fit.Beta` | τ², θ, β | stochastic kriging parameters |

## Verification

`verify/verify_against_paper.m`, the last step of `reproduce_paper.m`:

- **Exact mode** (the published runs): the curves plotted in the paper's vector figures
  (`verify/paper_*_points.csv`, 19 curves, 1,820 points) are compared with the summaries in
  `results/published_run/`. The axis tick labels of these figures are outlines, so the log-log
  axis scales of each panel are fitted to its curves jointly (two parameters per axis); every
  plotted point then lies within 0.0006 px of the summary (the coordinates are printed to three
  decimals).
- **Monte Carlo mode** (full reruns with this code and new random numbers): MSE_rerun and
  MSE_published are compared per period, target and estimator by
  z = (x − v)/√(se_x² + se_v²); |z| > 2.5 is flagged.
  Two full reruns (10,000 macro replications each, as in the paper; `results/rerun_seed2017`,
  `results/rerun_seed2018`), 2,192 comparisons each, about 27 flags expected by chance:

  | Rerun | Flagged | Largest \|z\| | Where |
  |---|---|---|---|
  | seed 2017 | 38 | 3.44 | 21 at x_mi (ILR, EWLR), in neighbouring periods of one curve; the rest scattered |
  | seed 2018 | 10 | 2.92 | scattered (CAT bond SMC, MLR, GSK at x_mi) |

  The x_mi ILR/EWLR flags of seed 2017 do not recur with seed 2018. MSEs of neighbouring periods
  are strongly correlated (the same macro replications), so flags come in runs. The GSK curves at
  the fixed states x_mi and x_hi are 2–4% above the published ones in both reruns (late periods;
  z between 1 and 2.7); GSK at the current state and at x_lo agrees within 1%. This
  stochastic kriging fit is an independent implementation (see below), and the likelihood
  optimizer can stop at different local optima than the 2016 code.

`tests/run_tests.m`: the estimators equal a direct implementation of Equations 2, 6, 9 and 15
(all outputs stored, every likelihood ratio computed) to 1e-15; stochastic kriging interpolates
noiseless data, reproduces a linear surface and is stable under warm starts; `credit_density.m`
equals `copulapdf` times the marginals and stays finite in the far tails; all outputs equal a
saved snapshot.

## Numerical details

- **Common random numbers**: all estimators of a macro replication use the same outputs. Losses
  are drawn by inversion (`poissinv`, `gaminv`) from the macro replication's own stream.
- **Density of the CAT bond loss**: a point mass e^{−λ} at 0 and the continuous part of
  Equation 27; likelihood ratios use this density with respect to (point mass at 0) +
  (Lebesgue measure).
- **EWLR**: an experiment whose likelihood-ratio terms all equal (zero sample variance, e.g., no
  high losses in the credit example) gets weight 0, as in the code that produced the credit risk
  results.
- **Credit risk**: only outputs with F = 1 are stored; the per-experiment sums and sums of
  squares give the same ILR, EWLR and MLR as all outputs. The first state of each path is
  x_1 = S_0 = (100, 90), as in the runs behind the figures.
- **t-copula density** (`credit_density.m`): computed as g(q)/(t(q_1) t(q_2)) with the t
  quantiles q_i from `betaincinv` on the smaller tail. It equals `copulapdf('t', ...)` to 1e-12
  in the bulk; for stored scenarios that are extreme under another period's state (uniforms below
  about 1e-14), `tinv` loses accuracy and returns −Inf near 1e-17, which made `copulapdf`
  return NaN and the MLR estimate undefined in a few macro replications of new seeds.
- **Stochastic kriging** (`sk_fit.m`): linear trend b(x) = [1, x]; design rescaled to [0, 1]²;
  Gaussian correlation exp(−Σ_p θ_p Δu_p²); Σ = τ² R + diag(V), V = sample variance / r; (τ², θ)
  by maximum profile likelihood (`fmincon`, lower bounds 0.001 τ²₀ and 0.001); cold start at
  n = 5, then warm-started from the previous period's fit. This is an independent implementation
  of the model used for the published GSK curves (the stochastic kriging package used in 2016 is
  distributed for personal use only and is not included).
- **Fixed-state SMC**: its variance (1 − p)² P(1 − P)/r is exact, so the SMC lines in Figure 3
  are computed, not simulated.

## Reproducibility

Macro replication m uses substream m of an `mrg32k3a` stream with the run's seed, so serial and
parallel runs give identical results; each run saves its parts and can be resumed.

## Provenance

The published figures were drawn from two 2016 runs whose saved workspaces were reduced to
`results/published_run/` by `tools/extract_published_runs.m`: the CAT bond run
`results_10KReps_with_Trend_EB.mat` (its GSK results come from the run with the Gaussian
correlation function; its SMC curve at the current state was recomputed in a separate pass) and
the credit risk run `results_with_GIBF.mat`. The summaries reproduce the plotted curves exactly.
The script versions that wrote these workspaces were not kept; `original_2016/` holds the latest
versions of the scripts (see its README). This code was rebuilt from the paper and the original
code, tested against direct implementations of the formulas, and verified by full reruns.

## Run times

`reproduce_paper` (default): a few seconds. Full reruns on 44 workers (4 × Intel Xeon Gold 6126,
MATLAB R2025a): CAT bond (with GSK) 8.2–8.5 minutes per seed; credit risk 9.0–10.0 minutes per seed. On one core of an
Intel Core i7-1280P (MATLAB R2023b): about 10 s per CAT bond macro replication with GSK (1.4 s
without), 4.5 s per credit risk macro replication.

## Change log

- **2026-11 (v1.0.0)**: archived version. Rebuilt the 2016 code: one entry point for all
  figures, functions in the paper's notation, both examples in one seeded and resumable driver,
  common random numbers for all estimators, an independent stochastic kriging implementation,
  summaries of the published runs, verification against the published figures, tests.
- **2016**: original code (`original_2016/`).

## License

MIT (see `LICENSE`).
