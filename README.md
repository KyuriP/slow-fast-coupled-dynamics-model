# When Context Looks Like Coupling in a Slow–Fast Framework for Psychological Dynamics

Analysis package for revised_manuscript_v5, checked on 20 October 2026. This package contains the current simulation code, aggregate empirical results, simulated trajectories, and figure sources. It is based on repository commit 11e3ebd06e5dcff11c8bedb3e5544b1ee445ae17 with the corrections recorded in AUDIT.txt. It makes no claim about publication or review status.

## Start here

Run commands from this directory in separate, fresh R sessions. The checked environment used R 4.4.2. Exact installed package versions are in session-info.txt. Install the packages loaded at the top of each script from CRAN before running it. The code does not install packages automatically.

The four simulations concern a nine-symptom system. The separate NCT check uses twelve symptoms. Files for older exploratory analyses may remain under R/revision_2026 and res/revision_2026, but only the analyses in the map below contribute to this manuscript. In particular, the observation-window analysis and its figures are not part of the submitted supplement.

## Figures and their sources

Paths in the third column are relative to R/revision_2026 unless stated otherwise. Figure filenames preserve their original names, so their embedded numbers do not always equal the manuscript figure numbers.

| Document figure | Figure file | Analysis and plotting source |
|---|---|---|
| Main 1 | figs/slow-fast-model-illustration_v4.png | Schematic copied from the current Overleaf project |
| Main 2 | figs/revision_2026/Figure2_context_baseline.pdf | 01_sim_context_baseline.R, then figures/fig2_context_baseline.R |
| Main 3 | figs/revision_2026/Figure3_perturbation_recovery.pdf | 02_sim_stress_recovery.R, then figures/fig3_perturbation_recovery.R |
| Main 4 | figs/revision_2026/Figure3_regimes.pdf | 03c_sim_feedback_shock_grid_extended.R and 05_supp_regime_history_dependence.R, then figures/fig3_regimes.R |
| Main 5 | figs/revision_2026/Figure4_network_trio.pdf and Figure4_metric_strip.pdf | 04_sim_network_estimation.R and 04b_sim_network_estimation_replicated.R, then figures/fig4_network_trio.R and figures/fig4_metric_strip.R |
| Main 6 | figs/revision_2026/Figure_helius_mechanisms.pdf | empirical/02_fit_context_models.R and empirical/03_sim_mechanisms.R, then figures/fig_helius_mechanisms.R |
| Main 7 | figs/revision_2026/Figure_helius_discriminating.pdf | empirical/04_discriminating_tests.R, then figures/fig_helius_discriminating.R |
| Supplement S1 | figs/revision_2026/fig_sim1_burnin_trace.pdf | 01_sim_context_baseline.R |
| Supplement S2 | figs/revision_2026/fig_sim2_coupling_sensitivity.pdf | 02b_sim_stress_recovery_coupling_sensitivity.R, then 02c_sim2_coupling_sensitivity_figure.R |
| Supplement S3 | figs/nct_structure_weighted.pdf and nct_strength_weighted.pdf | R/scripts/08_(supp)NCT_analysis.R |
| Supplement S4 | figs/nct_edge_fp_weighted.pdf | R/scripts/08_(supp)NCT_analysis.R |

The main Simulation 3 feedback-on versus feedback-off numerical comparison additionally uses 03_sim_feedback.R and res/revision_2026/sim3. The history check uses steps 1250 through 1500, saved every five steps. Simulation 4's network diagrams show the single-run example, whereas its numerical summaries and metric panels use 30 independent replicates.

## Simulation reproduction

For example, regenerate Simulation 1 and its main figure with

```sh
Rscript R/revision_2026/01_sim_context_baseline.R
Rscript R/revision_2026/figures/fig2_context_baseline.R
```

The Simulation 1 script now explicitly selects Mersenne-Twister, Inversion and Rejection before setting seed 2026. A fresh run reproduces the manuscript means and the supplementary burn-in table. Its corrected cache is included. Run the other simulation and plotting scripts listed above in fresh Rscript processes. Their computational demands differ substantially. The NCT simulation performs 100 replicates for each of 27 conditions with 250 permutations per test. These are the settings reported in the supplement, not a high-precision calibration study.

## HELIUS data and empirical reproduction

Restricted participant data are not included. Public empirical files contain aggregate estimates, covariance matrices of fitted parameters, or simulated summary statistics. The saved fits.rds was inspected and contains fitted parameter objects rather than individual records.

The exact preprocessing scripts and the aggregate inputs from the earlier perspective analysis are bundled under inputs/perspective. An authorized analyst may place helius.sav in inputs/perspective/data, or set SLOWFAST_PERSPECTIVE_DIR to a local perspective project containing its R scripts, data/helius.sav, and results/network_invariance_exact. The raw-data directory is excluded by .gitignore. Access to HELIUS data remains subject to the cohort's data-access requirements.

```sh
Rscript R/revision_2026/empirical/02_fit_context_models.R
Rscript R/revision_2026/empirical/03_sim_mechanisms.R
Rscript R/revision_2026/empirical/04_discriminating_tests.R
Rscript R/revision_2026/figures/fig_helius_mechanisms.R
Rscript R/revision_2026/figures/fig_helius_discriminating.R
```

Step 1 fits nine models both without and with demographic covariates and uses five common cross-validation folds. The unadjusted sample has 23,283 participants. One participant lacks ethnicity, leaving 23,282 for adjusted fits. All 18 saved likelihoods were recomputed against the local HELIUS data during the audit and matched to numerical precision. Cross-validation and the full Monte Carlo analyses were checked against their cached outputs and scripts, but were not rerun in full during this audit.

Step 2 uses the models with linear activation shifts for its context and context-plus-coupling mechanisms. The model with nonlinear activation and continuous coupling moderation is a separate fit comparison in Step 1. It must not be substituted when interpreting the Step 2 figures or the reported linear-model edge effects.

The age and ethnicity checks use the same participants with additional groupings. They are not independent or out-of-sample validation. Their result file is named b_additional_groupings.csv in this package. The LR_J_or_noncentrality column contains observed LR values for observed rows and population contrasts for model rows. Add 36 to the latter only when using the approximate mean finite-sample LR described in the manuscript. The plotted intervals come from replicated samples and require no such correction.

Set SLOWFAST_DATA_MODE=synthetic to exercise the empirical pipeline without restricted data. This mode uses the bundled aggregate inputs. It produces artificial data and must not be interpreted as a replication of the HELIUS findings. Run-size variables are documented at the top of the scripts.

## Scope of the findings

The empirical models evaluate the conditional cross-sectional distribution of symptoms. They do not estimate temporal recovery or symptom-to-context feedback. The age and ethnicity comparisons and the failures to reproduce some observed statistics are part of the reported results. NCT rejection under unequal thresholds also involves a departure from exchangeability, as discussed in Supplement S2.

## Verification and license

AUDIT.txt records corrections and verification limits. figure_manifest.csv maps manuscript assets and records SHA-256 hashes. validation/ contains numerical checks and the reference audit. LICENSE contains the original GPL-3.0 license.
