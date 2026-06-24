# Adaptive GO: Coordinate-Descent Estimation and Simulation Study

Reproducible `R` code for the computational study of the **adaptive Generalized
O'Sullivan (adaptive GO / Ad-GO)** penalized regression estimator. The Ad-GO
estimator is fitted with the **pathwise coordinate-descent algorithm** of
Section 5, and the regularization parameters are selected by **BIC** *and*
**cross-validation** (the two tuning strategies are provided for Ad-GO; the
competing methods are tuned by BIC, following the manuscript).

## The estimator

For a standardized design (columns centred, `(1/n) Σ x_ij² = 1`) the Ad-GO
coordinate update is the soft-thresholding rule (Eq. 5.4):

```
β_j = (1 / (1 + λ₂)) · S( (1 + κ λ₂) · b_ls_j , λ₁ w_j ),   S(z,t) = sign(z)(|z| − t)₊
```

where `b_ls_j` is the univariate LS coefficient of the partial residual on
column `j`, `w_j` are adaptive penalty weights, and `(λ₁, λ₂, κ)` are tuning
parameters. Setting `w_j ≡ 1` recovers the (non-adaptive) **GO** estimator.

Adaptive weights are `w_j = 1/|β̂_j^init|^γ` (γ = 1) with `β̂^init` from OLS
when `p < n` and ridge otherwise (Section 5).

## Methods compared

| Method        | Implementation                | Tuning |
|---------------|-------------------------------|--------|
| Lasso         | `glmnet` (α = 1)              | BIC    |
| ElasticNet    | `glmnet` (α = 0.5, config.)   | BIC    |
| Adaptive Lasso| `glmnet` + penalty weights    | BIC    |
| Adaptive ENet | `glmnet` + penalty weights    | BIC    |
| SCAD          | `ncvreg`                      | BIC    |
| GO            | coordinate descent (`R/adaptive_go.R`) | BIC |
| **Ad-GO**     | coordinate descent (`R/adaptive_go.R`) | **BIC and CV** |

All estimators are fitted on the same standardized design and converted back to
the original predictor scale before computing the MSE.

## Repository layout

```
R/adaptive_go.R    coordinate descent, lambda-path, BIC and CV tuning for Ad-GO/GO
R/competitors.R    glmnet/ncvreg wrappers, BIC path selection, adaptive weights
R/simulation.R     data generation, standardization, (CZ, IZ, MSE), per-config driver
run_simulation.R   Tables 6.1 / 6.2  (Scenario 1 and Scenario 2)
run_timing.R       Table 6.4  (computation times vs p, n = 500)
make_boxplots.R    MSE boxplots (Figures Boxplot1 / Boxplot2)
tests/test_ago.R   self-checks for the coordinate-descent core
install_deps.R     installs glmnet and ncvreg
```

## Data-generating model

`y = Xβ* + ε`, `ε ~ N(0, 6²)`, `X ~ N_p(0, Σ)` with `Σ_{jk} = ρ^{|j−k|}`.
The active set is `A = {1,…,s}`, `s = 3·⌊p/9⌋`; for `j ∈ A`,
`β*_j = ξ_j u_j` with `u_j ~ Unif(1,3)` and `ξ_j` a Rademacher sign.

* **Scenario 1** `p = ⌊4√n⌋ − 5`
* **Scenario 2** `p = ⌊4 n^{2/3}⌋ − 5`

for `n ∈ {200, 400, 800}` and `ρ ∈ {0, 0.5, 0.75}`, over 100 replications.

**Metrics.** `MSE = (β̂ − β*)ᵀ Σ (β̂ − β*)`; `CZ` = correctly identified zeros;
`IZ` = active coefficients wrongly set to zero.

## Running

```bash
Rscript install_deps.R          # one-time: glmnet, ncvreg
Rscript tests/test_ago.R        # sanity-check the estimator (fast)

Rscript run_simulation.R both 100   # full study -> results/scenario{1,2}_summary.csv
Rscript run_simulation.R 1 25       # quick: Scenario 1, 25 reps
Rscript make_boxplots.R 1           # MSE boxplots for Scenario 1
Rscript run_timing.R                # timing table -> results/timing.csv
```

Outputs are written to `results/` (git-ignored): per-configuration summaries
(`CZ`, `IZ`, mean `MSE`, `MSE_SD`), the raw per-replication MSEs for the
boxplots, and the timing table.

## Notes on the BIC

The selection criterion used is the standard regression BIC
`n·log(RSS/n) + log(n)·df` with `df` = number of nonzero coefficients
(Wang et al. 2007). This is the operational form of manuscript Eq. (5.5):
written literally as `log(RSS) + log(n)·df` (no `n` factor) the df penalty
dominates the fit term and selects the null model, so the conventional
`n`-scaled version is used; it gives the same ranking with sensible selection.

## Reproducibility / tuning knobs

`fit_all()` (in `R/simulation.R`) exposes the Ad-GO/GO tuning grids:
`l2seq` (default `c(0, 0.01, 0.1, 1)`), `kapseq` (default `c(0.3, 0.6, 0.9)`),
`nl1` (number of `λ₁` grid points), `nfolds` (CV folds), `enet_alpha`,
`gamma`, plus coordinate-descent `tol`/`maxit`. The full study is compute-heavy
in Scenario 2 (`p` up to 339); reduce `nreps`/`nl1` for a quick pass.
```
