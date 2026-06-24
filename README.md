<div align="center">

# Adaptive GO

### Coordinate-Descent Estimation, Simulation Study & Real-Data Application

*Reproducible `R` code for the **adaptive Generalized O'Sullivan (Ad-GO)** penalized regression estimator.*

[![R](https://img.shields.io/badge/R-%E2%89%A54.0-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Tuning](https://img.shields.io/badge/tuning-BIC%20%7C%20CV-success)](#-tuning)

</div>

---

The Ad-GO estimator is fitted with a **pathwise coordinate-descent** algorithm and tuned by both **BIC** and **cross-validation**. This repository reproduces the full computational study — the simulation tables, the timing table, the MSE boxplots — and a real-data application on two benchmark datasets.

## 📑 Contents

- [The estimator](#-the-estimator)
- [Methods compared](#-methods-compared)
- [Repository layout](#-repository-layout)
- [Quick start](#-quick-start)
- [Simulation study](#-simulation-study)
- [Real-data application](#-real-data-application)
- [Tuning](#-tuning)
- [Reproducibility & knobs](#-reproducibility--knobs)

## 🧮 The estimator

For a standardized design (columns centred, `(1/n) Σ xᵢⱼ² = 1`) the Ad-GO coordinate update is the soft-thresholding rule **(Eq. 5.4)**:

```
β̂ⱼ = 1/(1+λ₂) · S( (1+κλ₂)·b_lsⱼ , λ₁wⱼ ),     S(z,t) = sign(z)·(|z|−t)₊
```

where `b_lsⱼ` is the univariate LS coefficient of the partial residual on column `j`, the `wⱼ` are adaptive penalty weights, and `(λ₁, λ₂, κ)` are the tuning parameters.

| Symbol | Role |
|:------:|------|
| `λ₁`   | ℓ₁ (sparsity) penalty |
| `λ₂`   | ridge-type shrinkage toward `κ·b_ls` |
| `κ`    | shrinkage target multiplier |
| `wⱼ`   | adaptive weights `1/\|β̂ⱼ^init\|^γ` (γ = 1); `β̂^init` = OLS if `p<n`, ridge otherwise |

> Setting `wⱼ ≡ 1` recovers the non-adaptive **GO** estimator.

## ⚖️ Methods compared

| Method | Implementation | Tuning |
|--------|----------------|:------:|
| Lasso | `glmnet` (α = 1) | BIC |
| ElasticNet | `glmnet` (α = 0.5) | BIC |
| Adaptive Lasso | `glmnet` + penalty weights | BIC |
| Adaptive ElasticNet | `glmnet` + penalty weights | BIC |
| SCAD | `ncvreg` | BIC |
| GO | coordinate descent | BIC |
| **Ad-GO** | **coordinate descent** | **BIC *and* CV** |

All estimators are fitted on the same standardized design and converted back to the original scale before evaluation.

## 🗂️ Repository layout

```
R/
├── adaptive_go.R    coordinate descent, λ-path, BIC & CV tuning for Ad-GO / GO
├── competitors.R    glmnet / ncvreg wrappers, BIC path selection, adaptive weights
├── simulation.R     data generation, standardization, (CZ, IZ, MSE), shared fit core
└── realdata.R       dataset loaders, train/test split evaluation (MSPE, model size)

run_simulation.R     Tables 6.1 / 6.2   (Scenario 1 & 2)
run_realdata.R       real-data application (cookie NIR + diabetes)
run_timing.R         Table 6.4           (computation time vs p)
make_boxplots.R      MSE boxplots        (Figures Boxplot1 / Boxplot2)
tests/test_ago.R     self-checks for the coordinate-descent core
install_deps.R       installs glmnet, ncvreg, pls, lars
```

## 🚀 Quick start

```bash
Rscript install_deps.R          # one-time: glmnet, ncvreg, pls, lars
Rscript tests/test_ago.R        # sanity-check the estimator (fast)
```

| Goal | Command |
|------|---------|
| Full simulation (100 reps) | `Rscript run_simulation.R both 100` |
| Quick simulation (Scenario 1, 25 reps) | `Rscript run_simulation.R 1 25` |
| Real-data application | `Rscript run_realdata.R` |
| MSE boxplots (Scenario 1) | `Rscript make_boxplots.R 1` |
| Timing table | `Rscript run_timing.R` |

All outputs are written to `results/` (git-ignored).

## 🔬 Simulation study

Data are generated from `y = Xβ* + ε`, `ε ~ N(0, 6²)`, `X ~ Nₚ(0, Σ)` with `Σⱼₖ = ρ^|j−k|`. The active set is `A = {1,…,s}`, `s = 3⌊p/9⌋`; for `j ∈ A`, `β*ⱼ = ξⱼuⱼ` with `uⱼ ~ Unif(1,3)` and `ξⱼ` a Rademacher sign.

| Scenario | Dimension growth | `p` |
|:--------:|:----------------:|-----|
| **1** | `O(n^{1/2})` | `⌊4√n⌋ − 5` |
| **2** | `O(n^{2/3})` | `⌊4n^{2/3}⌋ − 5` |

with `n ∈ {200, 400, 800}`, `ρ ∈ {0, 0.5, 0.75}`, over 100 replications.

**Metrics** — `MSE = (β̂−β*)ᵀ Σ (β̂−β*)`; `CZ` = correctly identified zeros; `IZ` = active coefficients wrongly set to zero.

## 📈 Real-data application

The eight estimators are compared on two benchmark datasets via repeated random train/test splits. With no ground-truth coefficients, performance is the out-of-sample **MSPE** (mean squared prediction error) and the average **model size** (number of selected variables). Tuning uses the **training fold only**.

| Dataset | Source (R pkg) | `n` | `p` | Response |
|---------|----------------|:---:|:---:|----------|
| **cookie** (NIR spectroscopy) | `ppls` | 72 | 700 | fat content |
| **diabetes** | `lars` | 442 | 10 | disease progression |

> The cookie dough data carries four constituents (fat, sucrose, dry flour, water); `load_dataset("cookie", response=)` selects which one to model (default `"fat"`).

```bash
Rscript run_realdata.R          # 100 splits, 70/30 train/test  ->  results/realdata_summary.csv
Rscript run_realdata.R 25       # quick: 25 splits
```

> Swap in other datasets by editing `load_dataset()` in `R/realdata.R`.

## 🎯 Tuning

The selection criterion is the standard regression BIC

```
BIC = n·log(RSS/n) + log(n)·df,        df = number of nonzero coefficients
```

(Wang et al. 2007) — the operational form of manuscript **Eq. (5.5)**. Written literally as `log(RSS) + log(n)·df` (without the `n` factor) the df penalty dominates the fit term and the null model is always selected, so the conventional `n`-scaled version is used; it gives the same ranking with sensible selection. Ad-GO additionally supports **K-fold cross-validation** (`ago_cv`) over the same `(λ₁, λ₂, κ)` grid.

## 🔧 Reproducibility & knobs

`fit_methods_std()` / `fit_all()` (in `R/simulation.R`) expose the tuning grids and solver settings:

| Argument | Default | Meaning |
|----------|---------|---------|
| `l2seq` | `c(0, 0.01, 0.1, 1)` | candidate `λ₂` values |
| `kapseq` | `c(0.3, 0.6, 0.9)` | candidate `κ` values |
| `nl1` | `25` | number of `λ₁` grid points |
| `nfolds` | `5` | CV folds for Ad-GO (CV) |
| `enet_alpha` | `0.5` | ElasticNet mixing parameter |
| `gamma` | `1` | adaptive-weight exponent |
| `tol`, `maxit` | `1e-7`, `1000` | coordinate-descent convergence |

Replications run in parallel via `parallel::mclapply` (`run_config(..., ncores=)`, default: all-but-one core; serial on Windows; reproducible L'Ecuyer streams when a seed is given). Scenario 2 (`p` up to 339) is compute-heavy — lower `nreps`/`nl1` for a quick pass.

## 📜 License

[MIT](LICENSE) © 2026 R. K. Mishra
