# ---------------------------------------------------------------------------
# Data generation, standardisation, metrics, and the per-replication driver
# that fits all estimators and records (CZ, IZ, MSE).
# ---------------------------------------------------------------------------

#' Generate one data set:  y = X beta* + eps,  eps ~ N(0, sigma^2),
#' X ~ N_p(0, Sigma) with Sigma_{jk} = rho^|j-k|.
#' Active set A = {1,...,s}, |A| = 3*floor(p/9); for j in A,
#' beta*_j = xi_j * u_j with u_j ~ Unif(1,3), xi_j Rademacher.
generate_data <- function(n, p, rho, sigma = 6) {
  q <- floor(p / 9); s <- 3 * q
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  R <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% R
  beta <- numeric(p)
  if (s > 0) {
    u  <- runif(s, 1, 3)
    xi <- sample(c(-1, 1), s, replace = TRUE)
    beta[seq_len(s)] <- xi * u
  }
  y <- drop(X %*% beta) + rnorm(n, 0, sigma)
  list(X = X, y = y, beta = beta, Sigma = Sigma, active = seq_len(s))
}

#' Centre y and standardise X columns so that (1/n) sum_i x_ij^2 = 1.
#' Returns the column centres and y-mean too, for out-of-sample prediction.
standardize <- function(X, y) {
  cx <- colMeans(X)
  Xc <- sweep(X, 2, cx, "-")
  sx <- sqrt(colMeans(Xc^2)); sx[sx == 0] <- 1
  list(Xs = sweep(Xc, 2, sx, "/"), yc = y - mean(y),
       scale = sx, center = cx, ybar = mean(y))
}

#' MSE = (bhat - beta*)' Sigma (bhat - beta*).  bhat on ORIGINAL scale.
mse_sigma <- function(bhat, beta, Sigma) {
  d <- bhat - beta
  drop(crossprod(d, Sigma %*% d))
}

#' Variable-selection counts: CZ = correctly identified zeros,
#' IZ = active coefficients incorrectly set to zero.
selection_counts <- function(bhat, beta, tol = 1e-8) {
  est_zero  <- abs(bhat) < tol
  true_zero <- beta == 0
  c(CZ = sum(est_zero & true_zero), IZ = sum(est_zero & !true_zero))
}

# Every estimator is tuned by BOTH BIC and cross-validation.
BASE_METHODS <- c("Lasso", "Enet", "Ad-Lasso", "Ad-Enet", "SCAD", "GO", "Ad-GO")
METHODS <- as.vector(t(outer(BASE_METHODS, c("BIC", "CV"),
                             function(m, c) sprintf("%s (%s)", m, c))))

#' Fit every estimator on a standardised design Xs and centred response yc,
#' under BOTH tuning criteria (BIC and K-fold CV). Returns a named list of
#' coefficient vectors on the STANDARDISED scale (names = METHODS).
#' Shared by the simulation (fit_all) and the real-data analysis.
#' enet_alpha: ElasticNet mixing parameter (configurable; default 0.5).
fit_methods_std <- function(Xs, yc, enet_alpha = 0.5, gamma = 1,
                            l2seq = c(0, 0.01, 0.1, 1), kapseq = c(0.3, 0.6, 0.9),
                            nl1 = 25L, nfolds = 5L, tol = 1e-7, maxit = 1000L) {
  w_ad  <- adaptive_weights(init_coef(Xs, yc), gamma)
  w_one <- rep(1, ncol(Xs))
  list(
    "Lasso (BIC)"    = fit_glmnet_bic(Xs, yc, alpha = 1),
    "Lasso (CV)"     = fit_glmnet_cv (Xs, yc, alpha = 1,          nfolds = nfolds),
    "Enet (BIC)"     = fit_glmnet_bic(Xs, yc, alpha = enet_alpha),
    "Enet (CV)"      = fit_glmnet_cv (Xs, yc, alpha = enet_alpha, nfolds = nfolds),
    "Ad-Lasso (BIC)" = fit_glmnet_bic(Xs, yc, alpha = 1,          penalty.factor = w_ad),
    "Ad-Lasso (CV)"  = fit_glmnet_cv (Xs, yc, alpha = 1,          penalty.factor = w_ad, nfolds = nfolds),
    "Ad-Enet (BIC)"  = fit_glmnet_bic(Xs, yc, alpha = enet_alpha, penalty.factor = w_ad),
    "Ad-Enet (CV)"   = fit_glmnet_cv (Xs, yc, alpha = enet_alpha, penalty.factor = w_ad, nfolds = nfolds),
    "SCAD (BIC)"     = fit_scad_bic(Xs, yc),
    "SCAD (CV)"      = fit_scad_cv (Xs, yc, nfolds = nfolds),
    "GO (BIC)"       = ago_bic(Xs, yc, w_one, l2seq, kapseq, nl1, tol, maxit),
    "GO (CV)"        = ago_cv (Xs, yc, w_one, l2seq, kapseq, nl1, nfolds, tol, maxit),
    "Ad-GO (BIC)"    = ago_bic(Xs, yc, w_ad,  l2seq, kapseq, nl1, tol, maxit),
    "Ad-GO (CV)"     = ago_cv (Xs, yc, w_ad,  l2seq, kapseq, nl1, nfolds, tol, maxit)
  )
}

#' Fit every estimator on one simulated data set; return (CZ, IZ, MSE).
fit_all <- function(dat, ...) {
  st <- standardize(dat$X, dat$y)
  fits <- fit_methods_std(st$Xs, st$yc, ...)
  do.call(rbind, lapply(METHODS, function(m) {
    bhat <- fits[[m]] / st$scale                 # back to original scale
    sel  <- selection_counts(bhat, dat$beta)
    data.frame(Method = m, CZ = sel["CZ"], IZ = sel["IZ"],
               MSE = mse_sigma(bhat, dat$beta, dat$Sigma),
               row.names = NULL)
  }))
}

#' Run nreps replications for one (n, p, rho) configuration and aggregate.
#'
#' @param ncores  number of parallel workers. Default uses all but one core
#'   (1 on Windows, where forking via mclapply is unavailable).
run_config <- function(n, p, rho, nreps = 100L, seed = NULL, verbose = TRUE,
                       ncores = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(ncores)) {
    ncores <- if (.Platform$OS.type == "windows") 1L
              else max(1L, parallel::detectCores() - 1L)
  }
  one_rep <- function(r) fit_all(generate_data(n, p, rho), ...)

  if (ncores > 1L) {
    # Reproducible parallel streams (L'Ecuyer) seeded from the call above.
    if (!is.null(seed)) { RNGkind("L'Ecuyer-CMRG"); set.seed(seed) }
    acc <- parallel::mclapply(seq_len(nreps), one_rep,
                              mc.cores = ncores, mc.set.seed = TRUE)
    bad <- vapply(acc, inherits, logical(1), what = "try-error")
    if (any(bad)) stop("fit_all failed in ", sum(bad), " replication(s): ",
                       conditionMessage(attr(acc[[which(bad)[1]]], "condition")))
    if (verbose) message(sprintf("  n=%d p=%d rho=%.2f : %d reps on %d cores",
                                 n, p, rho, nreps, ncores))
  } else {
    acc <- vector("list", nreps)
    for (r in seq_len(nreps)) {
      acc[[r]] <- one_rep(r)
      if (verbose && r %% 10 == 0)
        message(sprintf("  n=%d p=%d rho=%.2f : rep %d/%d", n, p, rho, r, nreps))
    }
  }
  all <- do.call(rbind, acc)
  agg <- do.call(rbind, lapply(METHODS, function(m) {
    s <- all[all$Method == m, ]
    data.frame(n = n, p = p, rho = rho, Method = m,
               CZ = mean(s$CZ), IZ = mean(s$IZ),
               MSE = mean(s$MSE), MSE_SD = sd(s$MSE))
  }))
  list(summary = agg, mse = all[, c("Method", "MSE")])  # mse kept for boxplots
}
