# ---------------------------------------------------------------------------
# Real-data analysis. The same eight estimators are compared on two benchmark
# data sets by repeated random train/test splits:
#
#   * gasoline  (NIR spectroscopy)  : n = 60,  p = 401  -> octane number
#                                     (R package 'pls'; high-dimensional, p > n)
#   * diabetes  (lars)              : n = 442, p = 10    -> disease progression
#                                     (R package 'lars')
#
# With no ground-truth coefficients, performance is measured out-of-sample:
#   MSPE = mean squared prediction error on the held-out test set,
#   Size = number of selected (nonzero) coefficients.
# Tuning is done on the TRAINING fold only (BIC for all; Ad-GO also via CV).
# ---------------------------------------------------------------------------

#' Load a benchmark data set as list(X, y, name). 'which' in c("gasoline","diabetes").
load_dataset <- function(which = c("gasoline", "diabetes")) {
  which <- match.arg(which)
  if (which == "gasoline") {
    if (!requireNamespace("pls", quietly = TRUE))
      stop("Install the 'pls' package for the gasoline NIR data.")
    data("gasoline", package = "pls", envir = environment())
    g <- get("gasoline", envir = environment())
    list(X = unclass(g$NIR), y = as.numeric(g$octane), name = "gasoline (NIR)")
  } else {
    if (!requireNamespace("lars", quietly = TRUE))
      stop("Install the 'lars' package for the diabetes data.")
    data("diabetes", package = "lars", envir = environment())
    d <- get("diabetes", envir = environment())
    list(X = unclass(d$x), y = as.numeric(d$y), name = "diabetes")
  }
}

#' Fit all methods on a training set and predict a test set.
#' Standardisation and centring use TRAINING statistics only (no leakage).
#' @return data.frame(Method, MSPE, Size).
eval_split <- function(Xtr, ytr, Xte, yte, ...) {
  st  <- standardize(Xtr, ytr)
  fits <- fit_methods_std(st$Xs, st$yc, ...)
  Xte_c <- sweep(Xte, 2, st$center, "-")          # centre test X by train means
  do.call(rbind, lapply(METHODS, function(m) {
    bhat <- fits[[m]] / st$scale                  # original-scale coefficients
    yhat <- st$ybar + drop(Xte_c %*% bhat)        # intercept = train y-mean
    data.frame(Method = m, MSPE = mean((yte - yhat)^2),
               Size = sum(abs(bhat) > 1e-8), row.names = NULL)
  }))
}

#' Repeated random train/test splits on one data set; aggregate MSPE and Size.
run_realdata <- function(which, nsplits = 100L, train_frac = 0.7,
                         seed = 2024, verbose = TRUE, ...) {
  d <- load_dataset(which)
  X <- d$X; y <- d$y; n <- nrow(X)
  set.seed(seed)
  ntr <- floor(train_frac * n)
  acc <- vector("list", nsplits)
  for (b in seq_len(nsplits)) {
    idx <- sample(n, ntr)
    acc[[b]] <- eval_split(X[idx, , drop = FALSE], y[idx],
                           X[-idx, , drop = FALSE], y[-idx], ...)
    if (verbose && b %% 20 == 0)
      message(sprintf("  %s: split %d/%d", d$name, b, nsplits))
  }
  all <- do.call(rbind, acc)
  summ <- do.call(rbind, lapply(METHODS, function(m) {
    s <- all[all$Method == m, ]
    data.frame(Dataset = d$name, Method = m,
               MSPE = mean(s$MSPE), MSPE_SD = sd(s$MSPE),
               Size = mean(s$Size))
  }))
  list(summary = summ, raw = all, meta = list(n = n, p = ncol(X), name = d$name))
}
