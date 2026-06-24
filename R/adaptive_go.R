# ---------------------------------------------------------------------------
# Adaptive GO (Ad-GO) estimator via pathwise coordinate descent.
#
# Implements Eq. (5.4) of the manuscript:
#     beta_j = (1 / (1 + lambda2)) * S( (1 + kappa*lambda2) * b_ls_j , lambda1 * w_j )
# where b_ls_j is the univariate LS coefficient of the partial residual on x_j
# (columns of X standardised so that (1/n) sum_i x_ij^2 = 1), and
#     S(z, t) = sign(z) * (|z| - t)_+   is the soft-thresholding operator.
#
# Setting w_j = 1 for all j gives the (non-adaptive) GO estimator.
# ---------------------------------------------------------------------------

soft_threshold <- function(z, t) sign(z) * pmax(abs(z) - t, 0)

#' One Ad-GO coordinate-descent fit at a fixed (lambda1, lambda2, kappa).
#'
#' @param Xs  n x p standardised design, (1/n) sum_i x_ij^2 = 1, columns centred.
#' @param yc  centred response vector.
#' @param w   length-p penalty weights (w_j = 1 -> GO).
#' @param beta_init  optional warm start.
#' @return length-p coefficient vector (standardised scale).
ago_cd <- function(Xs, yc, lambda1, lambda2, kappa, w,
                   beta_init = NULL, tol = 1e-7, maxit = 1000L) {
  n <- nrow(Xs); p <- ncol(Xs)
  beta <- if (is.null(beta_init)) numeric(p) else beta_init
  r <- yc - drop(Xs %*% beta)            # full residual
  denom <- 1 + lambda2
  Lprev <- Inf
  for (it in seq_len(maxit)) {
    for (j in seq_len(p)) {
      # univariate LS coef of partial residual on x_j (c_j = 1 by standardisation)
      b_ls <- mean(Xs[, j] * r) + beta[j]
      bj <- soft_threshold((1 + kappa * lambda2) * b_ls, lambda1 * w[j]) / denom
      if (bj != beta[j]) r <- r - Xs[, j] * (bj - beta[j])
      beta[j] <- bj
    }
    L <- mean(r^2)                        # loss = (1/n) ||y - X beta||^2
    if (abs(L - Lprev) <= tol) break
    Lprev <- L
  }
  beta
}

#' Solve the Ad-GO path over a decreasing lambda1 grid (warm-started).
#' lambda2 and kappa fixed. lambda1 MUST be decreasing for warm starts to help.
ago_path <- function(Xs, yc, w, l1seq, lambda2, kappa,
                     tol = 1e-7, maxit = 1000L) {
  p <- ncol(Xs)
  B <- matrix(0, p, length(l1seq))
  beta <- numeric(p)
  for (i in seq_along(l1seq)) {
    beta <- ago_cd(Xs, yc, l1seq[i], lambda2, kappa, w,
                   beta_init = beta, tol = tol, maxit = maxit)
    B[, i] <- beta
  }
  B
}

#' Build a decreasing lambda1 grid from the data (KKT-style lambda_max / w).
make_l1seq <- function(Xs, yc, w, nl1 = 25L, ratio = 1e-3) {
  n <- nrow(Xs)
  z <- abs(drop(crossprod(Xs, yc))) / n
  l1max <- max(z / w) * 1.05             # smallest lambda1 giving all-zero fit
  exp(seq(log(l1max), log(l1max * ratio), length.out = nl1))
}

# --- Tuning strategies for Ad-GO -------------------------------------------

#' Select (lambda1, lambda2, kappa) by BIC. Standard regression BIC
#'   BIC = n*log(RSS/n) + log(n)*df    (operational form of Eq. (5.5);
#' df = number of nonzero coefficients).
ago_bic <- function(Xs, yc, w,
                    l2seq = c(0, 0.01, 0.1, 1),
                    kapseq = c(0.3, 0.6, 0.9),
                    nl1 = 25L, tol = 1e-7, maxit = 1000L) {
  n <- nrow(Xs)
  l1seq <- make_l1seq(Xs, yc, w, nl1)
  best <- list(crit = Inf, beta = numeric(ncol(Xs)))
  for (l2 in l2seq) for (kap in kapseq) {
    B <- ago_path(Xs, yc, w, l1seq, l2, kap, tol, maxit)
    rss <- colSums((yc - Xs %*% B)^2)
    rss[rss <= 0] <- .Machine$double.eps
    df  <- colSums(abs(B) > 0)
    bic <- n * log(rss / n) + log(n) * df
    j <- which.min(bic)
    if (bic[j] < best$crit) best <- list(crit = bic[j], beta = B[, j])
  }
  best$beta
}

#' Select (lambda1, lambda2, kappa) by K-fold cross-validation (prediction MSE).
#' Weights are computed once on the full data (standard adaptive-CV practice).
ago_cv <- function(Xs, yc, w,
                   l2seq = c(0, 0.01, 0.1, 1),
                   kapseq = c(0.3, 0.6, 0.9),
                   nl1 = 25L, nfolds = 5L, tol = 1e-7, maxit = 1000L) {
  n <- nrow(Xs)
  l1seq <- make_l1seq(Xs, yc, w, nl1)
  foldid <- sample(rep_len(seq_len(nfolds), n))
  best <- list(err = Inf, l1 = NA, l2 = NA, kap = NA)
  for (l2 in l2seq) for (kap in kapseq) {
    err <- numeric(length(l1seq))
    for (k in seq_len(nfolds)) {
      tr <- foldid != k; te <- !tr
      B <- ago_path(Xs[tr, , drop = FALSE], yc[tr], w, l1seq, l2, kap, tol, maxit)
      pred <- Xs[te, , drop = FALSE] %*% B
      err <- err + colSums((yc[te] - pred)^2)
    }
    j <- which.min(err)
    if (err[j] < best$err)
      best <- list(err = err[j], l1 = l1seq[j], l2 = l2, kap = kap)
  }
  ago_cd(Xs, yc, best$l1, best$l2, best$kap, w, tol = tol, maxit = maxit)
}
