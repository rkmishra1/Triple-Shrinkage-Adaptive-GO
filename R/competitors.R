# ---------------------------------------------------------------------------
# Competing estimators (Lasso, ElasticNet, adaptive Lasso, adaptive ElasticNet,
# SCAD) and shared helpers. All methods are tuned by BIC over their own
# regularisation path, matching the manuscript's protocol (Wang et al. 2007).
#
# Every estimator is fitted on the SAME standardised design Xs and centred
# response yc, and returns coefficients on the standardised scale; the driver
# converts back to the original predictor scale for MSE.
# ---------------------------------------------------------------------------

suppressMessages({
  library(glmnet)
  library(ncvreg)
})

#' BIC of a coefficient path B (p x nlambda). Standard regression BIC
#'   BIC = n*log(RSS/n) + log(n)*df,
#' the operational form of manuscript Eq. (5.5) (same ranking, with the
#' conventional n scaling so the fit term is not dominated by the df penalty;
#' df = number of nonzero coefficients, as in Wang et al. 2007).
bic_value <- function(B, Xs, yc) {
  n <- nrow(Xs)
  rss <- colSums((yc - Xs %*% B)^2)
  rss[rss <= 0] <- .Machine$double.eps
  df  <- colSums(abs(B) > 0)
  n * log(rss / n) + log(n) * df
}

bic_pick <- function(B, Xs, yc) B[, which.min(bic_value(B, Xs, yc))]

#' Preliminary estimator used to build adaptive weights:
#' OLS when p < n, ridge otherwise (manuscript Section 5).
init_coef <- function(Xs, yc, ridge_lambda = 1e-3) {
  p <- ncol(Xs); n <- nrow(Xs)
  if (p < n) {
    drop(qr.solve(crossprod(Xs), crossprod(Xs, yc)))
  } else {
    drop(solve(crossprod(Xs) + ridge_lambda * diag(p), crossprod(Xs, yc)))
  }
}

#' Adaptive penalty weights w_j = 1 / |b_init_j|^gamma.
adaptive_weights <- function(b_init, gamma = 1, eps = 1e-6) {
  1 / (abs(b_init) + eps)^gamma
}

fit_glmnet_bic <- function(Xs, yc, alpha, penalty.factor = rep(1, ncol(Xs))) {
  fit <- glmnet(Xs, yc, alpha = alpha, standardize = FALSE, intercept = FALSE,
                penalty.factor = penalty.factor)
  bic_pick(as.matrix(fit$beta), Xs, yc)
}

#' glmnet path tuned by K-fold CV (lambda.min). Same model, CV instead of BIC.
fit_glmnet_cv <- function(Xs, yc, alpha, penalty.factor = rep(1, ncol(Xs)),
                          nfolds = 5L) {
  cvf <- cv.glmnet(Xs, yc, alpha = alpha, standardize = FALSE, intercept = FALSE,
                   penalty.factor = penalty.factor, nfolds = nfolds)
  as.numeric(coef(cvf, s = "lambda.min"))[-1]    # drop (zero) intercept
}

fit_scad_bic <- function(Xs, yc) {
  fit <- ncvreg(Xs, yc, penalty = "SCAD")
  B <- as.matrix(fit$beta)[-1, , drop = FALSE]   # drop intercept row
  bic_pick(B, Xs, yc)
}

#' SCAD path tuned by K-fold CV (lambda.min).
fit_scad_cv <- function(Xs, yc, nfolds = 5L) {
  cvf <- cv.ncvreg(Xs, yc, penalty = "SCAD", nfolds = nfolds)
  as.numeric(coef(cvf))[-1]                       # at CV-optimal lambda, drop intercept
}
