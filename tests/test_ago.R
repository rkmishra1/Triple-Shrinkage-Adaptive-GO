# Minimal self-checks for the Ad-GO coordinate-descent core.
#   Rscript tests/test_ago.R
source("R/adaptive_go.R")

# 1) Orthonormal design: CD must equal the closed form Eq. (5.3),
#    beta_j = sign((1+k l2) b_ls) * ((1+k l2)|b_ls| - l1 w_j)_+ / (1+l2).
set.seed(1)
n <- 200; p <- 5
Q  <- qr.Q(qr(matrix(rnorm(n * p), n, p)))   # orthonormal columns
Xs <- Q * sqrt(n)                            # (1/n) sum x_ij^2 = 1
b  <- c(2, -1.5, 0, 0.8, 0)
yc <- drop(Xs %*% b) + rnorm(n, 0, 0.01)
yc <- yc - mean(yc)

l1 <- 0.1; l2 <- 0.05; kap <- 0.3; w <- rep(1, p)
bls <- drop(crossprod(Xs, yc)) / n           # univariate LS = projection
closed <- soft_threshold((1 + kap * l2) * bls, l1 * w) / (1 + l2)
cd <- ago_cd(Xs, yc, l1, l2, kap, w, tol = 1e-12, maxit = 5000)
stopifnot(max(abs(cd - closed)) < 1e-6)

# 2) GO (w=1) recovers the support of a sparse signal in a clean problem.
set.seed(2)
n <- 300; p <- 20
X  <- matrix(rnorm(n * p), n, p)
sx <- sqrt(colMeans(sweep(X, 2, colMeans(X))^2))
Xs <- sweep(sweep(X, 2, colMeans(X)), 2, sx, "/")
beta <- c(rep(3, 4), rep(0, p - 4))
yc <- drop(Xs %*% beta) + rnorm(n, 0, 0.5); yc <- yc - mean(yc)
bh <- ago_bic(Xs, yc, rep(1, p), nl1 = 30)
stopifnot(all(abs(bh[1:4]) > 2))             # strong signals retained near 3
stopifnot(all(abs(bh[5:p]) < 0.2))           # noise coefficients near zero
stopifnot(mean(abs(bh[5:p]) < 1e-8) > 0.5)   # majority of noise exactly zero

# 3) lambda1 large enough -> all coefficients shrunk to zero.
big <- ago_cd(Xs, yc, lambda1 = 1e6, lambda2 = 0.1, kappa = 0.3, w = rep(1, p))
stopifnot(all(big == 0))

cat("All Ad-GO self-checks passed.\n")
