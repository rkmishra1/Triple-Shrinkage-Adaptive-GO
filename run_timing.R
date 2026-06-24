# ===========================================================================
# Reproduce Table 6.4: computation time (seconds) of the final fit for each
# estimator as p grows, with n = 500 fixed and AR(1) correlation rho = 0.75.
# Coordinate-descent fits use tol = 1e-10, maxit = 5000 (manuscript protocol).
#
#   Rscript run_timing.R
# ===========================================================================

source("R/adaptive_go.R")
source("R/competitors.R")
source("R/simulation.R")

n     <- 500
P     <- c(100, 500, 1000, 3000, 5000)
rho   <- 0.75
TOL   <- 1e-10
MAXIT <- 5000L

time_fit <- function(expr) suppressWarnings(system.time(expr)["elapsed"])

rows <- list()
for (p in P) {
  set.seed(2024)
  dat <- generate_data(n, p, rho)
  st  <- standardize(dat$X, dat$y); Xs <- st$Xs; yc <- st$yc
  w   <- adaptive_weights(init_coef(Xs, yc))
  o   <- rep(1, p)
  l1  <- make_l1seq(Xs, yc, o, nl1 = 100L)          # 100-point lambda grid
  rows[[length(rows) + 1]] <- data.frame(
    p        = p,
    Lasso    = time_fit(fit_glmnet_bic(Xs, yc, alpha = 1)),
    Enet     = time_fit(fit_glmnet_bic(Xs, yc, alpha = 0.5)),
    `Ad-Lasso` = time_fit(fit_glmnet_bic(Xs, yc, alpha = 1,   penalty.factor = w)),
    `Ad-Enet`  = time_fit(fit_glmnet_bic(Xs, yc, alpha = 0.5, penalty.factor = w)),
    SCAD     = time_fit(fit_scad_bic(Xs, yc)),
    GO       = time_fit(ago_path(Xs, yc, o, l1, 0.1, 0.9, TOL, MAXIT)),
    `Ad-GO`  = time_fit(ago_path(Xs, yc, w, l1, 0.1, 0.9, TOL, MAXIT)),
    check.names = FALSE, row.names = NULL)
}
tab <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(tab, "results/timing.csv", row.names = FALSE)
print(tab, digits = 3)
