# ===========================================================================
# Reproduce the computational study (Tables 6.1 & 6.2):
#   Scenario 1: p = floor(4*sqrt(n)) - 5      (p = O(n^{1/2}))
#   Scenario 2: p = floor(4*n^{2/3})  - 5      (p = O(n^{2/3}))
# for n in {200,400,800} and rho in {0, 0.5, 0.75}.
#
# Usage:
#   Rscript run_simulation.R                 # full study, 100 reps (slow)
#   Rscript run_simulation.R 2 0.5 20        # scenario 2, only n>=... etc (see below)
#
# Optional positional args:  <scenario:1|2|both>  <nreps>
#   Rscript run_simulation.R both 100        # default
#   Rscript run_simulation.R 1 25            # quick check on scenario 1
# Results are written to results/scenarioX_summary.csv and ..._mse.csv.
# ===========================================================================

source("R/adaptive_go.R")
source("R/competitors.R")
source("R/simulation.R")

args      <- commandArgs(trailingOnly = TRUE)
scenario  <- if (length(args) >= 1) args[1] else "both"
nreps     <- if (length(args) >= 2) as.integer(args[2]) else 100L

p_scen <- list("1" = function(n) floor(4 * sqrt(n)) - 5,
               "2" = function(n) floor(4 * n^(2/3)) - 5)
N    <- c(200, 400, 800)
RHO  <- c(0, 0.5, 0.75)

dir.create("results", showWarnings = FALSE)

run_scenario <- function(s) {
  pf <- p_scen[[s]]
  summ <- list(); msel <- list()
  for (n in N) for (rho in RHO) {
    p <- pf(n)
    message(sprintf("[Scenario %s] n=%d, p=%d, rho=%.2f, reps=%d", s, n, p, rho, nreps))
    res <- run_config(n, p, rho, nreps = nreps, seed = 2024)
    summ[[length(summ) + 1]] <- res$summary
    m <- res$mse; m$n <- n; m$p <- p; m$rho <- rho
    msel[[length(msel) + 1]] <- m
  }
  summary_df <- do.call(rbind, summ)
  mse_df     <- do.call(rbind, msel)
  write.csv(summary_df, sprintf("results/scenario%s_summary.csv", s), row.names = FALSE)
  write.csv(mse_df,     sprintf("results/scenario%s_mse.csv", s),     row.names = FALSE)
  print(summary_df, digits = 3)
  summary_df
}

scenarios <- if (scenario == "both") c("1", "2") else scenario
for (s in scenarios) run_scenario(s)

message("Done. CSVs written to results/.")
