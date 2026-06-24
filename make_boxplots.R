# Boxplots of MSE across replications (Figures: Boxplot1 / Boxplot2).
# Run AFTER run_simulation.R has produced results/scenarioX_mse.csv.
#   Rscript make_boxplots.R 1     # or 2
source("R/simulation.R")  # for METHODS ordering
args <- commandArgs(trailingOnly = TRUE)
s    <- if (length(args) >= 1) args[1] else "1"

mse <- read.csv(sprintf("results/scenario%s_mse.csv", s))
mse$Method <- factor(mse$Method, levels = METHODS)
cfg <- unique(mse[, c("n", "p", "rho")])
cfg <- cfg[order(cfg$n, cfg$rho), ]

dir.create("results", showWarnings = FALSE)
png(sprintf("results/boxplot_scenario%s.png", s), width = 1400, height = 1200, res = 130)
op <- par(mfrow = c(3, 3), mar = c(6, 4, 2, 1))
for (i in seq_len(nrow(cfg))) {
  sub <- mse[mse$n == cfg$n[i] & mse$p == cfg$p[i] & mse$rho == cfg$rho[i], ]
  boxplot(MSE ~ Method, data = sub, las = 2, xlab = "", ylab = "MSE",
          main = sprintf("n=%d, p=%d, rho=%.2f", cfg$n[i], cfg$p[i], cfg$rho[i]),
          col = "grey85")
}
par(op); dev.off()
message(sprintf("Wrote results/boxplot_scenario%s.png", s))
