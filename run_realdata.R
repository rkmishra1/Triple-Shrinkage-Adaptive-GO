# ===========================================================================
# Real-data application: compare all eight estimators on two benchmark sets
#   - gasoline (NIR spectroscopy, p > n)
#   - diabetes (lars)
# via repeated random train/test splits. Reports out-of-sample MSPE and the
# average number of selected variables.
#
#   Rscript run_realdata.R              # both datasets, 100 splits
#   Rscript run_realdata.R 25           # quick: 25 splits
# Writes results/realdata_summary.csv and results/realdata_raw.csv.
# ===========================================================================

source("R/adaptive_go.R")
source("R/competitors.R")
source("R/simulation.R")
source("R/realdata.R")

args     <- commandArgs(trailingOnly = TRUE)
nsplits  <- if (length(args) >= 1) as.integer(args[1]) else 100L

datasets <- c("gasoline", "diabetes")
summ <- list(); raw <- list()
for (ds in datasets) {
  res <- run_realdata(ds, nsplits = nsplits)
  message(sprintf("[%s]  n=%d, p=%d, splits=%d",
                  res$meta$name, res$meta$n, res$meta$p, nsplits))
  summ[[ds]] <- res$summary
  r <- res$raw; r$Dataset <- res$meta$name; raw[[ds]] <- r
}
summary_df <- do.call(rbind, summ)
raw_df     <- do.call(rbind, raw)

dir.create("results", showWarnings = FALSE)
write.csv(summary_df, "results/realdata_summary.csv", row.names = FALSE)
write.csv(raw_df,     "results/realdata_raw.csv",     row.names = FALSE)
print(summary_df, digits = 4, row.names = FALSE)
message("Done. CSVs written to results/.")
