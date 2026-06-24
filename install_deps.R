# Install the two CRAN packages used for the competing estimators.
pkgs <- c("glmnet", "ncvreg")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
cat("Dependencies ready:", paste(pkgs, collapse = ", "), "\n")
