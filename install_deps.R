# Install CRAN packages: competing estimators (glmnet, ncvreg) and the two
# real-data sets (ppls -> cookie NIR spectroscopy, lars -> diabetes).
pkgs <- c("glmnet", "ncvreg", "ppls", "lars")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
cat("Dependencies ready:", paste(pkgs, collapse = ", "), "\n")
