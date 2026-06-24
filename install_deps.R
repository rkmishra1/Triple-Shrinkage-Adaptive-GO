# Install CRAN packages: competing estimators (glmnet, ncvreg) and the two
# real-data sets (pls -> gasoline NIR, lars -> diabetes).
pkgs <- c("glmnet", "ncvreg", "pls", "lars")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
cat("Dependencies ready:", paste(pkgs, collapse = ", "), "\n")
