#!/usr/bin/env Rscript
# Install CellChat from GitHub into the current R environment.
# Run after activating the conda/micromamba environment:
#   micromamba activate env_cellchat
#   Rscript install_cellchat.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

# remotes should already be present via environment.yml, but install if missing.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# NMF from conda-forge is too old (0.21.0); CellChat requires >= 0.23.0.
# Install it from CRAN first.
install.packages("NMF")

# Install CellChat v2.1.2 from the official repository.
# This matches the version currently used in the working env_cellchat environment.
remotes::install_github("jinworks/CellChat@v2.1.2", upgrade = "never")

# presto is required by CellChat's identifyOverExpressedGenes for the fast Wilcoxon test.
remotes::install_github("immunogenomics/presto")

cat("\n=== CellChat installed ===\n")
cat(paste0("Version: ", packageVersion("CellChat"), "\n"))
