#!/usr/bin/env Rscript
# 13_export_objects.R
# CellChat Procedure 2 Step 16: Export merged CellChat object and object.list
# Multi-group version.
#
# Usage:
#   Rscript 13_export_objects.R
#
# Configuration is read from config.R.

options(stringsAsFactors = FALSE)

# ============================================
# CONFIGURATION
# ============================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("--file=", args)]
if (length(file_arg) == 1) {
  script_dir <- dirname(normalizePath(sub("--file=", "", file_arg)))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "config.R"))

base_dir <- BASE_DIR
merge_prefix <- MERGE_PREFIX
merged_dir <- MERGED_DIR_NAME

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))

rdata_merged <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                          merged_dir, merge_prefix,
                          paste0("cellchat_merged_", merge_prefix, ".RData"))

rdata_deg <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                       merged_dir, STEP12_DIR,
                       "cellchat_deg.RData")

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, STEP16_DIR)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 16: Export CellChat objects ===\n")

# ============================================
# Load object.list
# ============================================

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("  Loading object.list...\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("    Datasets: ", paste(names(object.list), collapse = ", "), "\n"))

# Save to step16 directory
cat("  Saving object.list...\n")
save(object.list,
     file = file.path(out_dir, paste0("cellchat_object.list_", merge_prefix, ".RData")))

# ============================================
# Load merged cellchat (prefer DEG version)
# ============================================

if (file.exists(rdata_deg)) {
  cat("  Loading cellchat (with DEG results)...\n")
  load(rdata_deg)  # loads 'cellchat'
} else if (file.exists(rdata_merged)) {
  cat("  Loading merged cellchat...\n")
  load(rdata_merged)  # loads 'cellchat'
} else {
  stop(paste0("Merged cellchat not found. Run previous steps first.\n"))
}

cat("  Saving cellchat...\n")
save(cellchat,
     file = file.path(out_dir, paste0("cellchat_merged_", merge_prefix, ".RData")))

cat("\n--- Procedure 2 Step 16 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
cat(paste0("  ", paste0("cellchat_object.list_", merge_prefix, ".RData"), "\n"))
cat(paste0("  ", paste0("cellchat_merged_", merge_prefix, ".RData"), "\n"))
