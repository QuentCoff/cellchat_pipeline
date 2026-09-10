#!/usr/bin/env Rscript
# 13_export_objects.R
# CellChat Procedure 2 Step 16: Export merged CellChat object and object.list
#
# What it does:
#   Loads the final multi-group object list and the merged CellChat object (preferring
#   the DEG-updated version if available) and copies them into the export folder.
#   This provides a single location for downstream use or sharing.
#
# Inputs:
#   - Multi-group object list RData:
#     Result/data/merged/<MERGED_DIR_NAME>/<MERGE_PREFIX>/cellchat_object.list_<MERGE_PREFIX>.RData
#   - Merged or DEG CellChat RData:
#     Result/data/merged/<MERGED_DIR_NAME>/<MERGE_PREFIX>/cellchat_merged_<MERGE_PREFIX>.RData
#     Result/plot/<MERGED_DIR_NAME>/step12/cellchat_deg.RData (preferred)
#   - Configuration file: config.R
#
# Outputs:
#   - Exported RData files in Result/plot/<MERGED_DIR_NAME>/step16/:
#       cellchat_object.list_<MERGE_PREFIX>.RData
#       cellchat_merged_<MERGE_PREFIX>.RData
#
# Previous step: 01_merge_cellchat.R (or 09_bubble_dysfunctional.R for DEG object)
# Next step: none (final export)
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
# Load and export object.list
# ============================================

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("  Loading object.list...\n")
load(rdata_list)  # loads 'object.list' (one CellChat object per condition)
cat(paste0("    Datasets: ", paste(names(object.list), collapse = ", "), "\n"))

# Copy the object list to the export folder.
cat("  Saving object.list...\n")
save(object.list,
     file = file.path(out_dir, paste0("cellchat_object.list_", merge_prefix, ".RData")))

# ============================================
# Load and export merged cellchat (prefer DEG version)
# ============================================

# Prefer the DEG-updated object from step 12; fall back to the plain merged object.
if (file.exists(rdata_deg)) {
  cat("  Loading cellchat (with DEG results)...\n")
  load(rdata_deg)  # loads 'cellchat'
} else if (file.exists(rdata_merged)) {
  cat("  Loading merged cellchat...\n")
  load(rdata_merged)  # loads 'cellchat'
} else {
  stop(paste0("Merged cellchat not found. Run previous steps first.\n"))
}

# Copy the merged cellchat object to the export folder.
cat("  Saving cellchat...\n")
save(cellchat,
     file = file.path(out_dir, paste0("cellchat_merged_", merge_prefix, ".RData")))

cat("\n--- Procedure 2 Step 16 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
cat(paste0("  ", paste0("cellchat_object.list_", merge_prefix, ".RData"), "\n"))
cat(paste0("  ", paste0("cellchat_merged_", merge_prefix, ".RData"), "\n"))
