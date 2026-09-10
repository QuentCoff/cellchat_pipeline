#!/usr/bin/env Rscript
# 02_compare_interactions.R
# CellChat Procedure 2 Step 5: Compare total number of interactions and strength
# Multi-group version: one bar per condition.
#
# Usage:
#   Rscript 02_compare_interactions.R
#
# Configuration is read from config.R.

suppressPackageStartupMessages({
  library(CellChat)
})

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

rdata_merged <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                          merged_dir, merge_prefix,
                          paste0("cellchat_merged_", merge_prefix, ".RData"))

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, STEP5_DIR)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading merged CellChat object ===\n")
load(rdata_merged)  # loads 'cellchat' object
cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

# ============================================
# Step 5A: Comparing total number of interactions
# ============================================

cat("=== Step 5A: Comparing number of interactions ===\n")

png(file.path(out_dir, "compare_interactions_count.png"),
    width = COMPARE_WIDTH, height = COMPARE_HEIGHT, res = COMPARE_RES)
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = seq_along(CONDITIONS), color.use = unname(CONDITION_COLORS[CONDITIONS])) +
  ggtitle(paste0("Number of Interactions\n", paste(CONDITIONS, collapse = " vs "))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
print(gg1)
dev.off()
cat(paste0("  Saved: compare_interactions_count.png\n"))

# ============================================
# Step 5B: Comparing total interaction strength
# ============================================

cat("\n=== Step 5B: Comparing interaction strength ===\n")

png(file.path(out_dir, "compare_interactions_weight.png"),
    width = COMPARE_WIDTH, height = COMPARE_HEIGHT, res = COMPARE_RES)
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = seq_along(CONDITIONS), measure = "weight", color.use = unname(CONDITION_COLORS[CONDITIONS])) +
  ggtitle(paste0("Interaction Strength\n", paste(CONDITIONS, collapse = " vs "))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
print(gg2)
dev.off()
cat(paste0("  Saved: compare_interactions_weight.png\n\n"))

cat("--- Procedure 2 Step 5 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
