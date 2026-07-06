#!/usr/bin/env Rscript
# 02_compare_interactions.R
# CellChat Procedure 2 Step 5: Identify altered interactions and cell populations
# Compare total number of interactions (Option A) and interaction strength (Option B)
#
# Usage:
#   Rscript 02_compare_interactions.R <prefix1> <prefix2>
#
# Example:
#   Rscript 02_compare_interactions.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 02_compare_interactions.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# Read merged object from Vtruncated_5Clusters/Result/data/merged
rdata_merged <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

# Write plots to Vtruncated_5Clusters/Result/plot
out_dir <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step5")

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
    width = 600, height = 500, res = 120)
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
print(gg1)
dev.off()
cat(paste0("  Saved: compare_interactions_count.png\n"))

# ============================================
# Step 5B: Comparing total interaction strength
# ============================================

cat("\n=== Step 5B: Comparing interaction strength ===\n")

png(file.path(out_dir, "compare_interactions_weight.png"),
    width = 600, height = 500, res = 120)
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
print(gg2)
dev.off()
cat(paste0("  Saved: compare_interactions_weight.png\n\n"))

cat("--- Procedure 2 Step 5 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
