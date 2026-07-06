#!/usr/bin/env Rscript
# 03_diff_interactions.R
# CellChat Procedure 2 Step 6:
# Identify substantially altered interactions between cell populations
# by comparing the number of interactions and interaction strength.
# Visualize with circle plot (Option A) or heat map (Option B).
#
# Usage:
#   Rscript 03_diff_interactions.R <prefix1> <prefix2>
#
# Example:
#   Rscript 03_diff_interactions.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
})

cat(paste0("=== CellChat version: ", packageVersion("CellChat"), " ===\n\n"))

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 03_diff_interactions.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# Read merged object from Vtrimean/Result/data/merged
rdata_merged <- file.path(base_dir, "Vtrimean", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

# Write plots to Vtrimean/Result/plot
out_dir <- file.path(base_dir, "Vtrimean", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step6")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading merged CellChat object ===\n")
load(rdata_merged)  # loads 'cellchat'
cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

# ============================================
# Option A: Circle plots (differential)
# ============================================

cat("=== Option A: Circle plots ===\n")

# (i) Examine the differential number of interactions
png(file.path(out_dir, "diff_circle_count.png"), width = 900, height = 800, res = 120)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()
cat("  Saved: diff_circle_count.png\n")

# (ii) Examine the differential interaction strengths
png(file.path(out_dir, "diff_circle_weight.png"), width = 900, height = 800, res = 120)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()
cat("  Saved: diff_circle_weight.png\n\n")

# ============================================
# Option B: Heat maps (differential)
# ============================================

cat("=== Option B: Heat maps ===\n")

# (i) Examine the differential number of interactions
png(file.path(out_dir, "diff_heatmap_count.png"), width = 900, height = 700, res = 120)
gg1 <- netVisual_heatmap(cellchat)
ComplexHeatmap::draw(gg1)
dev.off()
cat("  Saved: diff_heatmap_count.png\n")

# (ii) Examine the differential interaction strengths
png(file.path(out_dir, "diff_heatmap_weight.png"), width = 900, height = 700, res = 120)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
ComplexHeatmap::draw(gg2)
dev.off()
cat("  Saved: diff_heatmap_weight.png\n\n")

cat("--- Procedure 2 Step 6 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
