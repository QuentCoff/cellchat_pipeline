#!/usr/bin/env Rscript
# 03_diff_interactions.R
# CellChat Procedure 2 Step 6:
# Pairwise differential interactions between cell populations.
# Loops through all PAIRWISE comparisons defined in config.R.
#
# Usage:
#   Rscript 03_diff_interactions.R

suppressPackageStartupMessages({
  library(CellChat)
})

cat(paste0("=== CellChat version: ", packageVersion("CellChat"), " ===\n\n"))

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

pair_paths <- function(pair) {
  pair_prefix <- paste(tolower(pair), collapse = "_vs_")
  list(
    prefix = pair_prefix,
    rdata_merged = file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                             MERGED_DIR_NAME, pair_prefix,
                             paste0("cellchat_merged_", pair_prefix, ".RData")),
    out_dir = file.path(base_dir, PROJECT_NAME, "Result", "plot",
                        MERGED_DIR_NAME, pair_prefix, "step6")
  )
}

run_pair <- function(pair) {
  paths <- pair_paths(pair)
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(paths$rdata_merged)) {
    cat(paste0("  SKIPPED: merged RData not found: ", paths$rdata_merged, "\n"))
    return(invisible(NULL))
  }

  cat("=== Loading merged CellChat object ===\n")
  load(paths$rdata_merged)  # loads 'cellchat'
  cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

  # ============================================
  # Option A: Circle plots (differential) -- DISABLED
  # ============================================
  #
  # cat("=== Option A: Circle plots ===\n")
  #
  # png(file.path(paths$out_dir, "diff_circle_count.png"), width = 900, height = 800, res = 120)
  # netVisual_diffInteraction(cellchat, weight.scale = TRUE)
  # dev.off()
  # cat("  Saved: diff_circle_count.png\n")
  #
  # png(file.path(paths$out_dir, "diff_circle_weight.png"), width = 900, height = 800, res = 120)
  # netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
  # dev.off()
  # cat("  Saved: diff_circle_weight.png\n\n")

  # ============================================
  # Option B: Heat maps (differential) -- DISABLED
  # ============================================
  #
  # cat("=== Option B: Heat maps ===\n")
  #
  # png(file.path(paths$out_dir, "diff_heatmap_count.png"), width = 900, height = 700, res = 120)
  # gg1 <- netVisual_heatmap(cellchat)
  # ComplexHeatmap::draw(gg1)
  # dev.off()
  # cat("  Saved: diff_heatmap_count.png\n")
  #
  # png(file.path(paths$out_dir, "diff_heatmap_weight.png"), width = 900, height = 700, res = 120)
  # gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  # ComplexHeatmap::draw(gg2)
  # dev.off()
  # cat("  Saved: diff_heatmap_weight.png\n\n")

  cat("--- Procedure 2 Step 6 complete ---\n")
  cat(paste0("Output: ", paths$out_dir, "\n"))
}

# Run each pairwise comparison -- DISABLED to avoid creating pairwise output directories
# for (pair in PAIRWISE) {
#   cat(paste0("\n##############################################\n"))
#   cat(paste0("# Pairwise comparison: ", paste(pair, collapse = " vs "), "\n"))
#   cat(paste0("##############################################\n"))
#   run_pair(pair)
# }

# ============================================
# Option C: Combined side-by-side diff circle plots
# ============================================

cat("\n=== Option C: Combined side-by-side diff circle plots ===\n")

hc_rdata <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                      MERGED_DIR_NAME, "healthy_vs_crypto",
                      "cellchat_merged_healthy_vs_crypto.RData")
hi_rdata <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                      MERGED_DIR_NAME, "healthy_vs_immuno",
                      "cellchat_merged_healthy_vs_immuno.RData")

if (file.exists(hc_rdata) && file.exists(hi_rdata)) {
  load(hc_rdata)
  cellchat_hc <- cellchat

  load(hi_rdata)
  cellchat_hi <- cellchat

  # Keep only colors for the cell types present after subsetting.
  # Merged CellChat objects do not have levels(cellchat@idents), so we use
  # the rownames of the per-condition count matrix as the cell-type order.
  cell_types <- rownames(cellchat_hc@net[[1]]$count)
  celltype_colors <- CELLTYPE_COLORS[cell_types]
  celltype_colors <- celltype_colors[!is.na(celltype_colors)]

  combined_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                            MERGED_DIR_NAME, "step6")
  dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

  png(file.path(combined_dir, "diff_circle_count_combined.png"),
      width = 1800, height = 800, res = 120)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), xpd = TRUE)
  netVisual_diffInteraction(cellchat_hc, weight.scale = TRUE, title.name = "Healthy vs crypto", color.use = celltype_colors)
  netVisual_diffInteraction(cellchat_hi, weight.scale = TRUE, title.name = "Healthy vs Immuno", color.use = celltype_colors)
  mtext("Differential number of interactions", side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  cat("  Saved: diff_circle_count_combined.png\n")

  png(file.path(combined_dir, "diff_circle_weight_combined.png"),
      width = 1800, height = 800, res = 120)
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), xpd = TRUE)
  netVisual_diffInteraction(cellchat_hc, weight.scale = TRUE, measure = "weight", title.name = "Healthy vs crypto", color.use = celltype_colors)
  netVisual_diffInteraction(cellchat_hi, weight.scale = TRUE, measure = "weight", title.name = "Healthy vs Immuno", color.use = celltype_colors)
  mtext("Differential interaction strength", side = 3, line = 1, outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  cat("  Saved: diff_circle_weight_combined.png\n")
} else {
  cat("  SKIPPED: combined diff circle plots (one or both pairwise merged objects missing)\n")
}
