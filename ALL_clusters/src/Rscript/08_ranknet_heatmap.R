#!/usr/bin/env Rscript
# 08_ranknet_heatmap.R
# CellChat Procedure 2 Step 11: Identify altered signaling with distinct interaction strength
# Multi-group version.
#
# Usage:
#   Rscript 08_ranknet_heatmap.R
#
# Configuration is read from config.R.

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
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

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))

rdata_deg <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                       merged_dir, "step12",
                       "cellchat_deg.RData")

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, "step11")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged, "\nRun 01_merge_cellchat.R first.\n"))
}
if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading CellChat objects ===\n")
load(rdata_merged)  # loads 'cellchat'
load(rdata_list)    # loads 'object.list'
cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

comparison_idx <- seq_along(CONDITIONS)

# ============================================
# Dynamic plot sizing based on pathway / L-R counts
# ============================================
pathway.union <- Reduce(union, lapply(object.list, function(x) x@netP$pathways))
n_pathways <- length(pathway.union)

# L-R pairs count (union across datasets, slot 'net')
lr_all <- lapply(object.list, function(x) {
  if (!is.null(dimnames(x@net$prob))) dimnames(x@net$prob)[[3]] else character(0)
})
n_lr <- length(Reduce(union, lr_all))

cat(paste0("Pathways union: ", n_pathways, " | L-R pairs union: ", n_lr, "\n\n"))

# Width per heatmap
n_datasets <- length(object.list)
width_heatmap <- max(7, 7 * n_datasets * 0.6)

# Helper: ~0.25 in per pathway (min 6, capped at 40)
height_pathways <- min(40, max(6, n_pathways * 0.25))
# L-R pairs are usually 5-10x more numerous
height_lr       <- min(60, max(8, n_lr * 0.18))
# Heatmap row height ~0.3 in per pathway
height_heatmap  <- min(40, max(8, n_pathways * 0.30))

cat(sprintf("Plot heights: pathways=%.1f in | LR=%.1f in | heatmap=%.1f in\n\n",
            height_pathways, height_lr, height_heatmap))

# ============================================
# Step 11A: rankNet overall information flow
# ============================================

# --- 11A-i: Stacked bar chart (pathways, no stat) [DISABLED] ---
# cat("=== Step 11A-i: rankNet pathways stacked (no stat) ===\n")
#
# tryCatch({
#   gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
#                 measure = "weight", sources.use = NULL, targets.use = NULL,
#                 comparison = comparison_idx,
#                 stacked = TRUE, do.stat = FALSE)
#   ggplot2::ggsave(
#     filename = file.path(out_dir, "ranknet_pathways_stacked.png"),
#     plot     = gg,
#     width    = 10,
#     height   = height_pathways,
#     dpi      = 150,
#     limitsize = FALSE
#   )
#   cat("  Saved: ranknet_pathways_stacked.png\n")
# }, error = function(e) {
#   cat(paste0("  FAILED: ", e$message, "\n"))
# })

# --- 11A-ii: Stacked bar chart (pathways, no stat) per pairwise comparison ---
cat("\n=== Step 11A-ii: rankNet pathways stacked (no stat) per pairwise comparison ===\n")

for (pair in PAIRWISE) {
  pair_prefix <- paste(tolower(pair), collapse = "_vs_")
  rdata_pair <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                          MERGED_DIR_NAME, pair_prefix,
                          paste0("cellchat_merged_", pair_prefix, ".RData"))
  if (!file.exists(rdata_pair)) {
    cat(paste0("  SKIPPED: merged RData not found: ", rdata_pair, "\n"))
    next
  }
  load(rdata_pair)  # loads 'cellchat' for this pairwise comparison
  cat(paste0("  Pairwise comparison: ", paste(pair, collapse = " vs "), "\n"))

  # Recompute dynamic height for this pair
  pathway.union.pair <- Reduce(union, lapply(cellchat@netP$similarity[[1]]$dr, function(x) rownames(x)))
  if (length(pathway.union.pair) == 0) {
    pathway.union.pair <- cellchat@netP$pathways
  }
  height_pathways_pair <- min(40, max(6, length(pathway.union.pair) * 0.25))

  tryCatch({
    gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
                  measure = "weight", comparison = c(1, 2),
                  stacked = TRUE, do.stat = FALSE,
                  color.use = CONDITION_COLORS[pair]) +
      ggtitle(paste0("RankNet: ", paste(pair, collapse = " vs "))) +
      theme(
        axis.text.y  = element_text(size = 10),
        legend.text  = element_text(size = 11),
        legend.title = element_text(size = 12)
      )
    filename <- paste0("ranknet_pathways_stacked_", pair_prefix, ".png")
    ggplot2::ggsave(
      filename = file.path(out_dir, filename),
      plot     = gg,
      width    = 10,
      height   = height_pathways_pair,
      dpi      = 150,
      limitsize = FALSE
    )
    cat(paste0("    Saved: ", filename, "\n"))
  }, error = function(e) {
    cat(paste0("    FAILED: ", e$message, "\n"))
  })
}

# Reload multi-group objects for the remaining steps
cat("  Reloading multi-group merged object\n")
load(rdata_merged)  # loads 'cellchat'
load(rdata_list)    # loads 'object.list'

# --- 11A-iii: Stacked bar chart with paired Wilcoxon test (L-R pairs) [DISABLED] ---
# cat("\n=== Step 11A-iii: rankNet L-R pairs stacked (Wilcoxon) ===\n")
#
# tryCatch({
#   gg <- rankNet(cellchat, slot.name = "net", mode = "comparison",
#                 measure = "weight", comparison = comparison_idx,
#                 stacked = TRUE, do.stat = TRUE)
#   ggplot2::ggsave(
#     filename = file.path(out_dir, "ranknet_lr_stacked_stat.png"),
#     plot     = gg,
#     width    = 12,
#     height   = height_lr,
#     dpi      = 150,
#     limitsize = FALSE
#   )
#   cat("  Saved: ranknet_lr_stacked_stat.png\n")
# }, error = function(e) {
#   cat(paste0("  FAILED: ", e$message, "\n"))
# })

# --- 11A-iv: Grouped bar chart (pathways, no stat) ---
cat("\n=== Step 11A-iv: rankNet pathways grouped (no stat) ===\n")

tryCatch({
  gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
                measure = "weight", comparison = comparison_idx,
                stacked = FALSE, do.stat = FALSE,
                color.use = CONDITION_COLORS) +
    ggtitle("Global Signaling Information Flow") +
    theme(
      axis.text.y  = element_text(size = 10),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12)
    )
  ggplot2::ggsave(
    filename = file.path(out_dir, "ranknet_pathways_grouped.png"),
    plot     = gg,
    width    = 10,
    height   = height_pathways,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: ranknet_pathways_grouped.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# ============================================
# Step 11B: ComplexHeatmap outgoing/incoming
# ============================================

cat("\n=== Step 11B: ComplexHeatmap signaling patterns ===\n")

# Pre-requisite: compute network centrality scores for each object
# (required by netAnalysis_signalingRole_heatmap)
cat("  Pre-processing: netAnalysis_computeCentrality on object.list\n")
object.list <- lapply(object.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})
cat("  Done.\n")

has_complexheatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)

if (!has_complexheatmap) {
  cat("  SKIPPED: ComplexHeatmap not installed.\n")
  cat("  Install with: BiocManager::install('ComplexHeatmap')\n")
} else {
  suppressPackageStartupMessages(library(ComplexHeatmap))

  cat(paste0("  Union pathways: ", length(pathway.union), "\n"))

  # Helper to draw heatmaps for all datasets in one PDF
  draw_heatmaps <- function(pattern, filename) {
    cat(paste0("\n=== Step 11B: ", pattern, " signaling heatmaps ===\n"))
    tryCatch({
      hts <- lapply(seq_along(object.list), function(i) {
        netAnalysis_signalingRole_heatmap(object.list[[i]],
          pattern = pattern, signaling = pathway.union,
          title = names(object.list)[i], width = 7, height = height_heatmap)
      })
      ht_list <- hts[[1]]
      if (length(hts) > 1) {
        for (i in 2:length(hts)) {
          ht_list <- ht_list + hts[[i]]
        }
      }
      pdf(file.path(out_dir, filename),
          width = width_heatmap, height = height_heatmap)
      draw(ht_list, ht_gap = unit(0.5, "cm"))
      dev.off()
      cat(paste0("  Saved: ", filename, "\n"))
    }, error = function(e) {
      cat(paste0("  FAILED ", pattern, ": ", e$message, "\n"))
    })
  }

  draw_heatmaps("outgoing", "heatmap_outgoing.pdf")
  draw_heatmaps("incoming", "heatmap_incoming.pdf")
  draw_heatmaps("all", "heatmap_all.pdf")
}

cat("\n--- Procedure 2 Step 11 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
