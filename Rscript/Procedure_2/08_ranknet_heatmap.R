#!/usr/bin/env Rscript
# 08_ranknet_heatmap.R
# CellChat Procedure 2 Step 11: Identify altered signaling with distinct interaction strength
# Option A: rankNet overall information flow (pathways + L-R pairs)
# Option B: ComplexHeatmap outgoing/incoming signaling patterns
#
# Usage:
#   Rscript 08_ranknet_heatmap.R <prefix1> <prefix2>
#
# Example:
#   Rscript 08_ranknet_heatmap.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 08_ranknet_heatmap.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_merged <- file.path(base_dir, "results", "Procedure_2",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

rdata_list <- file.path(base_dir, "results", "Procedure_2",
                        paste0(prefix1, "_vs_", prefix2),
                        paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "results", "Procedure_2",
                     paste0(prefix1, "_vs_", prefix2), "comparison", "step11")

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

# ============================================
# Step 11A: rankNet overall information flow
# ============================================

# --- 11A-i: Stacked bar chart (pathways, no stat) ---
cat("=== Step 11A-i: rankNet pathways stacked (no stat) ===\n")

tryCatch({
  gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
                measure = "weight", sources.use = NULL, targets.use = NULL,
                stacked = TRUE, do.stat = FALSE)
  ggplot2::ggsave(
    filename = file.path(out_dir, "ranknet_pathways_stacked.png"),
    plot     = gg,
    width    = 10,
    height   = 8,
    dpi      = 150
  )
  cat("  Saved: ranknet_pathways_stacked.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# --- 11A-ii: Stacked bar chart with paired Wilcoxon test (pathways) ---
cat("\n=== Step 11A-ii: rankNet pathways stacked (Wilcoxon) ===\n")

tryCatch({
  gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
                measure = "weight", stacked = TRUE, do.stat = TRUE)
  ggplot2::ggsave(
    filename = file.path(out_dir, "ranknet_pathways_stacked_stat.png"),
    plot     = gg,
    width    = 10,
    height   = 8,
    dpi      = 150
  )
  cat("  Saved: ranknet_pathways_stacked_stat.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# --- 11A-iii: Stacked bar chart with paired Wilcoxon test (L-R pairs) ---
cat("\n=== Step 11A-iii: rankNet L-R pairs stacked (Wilcoxon) ===\n")

tryCatch({
  gg <- rankNet(cellchat, slot.name = "net", mode = "comparison",
                measure = "weight", stacked = TRUE, do.stat = TRUE)
  ggplot2::ggsave(
    filename = file.path(out_dir, "ranknet_lr_stacked_stat.png"),
    plot     = gg,
    width    = 12,
    height   = 10,
    dpi      = 150
  )
  cat("  Saved: ranknet_lr_stacked_stat.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# --- 11A-iv: Grouped bar chart (pathways, no stat) ---
cat("\n=== Step 11A-iv: rankNet pathways grouped (no stat) ===\n")

tryCatch({
  gg <- rankNet(cellchat, slot.name = "netP", mode = "comparison",
                measure = "weight", stacked = FALSE, do.stat = FALSE)
  ggplot2::ggsave(
    filename = file.path(out_dir, "ranknet_pathways_grouped.png"),
    plot     = gg,
    width    = 10,
    height   = 8,
    dpi      = 150
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

  # Union of all pathways across datasets
  pathway.union <- union(object.list[[1]]@netP$pathways,
                         object.list[[2]]@netP$pathways)
  cat(paste0("  Union pathways: ", length(pathway.union), "\n"))

  # --- Outgoing ---
  cat("\n=== Step 11B-i/ii: Outgoing signaling heatmaps ===\n")
  tryCatch({
    ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]],
      pattern = "outgoing", signaling = pathway.union,
      title = names(object.list)[1], width = 7, height = 14)
    ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]],
      pattern = "outgoing", signaling = pathway.union,
      title = names(object.list)[2], width = 7, height = 14)

    pdf(file.path(out_dir, "heatmap_outgoing.pdf"),
        width = 15, height = 14)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    cat("  Saved: heatmap_outgoing.pdf\n")
  }, error = function(e) {
    cat(paste0("  FAILED outgoing: ", e$message, "\n"))
  })

  # --- Incoming ---
  cat("\n=== Step 11B-iii/iv: Incoming signaling heatmaps ===\n")
  tryCatch({
    ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]],
      pattern = "incoming", signaling = pathway.union,
      title = names(object.list)[1], width = 7, height = 14)
    ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]],
      pattern = "incoming", signaling = pathway.union,
      title = names(object.list)[2], width = 7, height = 14)

    pdf(file.path(out_dir, "heatmap_incoming.pdf"),
        width = 15, height = 14)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    cat("  Saved: heatmap_incoming.pdf\n")
  }, error = function(e) {
    cat(paste0("  FAILED incoming: ", e$message, "\n"))
  })

  # --- All (aggregated outgoing + incoming) ---
  cat("\n=== Step 11B-v/vi: All signaling heatmaps ===\n")
  tryCatch({
    ht1 <- netAnalysis_signalingRole_heatmap(object.list[[1]],
      pattern = "all", signaling = pathway.union,
      title = names(object.list)[1], width = 7, height = 14)
    ht2 <- netAnalysis_signalingRole_heatmap(object.list[[2]],
      pattern = "all", signaling = pathway.union,
      title = names(object.list)[2], width = 7, height = 14)

    pdf(file.path(out_dir, "heatmap_all.pdf"),
        width = 15, height = 14)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    dev.off()
    cat("  Saved: heatmap_all.pdf\n")
  }, error = function(e) {
    cat(paste0("  FAILED all: ", e$message, "\n"))
  })
}

cat("\n--- Procedure 2 Step 11 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
