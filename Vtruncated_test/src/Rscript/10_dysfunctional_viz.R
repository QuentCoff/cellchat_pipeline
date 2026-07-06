#!/usr/bin/env Rscript
# 10_dysfunctional_viz.R
# CellChat Procedure 2 Step 13: Visualize up/down regulated signaling events
# Option A: bubble plot
# Option B: chord diagram
# Option C: wordcloud plot
#
# Usage:
#   Rscript 10_dysfunctional_viz.R <prefix1> <prefix2>
#
# Example:
#   Rscript 10_dysfunctional_viz.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 10_dysfunctional_viz.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list <- file.path(base_dir, "Vtruncated_test", "Result", "data", "merged",
                        paste0(prefix1, "_vs_", prefix2),
                        paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

rdata_merged <- file.path(base_dir, "Vtruncated_test", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

step12_dir <- file.path(base_dir, "Vtruncated_test", "Result", "plot",
                        paste0(prefix1, "_vs_", prefix2), "step12")

out_dir <- file.path(base_dir, "Vtruncated_test", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step13")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}
if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading CellChat objects ===\n")
load(rdata_list)    # loads 'object.list'
load(rdata_merged)  # loads 'cellchat'

dataset_names <- names(object.list)
cat(paste0("Datasets: ", paste(dataset_names, collapse = ", "), "\n"))

# ============================================
# Step 13: Dysfunctional signaling visualization
# ============================================

# Load DEG results from Step 12 if available
net_up_file   <- file.path(step12_dir, "net_up.csv")
net_down_file <- file.path(step12_dir, "net_down.csv")

cat("\n=== Step 13: Dysfunctional signaling visualization ===\n")

# --- Recompute DEG if CSVs missing ---
if (!file.exists(net_up_file) || !file.exists(net_down_file)) {
  cat("  Step 12 results not found — recomputing DEG...\n")
  pos.dataset   <- dataset_names[2]
  features.name <- pos.dataset

  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset   = pos.dataset,
    features.name = features.name,
    only.pos      = FALSE,
    thresh.pc     = 0.1,
    thresh.fc     = 0.05,
    do.fast       = TRUE
  )
  net <- netMappingDEG(cellchat, features.name = features.name)
  net.up   <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,
                                 ligand.logFC = 0.05, receptor.logFC = NULL)
  net.down <- subsetCommunication(cellchat, net = net, datasets = dataset_names[1],
                                 ligand.logFC = -0.05, receptor.logFC = NULL)
} else {
  cat("  Loading Step 12 results...\n")
  net.up   <- read.csv(net_up_file,   stringsAsFactors = FALSE)
  net.down <- read.csv(net_down_file, stringsAsFactors = FALSE)
}

cat(paste0("  Up-regulated interactions: ", nrow(net.up), "\n"))
cat(paste0("  Down-regulated interactions: ", nrow(net.down), "\n\n"))

# ============================================
# Option A: Bubble plot
# ============================================

cat("=== Step 13A: Bubble plots ===\n")

# --- 13A-i/ii: Upregulated ---
tryCatch({
  pairLR.use.up <- net.up[, "interaction_name", drop = FALSE]
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up,
                          sources.use = c("Sertoli", "SSC"),
                          targets.use = c("Sertoli", "SSC"),
                          comparison = c(1, 2), angle.x = 90,
                          remove.isolate = TRUE,
                          title.name = paste0("Up-regulated signaling in ", dataset_names[2]))
  ggplot2::ggsave(
    filename = file.path(out_dir, "bubble_up_dysfunctional.png"),
    plot     = gg1,
    width    = 14,
    height   = 12,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: bubble_up_dysfunctional.png\n")
}, error = function(e) {
  cat(paste0("  FAILED up bubble: ", e$message, "\n"))
})

# --- 13A-iii/iv: Downregulated ---
tryCatch({
  pairLR.use.down <- net.down[, "interaction_name", drop = FALSE]
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                          sources.use = c("Sertoli", "SSC"),
                          targets.use = c("Sertoli", "SSC"),
                          comparison = c(1, 2), angle.x = 90,
                          remove.isolate = TRUE,
                          title.name = paste0("Down-regulated signaling in ", dataset_names[2]))
  ggplot2::ggsave(
    filename = file.path(out_dir, "bubble_down_dysfunctional.png"),
    plot     = gg2,
    width    = 14,
    height   = 12,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: bubble_down_dysfunctional.png\n")
}, error = function(e) {
  cat(paste0("  FAILED down bubble: ", e$message, "\n"))
})

# ============================================
# Option B: Chord diagram
# ============================================

cat("\n=== Step 13B: Chord diagrams ===\n")

# --- 13B-i: Upregulated in dataset 2 ---
tryCatch({
  pdf(file.path(out_dir, "chord_up.pdf"),
      width = 10, height = 10)
  netVisual_chord_gene(object.list[[2]],
    sources.use = c("Sertoli", "SSC"),
    targets.use = c("Sertoli", "SSC"),
    slot.name = "net", net = net.up,
    lab.cex = 0.8, small.gap = 3.5,
    title.name = paste0("Up-regulated signaling in ", dataset_names[2]))
  dev.off()
  cat("  Saved: chord_up.pdf\n")
}, error = function(e) {
  cat(paste0("  FAILED up chord: ", e$message, "\n"))
})

# --- 13B-ii: Downregulated in dataset 2 (visualized in dataset 1) ---
tryCatch({
  pdf(file.path(out_dir, "chord_down.pdf"),
      width = 10, height = 10)
  netVisual_chord_gene(object.list[[1]],
    sources.use = c("Sertoli", "SSC"),
    targets.use = c("Sertoli", "SSC"),
    slot.name = "net", net = net.down,
    lab.cex = 0.8, small.gap = 3.5,
    title.name = paste0("Down-regulated signaling in ", dataset_names[2]))
  dev.off()
  cat("  Saved: chord_down.pdf\n")
}, error = function(e) {
  cat(paste0("  FAILED down chord: ", e$message, "\n"))
})

# ============================================
# Option C: Wordcloud plot
# ============================================

cat("\n=== Step 13C: Wordcloud plots ===\n")

# --- 13C-i: Enriched ligands in dataset 2 (up) ---
tryCatch({
  pdf(file.path(out_dir, "wordcloud_up.pdf"),
      width = 10, height = 8)
  computeEnrichmentScore(net.up, species = "human")
  dev.off()
  cat("  Saved: wordcloud_up.pdf\n")
}, error = function(e) {
  cat(paste0("  FAILED up wordcloud: ", e$message, "\n"))
})

# --- 13C-ii: Enriched ligands in dataset 1 (down) ---
tryCatch({
  pdf(file.path(out_dir, "wordcloud_down.pdf"),
      width = 10, height = 8)
  computeEnrichmentScore(net.down, species = "human")
  dev.off()
  cat("  Saved: wordcloud_down.pdf\n")
}, error = function(e) {
  cat(paste0("  FAILED down wordcloud: ", e$message, "\n"))
})

cat("\n--- Procedure 2 Step 13 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
