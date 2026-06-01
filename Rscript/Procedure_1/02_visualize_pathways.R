#!/usr/bin/env Rscript
# 02_visualize_pathways.R
# CellChat Procedure 1 Step 16: Visualization of individual signaling pathways
#
# Usage:
#   Rscript 02_visualize_pathways.R <prefix> [receiver_indices]
#
# Example:
#   Rscript 02_visualize_pathways.R healthy
#   Rscript 02_visualize_pathways.R healthy 1,2,3,4

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 02_visualize_pathways.R <prefix> [receiver_indices]\n")
}

prefix <- tolower(args[1])
base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

in_rds  <- file.path(base_dir, "results", "Procedure_1", prefix, paste0("cellchat_", prefix, ".rds"))
out_dir <- file.path(base_dir, "results", "Procedure_1", prefix, "visualization", "pathways")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading CellChat object ===\n")
cellchat <- readRDS(in_rds)

pathways.show.all <- cellchat@netP$pathways
cat(paste0("Significant pathways: ", length(pathways.show.all), "\n"))

# Optional: parse vertex.receiver for hierarchy plot (comma-separated indices)
vertex.receiver <- NULL
if (length(args) >= 2) {
  vertex.receiver <- as.numeric(strsplit(args[2], ",")[[1]])
  cat(paste0("Hierarchy vertex.receiver: ", paste(vertex.receiver, collapse = ", "), "\n"))
} else {
  cat("Note: No receiver indices provided, hierarchy plot will be skipped.\n")
}

# ============================================
# STEP 16: Visualize each signaling pathway
# ============================================

for (pathway in pathways.show.all) {
  
  safe_name <- gsub("[ /\\:]", "_", pathway)
  cat(paste0("Processing: ", pathway, "\n"))
  
  # --- A. Circle plot ---
  png(file.path(out_dir, paste0(safe_name, "_circle.png")),
      width = 800, height = 600, res = 100)
  tryCatch({
    netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message))
  })
  dev.off()
  
  # --- B. Hierarchy plot ---
  if (!is.null(vertex.receiver)) {
    png(file.path(out_dir, paste0(safe_name, "_hierarchy.png")),
        width = 800, height = 600, res = 100)
    tryCatch({
      netVisual_aggregate(cellchat, signaling = pathway, layout = "hierarchy",
                          vertex.receiver = vertex.receiver)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message))
    })
    dev.off()
  }
  
  # --- C. Chord diagram ---
  png(file.path(out_dir, paste0(safe_name, "_chord.png")),
      width = 800, height = 600, res = 100)
  tryCatch({
    par(mfrow = c(1, 1))
    netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message))
  })
  dev.off()
  
  # --- D. Heat map ---
  png(file.path(out_dir, paste0(safe_name, "_heatmap.png")),
      width = 800, height = 600, res = 100)
  tryCatch({
    par(mfrow = c(1, 1))
    netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error:", e$message))
  })
  dev.off()
}

cat(paste0("\n--- Step 16 complete ---\n"))
cat(paste0("Output: ", out_dir, "\n"))
