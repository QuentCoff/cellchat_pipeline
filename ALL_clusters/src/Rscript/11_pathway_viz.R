#!/usr/bin/env Rscript
# 11_pathway_viz.R
# CellChat Procedure 2 Step 14: Visually compare inferred cell-cell communication networks
# Multi-group version: circle plots across all conditions.
#
# Usage:
#   Rscript 11_pathway_viz.R [pathway1,pathway2,...]
#
# Examples:
#   Rscript 11_pathway_viz.R
#   Rscript 11_pathway_viz.R CXCL,BMP,WNT

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

# Default: all pathways from the merged object (will be set after loading)
# User can override with: Rscript 11_pathway_viz.R CXCL,BMP
args <- commandArgs(trailingOnly = TRUE)
user_pathways <- NULL
if (length(args) >= 1) {
  user_pathways <- strsplit(args[1], ",")[[1]]
}

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, STEP14_DIR)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading CellChat objects ===\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("Datasets: ", paste(names(object.list), collapse = ", "), "\n"))

# Keep only colors for the cell types present after subsetting
celltype_colors <- CELLTYPE_COLORS[levels(object.list[[1]]@idents)]
celltype_colors <- celltype_colors[!is.na(celltype_colors)]

# Resolve pathway list: user-specified OR union of all significant pathways
if (!is.null(user_pathways)) {
  pathways.show <- user_pathways
} else {
  pathways.show <- Reduce(union, lapply(object.list, function(x) x@netP$pathways))
}
cat(paste0("Pathways to visualize: ", length(pathways.show), " (",
           paste(head(pathways.show, 5), collapse = ", "),
           ifelse(length(pathways.show) > 5, "...)", ")"), "\n\n"))

# ============================================
# Option A: Circle plots
# ============================================

cat("=== Step 14A: Circle plots ===\n")

for (pw in pathways.show) {
  cat(paste0("  Pathway: ", pw, "\n"))

  # Determine which datasets actually contain the pathway
  idx_valid <- which(sapply(object.list, function(x) pw %in% x@netP$pathways))
  if (length(idx_valid) == 0) {
    cat(paste0("    SKIPPED: pathway '", pw, "' not found in any dataset\n"))
    next
  }
  cat(paste0("    Found in: ", paste(names(object.list)[idx_valid], collapse = ", "), "\n"))

  tryCatch({
    # Compute max weight across datasets that have the pathway
    obj_subset <- object.list[idx_valid]
    weight.max <- getMaxWeight(obj_subset, slot.name = c("netP"), attribute = pw)

    nplots <- length(obj_subset)
    if (nplots == 1) {
      plot_w <- PATHWAY_VIZ_WIDTH_1
      plot_h <- PATHWAY_VIZ_HEIGHT_1
      mfrow_layout <- c(1, 1)
    } else if (nplots == 2) {
      plot_w <- PATHWAY_VIZ_WIDTH_2
      plot_h <- PATHWAY_VIZ_HEIGHT_2
      mfrow_layout <- c(1, 2)
    } else {
      # 3 or more: use 2 rows, one panel width per plot
      ncols <- ceiling(nplots / 2)
      plot_w <- PATHWAY_VIZ_WIDTH_N * ncols
      plot_h <- PATHWAY_VIZ_HEIGHT_N
      mfrow_layout <- c(2, ncols)
    }

    png(file.path(out_dir, paste0("circle_", pw, ".png")),
        width = plot_w, height = plot_h, res = PATHWAY_VIZ_RES)
    par(mfrow = mfrow_layout, xpd = TRUE, oma = c(0, 0, 2, 0), mar = c(2, 2, 2, 2))

    for (i in idx_valid) {
      netVisual_aggregate(
        object.list[[i]],
        signaling        = pw,
        layout           = "circle",
        edge.weight.max  = weight.max[1],
        edge.width.max   = 10,
        vertex.label.cex = 1.6,
        signaling.name   = names(object.list)[i],
        color.use        = celltype_colors
      )

      par(xpd = TRUE)
      usr <- par("usr")
      rect(usr[1] - 100, usr[4], usr[2] + 100, usr[4] + 100, col = "white", border = NA)
      title(main = names(object.list)[i], line = 0.5)
    }

    mtext(paste0(pw, " signaling pathway network"), side = 3, line = 2, outer = TRUE, cex = 1.2, font = 2)

    dev.off()
    cat(paste0("    Saved: circle_", pw, ".png\n"))
  }, error = function(e) {
    cat(paste0("    FAILED: ", e$message, "\n"))
  })
}

cat("\n--- Procedure 2 Step 14 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
