#!/usr/bin/env Rscript
# 11_pathway_viz.R
# CellChat Procedure 2 Step 14: Visually compare inferred cell-cell communication networks
# Option A: circle plots
#
# Usage:
#   Rscript 11_pathway_viz.R <prefix1> <prefix2> [pathway1,pathway2,...]
#
# Examples:
#   Rscript 11_pathway_viz.R healthy crypto
#   Rscript 11_pathway_viz.R healthy crypto CXCL,BMP,WNT

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 11_pathway_viz.R <prefix1> <prefix2> [pathway1,pathway2,...]\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

# Default: all pathways from the merged object (set after loading)
# User can override with: Rscript 11_pathway_viz.R healthy crypto CXCL,BMP
user_pathways <- NULL
if (length(args) >= 3) {
  user_pathways <- strsplit(args[3], ",")[[1]]
}

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list <- file.path(base_dir, "Vtruncated_5Clusters_test_020", "Result", "data", "merged",
                        paste0(prefix1, "_vs_", prefix2),
                        paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "Vtruncated_5Clusters_test_020", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step14")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading CellChat objects ===\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("Datasets: ", paste(names(object.list), collapse = ", "), "\n"))

# Resolve pathway list: user-specified OR union of all significant pathways
if (!is.null(user_pathways)) {
  pathways.show <- user_pathways
} else {
  pathways.show <- union(object.list[[1]]@netP$pathways,
                         object.list[[2]]@netP$pathways)
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
      plot_w <- 900
      plot_h <- 900
      mfrow_layout <- c(1, 1)
    } else {
      plot_w <- 1800
      plot_h <- 900
      mfrow_layout <- c(1, 2)
    }

    png(file.path(out_dir, paste0("circle_", pw, ".png")),
        width = plot_w, height = plot_h, res = 150)
    par(mfrow = mfrow_layout, xpd = TRUE)

    for (i in idx_valid) {
      netVisual_aggregate(
        object.list[[i]],
        signaling        = pw,
        layout           = "circle",
        edge.weight.max  = weight.max[1],
        edge.width.max   = 10,
        signaling.name   = paste(pw, names(object.list)[i])
      )
    }

    dev.off()
    cat(paste0("    Saved: circle_", pw, ".png\n"))
  }, error = function(e) {
    cat(paste0("    FAILED: ", e$message, "\n"))
  })
}

cat("\n--- Procedure 2 Step 14 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
