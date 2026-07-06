#!/usr/bin/env Rscript
# 12_gene_expression_viz.R
# CellChat Procedure 2 Step 15: Plot gene expression distribution of signaling genes
# using plotGeneExpression (Seurat wrapper)
#
# Usage:
#   Rscript 12_gene_expression_viz.R <prefix1> <prefix2> [pathway1,pathway2,...]
#
# Examples:
#   Rscript 12_gene_expression_viz.R healthy crypto
#   Rscript 12_gene_expression_viz.R healthy crypto CXCL,BMP,NOTCH

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 12_gene_expression_viz.R <prefix1> <prefix2> [pathway1,pathway2,...]\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

# Default pathways or user-specified
if (length(args) >= 3) {
  pathways.show <- strsplit(args[3], ",")[[1]]
} else {
  pathways.show <- c("PROS")
}

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_merged <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

rdata_deg <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "plot",
                       paste0(prefix1, "_vs_", prefix2), "step12",
                       "cellchat_deg.RData")

out_dir <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step15")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Prefer DEG-updated object if available (has expression data embedded)
if (file.exists(rdata_deg)) {
  cat("=== Loading DEG-updated CellChat object ===\n")
  load(rdata_deg)  # loads 'cellchat'
} else if (file.exists(rdata_merged)) {
  cat("=== Loading merged CellChat object ===\n")
  load(rdata_merged)  # loads 'cellchat'
} else {
  stop(paste0("No CellChat object found. Run previous steps first.\n"))
}

dataset_names <- levels(cellchat@meta$datasets)
cat(paste0("Datasets in object: ", paste(dataset_names, collapse = ", "), "\n"))
cat(paste0("Pathways to visualize: ", paste(pathways.show, collapse = ", "), "\n\n"))

# ============================================
# Step 15: Gene expression distribution
# ============================================

cat("=== Step 15: Gene expression violin plots ===\n")

# Optional: ensure dataset order matches comparison order
# (only if user wants; we keep natural order by default)
# cellchat@meta$datasets <- factor(cellchat@meta$datasets,
#                                   levels = c(prefix1, prefix2))

for (pw in pathways.show) {
  cat(paste0("  Pathway: ", pw, "\n"))

  # Check pathway exists in the merged object DB
  if (!(pw %in% cellchat@DB$interaction$pathway_name)) {
    cat(paste0("    SKIPPED: pathway '", pw, "' not found in CellChat DB\n"))
    next
  }

  # --- Violin plot ---
  tryCatch({
    gg_violin <- plotGeneExpression(
      cellchat,
      signaling   = pw,
      split.by    = "datasets",
      colors.ggplot = TRUE,
      type        = "violin"
    )
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0("violin_", pw, ".png")),
      plot     = gg_violin,
      width    = 14,
      height   = 10,
      dpi      = 150,
      limitsize = FALSE
    )
    cat(paste0("    Saved: violin_", pw, ".png\n"))
  }, error = function(e) {
    cat(paste0("    FAILED violin: ", e$message, "\n"))
  })

  # --- Dot plot (bonus) ---
  tryCatch({
    gg_dot <- plotGeneExpression(
      cellchat,
      signaling   = pw,
      split.by    = "datasets",
      colors.ggplot = TRUE,
      type        = "dot"
    )
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0("dot_", pw, ".png")),
      plot     = gg_dot,
      width    = 14,
      height   = 10,
      dpi      = 150,
      limitsize = FALSE
    )
    cat(paste0("    Saved: dot_", pw, ".png\n"))
  }, error = function(e) {
    cat(paste0("    FAILED dot: ", e$message, "\n"))
  })
}

cat("\n--- Procedure 2 Step 15 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
