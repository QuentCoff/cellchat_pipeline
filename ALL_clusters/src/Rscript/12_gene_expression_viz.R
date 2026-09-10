#!/usr/bin/env Rscript
# 12_gene_expression_viz.R
# CellChat Procedure 2 Step 15: Plot gene expression distribution of signaling genes
# Multi-group version: split by dataset.
#
# Usage:
#   Rscript 12_gene_expression_viz.R [pathway1,pathway2,...]
#
# Examples:
#   Rscript 12_gene_expression_viz.R
#   Rscript 12_gene_expression_viz.R CXCL,BMP,NOTCH

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(patchwork)
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

# Default: all pathways from the merged object (set after loading)
# User can override with: Rscript 12_gene_expression_viz.R CXCL,BMP
args <- commandArgs(trailingOnly = TRUE)
user_pathways <- NULL
if (length(args) >= 1) {
  user_pathways <- strsplit(args[1], ",")[[1]]
}

rdata_merged <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                          merged_dir, merge_prefix,
                          paste0("cellchat_merged_", merge_prefix, ".RData"))

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))

rdata_deg <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                       merged_dir, STEP12_DIR,
                       "cellchat_deg.RData")

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, STEP15_DIR)

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

# Load object.list for pathway resolution (DEG object may have empty netP$pathways)
if (file.exists(rdata_list)) {
  load(rdata_list)  # loads 'object.list'
}

dataset_names <- levels(cellchat@meta$datasets)
cat(paste0("Datasets in object: ", paste(dataset_names, collapse = ", "), "\n"))

# Resolve pathway list: user-specified OR union of all datasets' pathways
if (!is.null(user_pathways)) {
  pathways.show <- user_pathways
} else if (exists("object.list") && length(object.list) >= 1) {
  pathways.show <- Reduce(union, lapply(object.list, function(x) x@netP$pathways))
} else {
  pathways.show <- cellchat@netP$pathways
}
cat(paste0("Pathways to visualize: ", length(pathways.show), " (",
           paste(head(pathways.show, 5), collapse = ", "),
           ifelse(length(pathways.show) > 5, "...)", ")"), "\n\n"))

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
      color.use   = CONDITION_COLORS,
      type        = GENE_EXPR_VIZ_TYPE
    ) +
      plot_annotation(
        title = paste0(pw, " - Gene Expression"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      ) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "right",
        axis.text.y     = element_text(size = 15),
        axis.title.y    = element_text(size = 15, angle = 0),
        legend.text     = element_text(size = 17, margin = margin(t = 2, b = 2)),
        legend.title    = element_text(size = 17, margin = margin(b = 4)),
        legend.spacing.y = unit(0.15, "cm")
      ) &
      labs(fill = "Condition", color = "Condition")

    # Show cell type labels only on the bottom panel
    if (!is.null(gg_violin$patches) && length(gg_violin$patches$plots) > 0) {
      n_plots <- length(gg_violin$patches$plots)
      for (i in seq_len(n_plots - 1)) {
        gg_violin$patches$plots[[i]] <- gg_violin$patches$plots[[i]] +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      gg_violin$patches$plots[[n_plots]] <- gg_violin$patches$plots[[n_plots]] +
        theme(
          axis.text.x = element_text(size = 19, angle = 0, hjust = 1, vjust = 1,
                                       margin = margin(t = 4)),
          axis.ticks.x = element_line()
        )
    }

    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(GENE_EXPR_VIZ_TYPE, "_", pw, ".png")),
      plot     = gg_violin,
      width    = GENE_EXPR_VIZ_WIDTH,
      height   = GENE_EXPR_VIZ_HEIGHT,
      dpi      = GENE_EXPR_VIZ_DPI,
      limitsize = FALSE
    )
    cat(paste0("    Saved: ", GENE_EXPR_VIZ_TYPE, "_", pw, ".png\n"))
  }, error = function(e) {
    cat(paste0("    FAILED ", GENE_EXPR_VIZ_TYPE, ": ", e$message, "\n"))
  })

  # --- Dot plot disabled ---
  # tryCatch({
  #   gg_dot <- plotGeneExpression(
  #     cellchat,
  #     signaling   = pw,
  #     split.by    = "datasets",
  #     type        = "dot",
  #     color.use   = CONDITION_COLORS
  #   )
  #   ggplot2::ggsave(
  #     filename = file.path(out_dir, paste0("dot_", pw, ".png")),
  #     plot     = gg_dot,
  #     width    = 14,
  #     height   = 10,
  #     dpi      = 150,
  #     limitsize = FALSE
  #   )
  #   cat(paste0("    Saved: dot_", pw, ".png\n"))
  # }, error = function(e) {
  #   cat(paste0("    FAILED dot: ", e$message, "\n"))
  # })
}

cat("\n--- Procedure 2 Step 15 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
