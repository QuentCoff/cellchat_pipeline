#!/usr/bin/env Rscript
# 06_signaling_role_scatter.R
# CellChat Procedure 2 Step 9: Compare major sources and targets in 2D space
# Multi-group version for 9A; 9B is only generated when exactly two groups are present.
#
# Usage:
#   Rscript 06_signaling_role_scatter.R
#
# Configuration is read from config.R.

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
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

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))
rdata_merged <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                          merged_dir, merge_prefix,
                          paste0("cellchat_merged_", merge_prefix, ".RData"))

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, "step9")

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
cat(paste0("Conditions: ", paste(names(object.list), collapse = ", "), "\n\n"))

# Keep only colors for the cell types present after subsetting
celltype_colors <- CELLTYPE_COLORS[levels(object.list[[1]]@idents)]
celltype_colors <- celltype_colors[!is.na(celltype_colors)]

# ============================================
# Step 9A: Notable changes in sending/receiving
# ============================================

cat("=== Step 9A: Signaling role scatter (all cell populations) ===\n")

# Pre-requisite: compute network centrality scores for each object
cat("=== Pre-processing: netAnalysis_computeCentrality ===\n")
object.list <- lapply(object.list, function(x) {
  netAnalysis_computeCentrality(x, slot.name = "netP")
})
cat("  Done.\n\n")

# (i) Compute min/max number of interactions across all datasets
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))
cat(paste0("Interaction range: [", weight.MinMax[1], ", ", weight.MinMax[2], "]\n\n"))

# (ii) Scatter plot per dataset
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
    title = paste0("Signaling Role - ", names(object.list)[i]),
    weight.MinMax = weight.MinMax,
    color.use = celltype_colors) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 13, colour = "black"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# Use the same x/y axis limits across all scatter plots
get_axis_range <- function(built, axis = "x") {
  scales <- if (axis == "x") built$layout$panel_scales_x else built$layout$panel_scales_y
  if (!is.null(scales[[1]]$range$range)) {
    return(scales[[1]]$range$range)
  }
  # Fallback for newer ggplot2 versions
  params <- built$layout$panel_params[[1]]
  if (axis == "x") params$x.range else params$y.range
}

x_ranges <- lapply(gg, function(g) get_axis_range(ggplot2::ggplot_build(g), "x"))
y_ranges <- lapply(gg, function(g) get_axis_range(ggplot2::ggplot_build(g), "y"))
x_limits <- c(min(sapply(x_ranges, `[`, 1)), max(sapply(x_ranges, `[`, 2)))
y_limits <- c(min(sapply(y_ranges, `[`, 1)), max(sapply(y_ranges, `[`, 2)))

for (i in seq_along(gg)) {
  gg[[i]] <- gg[[i]] + ggplot2::coord_cartesian(xlim = x_limits, ylim = y_limits)
}

gg_combined <- patchwork::wrap_plots(plots = gg) +
  plot_annotation(
    title = "Comparison of Signaling Roles",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
    )
  )
ggplot2::ggsave(
  filename = file.path(out_dir, "signaling_role_scatter.png"),
  plot     = gg_combined,
  width    = 5 * length(object.list),
  height   = 5,
  dpi      = 150
)
cat(paste0("  Saved: signaling_role_scatter.png\n\n"))

# ============================================
# Step 9B: Signaling changes of specific cell populations
# ============================================

cat("=== Step 9B: Signaling changes for specific cell types ===\n")

cell_types <- c("Sertoli", "Leydig", "SSC")

# If the loaded object already contains exactly two groups, run 9B directly.
# Otherwise, run 9B for each pairwise comparison defined in config.R.
if (length(object.list) == 2) {
  pairwise_runs <- list(
    list(prefix = paste(tolower(names(object.list)), collapse = "_vs_"), obj_list = object.list)
  )
} else {
  pairwise_runs <- lapply(PAIRWISE, function(pair) {
    pair_prefix <- paste(tolower(pair), collapse = "_vs_")
    rdata_pair <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                            MERGED_DIR_NAME, pair_prefix,
                            paste0("cellchat_object.list_", pair_prefix, ".RData"))
    if (!file.exists(rdata_pair)) {
      cat(paste0("  SKIPPED: pairwise object list not found: ", rdata_pair, "\n"))
      return(NULL)
    }
    load(rdata_pair)  # loads 'object.list'
    obj_list <- lapply(object.list, netAnalysis_computeCentrality, slot.name = "netP")
    list(prefix = pair_prefix, obj_list = obj_list)
  })
}

for (run in pairwise_runs) {
  if (is.null(run)) next

  pair_dir <- file.path(out_dir, run$prefix)
  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  cat(paste0("  Pairwise comparison: ", run$prefix, "\n"))
  cellchat_b <- mergeCellChat(run$obj_list, add.names = names(run$obj_list))

  pair_colors <- unname(c(
    "grey10",
    CONDITION_COLORS[names(run$obj_list)[1]],
    CONDITION_COLORS[names(run$obj_list)[2]]
  ))

  for (ct in cell_types) {
    cat(paste0("    Processing: ", ct, "\n"))
    gg_b <- netAnalysis_signalingChanges_scatter(cellchat_b, idents.use = ct, color.use = pair_colors)
    ggplot2::ggsave(
      filename = file.path(pair_dir, paste0("signaling_changes_", gsub(" ", "_", ct), ".png")),
      plot     = gg_b,
      width    = 7,
      height   = 6,
      dpi      = 150
    )
    cat(paste0("      Saved: ", run$prefix, "/signaling_changes_", gsub(" ", "_", ct), ".png\n"))
  }
}

cat("\n")

cat("--- Procedure 2 Step 9 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
