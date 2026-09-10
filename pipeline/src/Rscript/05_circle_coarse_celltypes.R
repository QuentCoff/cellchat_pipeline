#!/usr/bin/env Rscript
# 05_circle_coarse_celltypes.R
# CellChat Procedure 2 Step 8: Circle plots at coarse cell type level
#
# What it does:
#   Loads the multi-group object list, merges fine cell types into coarse groups
#   defined in config.R (COARSE_GROUP_MAP), and produces per-dataset circle plots.
#   Differential coarse plots are generated only when exactly two groups are present.
#
# Inputs:
#   - Multi-group object list RData:
#     Result/data/merged/<MERGED_DIR_NAME>/<MERGE_PREFIX>/cellchat_object.list_<MERGE_PREFIX>.RData
#   - Configuration file: config.R
#
# Outputs:
#   - Circle plots in Result/plot/<MERGED_DIR_NAME>/step8/:
#       circle_coarse_per_dataset_count.png
#       diff_coarse_count.png  (if exactly 2 groups)
#       diff_coarse_weight.png (if exactly 2 groups)
#
# Previous step: 01_merge_cellchat.R
# Next step: 06_signaling_role_scatter.R
#
# Usage:
#   Rscript 05_circle_coarse_celltypes.R
#
# Configuration is read from config.R.

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

rdata_list <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                        merged_dir, merge_prefix,
                        paste0("cellchat_object.list_", merge_prefix, ".RData"))

out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "plot",
                     merged_dir, STEP8_DIR)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading object.list ===\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("Conditions: ", paste(names(object.list), collapse = ", "), "\n"))

# ============================================
# Step 8-i: Define coarse cell type groups
# ============================================

cat("\n=== Step 8-i: Define coarse cell type groups ===\n")

# Cell types must match the order of levels(object.list[[1]]@idents).
cell_levels <- levels(object.list[[1]]@idents)
cat(paste0("Cell type order: ", paste(cell_levels, collapse = ", "), "\n"))

# Map each fine cell type to its coarse group using config.R; unknown types go to default.
group.cellType <- ifelse(
  cell_levels %in% names(COARSE_GROUP_MAP),
  COARSE_GROUP_MAP[cell_levels],
  COARSE_GROUP_DEFAULT
)
# Enforce the desired coarse-group ordering for plotting.
group.cellType <- factor(group.cellType, levels = COARSE_GROUP_LEVELS)

cat(paste0("Grouping:\n"))
for (i in seq_along(cell_levels)) {
  cat(paste0("  ", cell_levels[i], " -> ", group.cellType[i], "\n"))
}

# ============================================
# Step 8-ii: Remerge based on coarse cell types
# ============================================

cat("\n=== Step 8-ii: mergeInteractions + mergeCellChat ===\n")

# Collapse fine cell types into coarse groups for every per-condition object.
object.list <- lapply(object.list, function(x) {
  mergeInteractions(x, group.cellType)
})

# Merge the coarse-grained objects into a multi-group CellChat object.
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# ============================================
# Step 8-iii: Circle plots per dataset (count.merged)
# ============================================

cat("\n=== Step 8-iii: Circle plots per dataset (coarse) ===\n")

# Common scale for the coarse-grained count matrices.
weight.max <- getMaxWeight(
  object.list,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "count", "count.merged")
)

n <- length(object.list)

png(file.path(out_dir, "circle_coarse_per_dataset_count.png"),
    width = COARSE_WIDTH * n, height = COARSE_HEIGHT, res = COARSE_RES)
par(mfrow = c(1, n), xpd = TRUE)

for (i in 1:n) {
  netVisual_circle(
    object.list[[i]]@net$count.merged,
    weight.scale    = TRUE,
    label.edge      = TRUE,
    edge.weight.max = weight.max[3],
    edge.width.max  = COARSE_EDGE_WIDTH,
    title.name      = paste0("Number of interactions - ", names(object.list)[i])
  )
}
dev.off()
cat(paste0("  Saved: circle_coarse_per_dataset_count.png\n"))

# ============================================
# Step 8-iv/v: Differential plots only for pairwise merges
# ============================================

# Differential coarse plots are only supported when exactly two groups are present.
if (length(object.list) == 2) {
  cat("\n=== Step 8-iv: Differential number of interactions (coarse) ===\n")

  png(file.path(out_dir, "diff_coarse_count.png"), width = COARSE_WIDTH, height = COARSE_HEIGHT, res = COARSE_RES)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "count.merged", label.edge = TRUE)
  dev.off()
  cat(paste0("  Saved: diff_coarse_count.png\n"))

  cat("\n=== Step 8-v: Differential interaction strength (coarse) ===\n")

  png(file.path(out_dir, "diff_coarse_weight.png"), width = COARSE_WIDTH, height = COARSE_HEIGHT, res = COARSE_RES)
  netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight.merged", label.edge = TRUE)
  dev.off()
  cat(paste0("  Saved: diff_coarse_weight.png\n\n"))
} else {
  cat("\n  Skipped differential coarse plots: only available for exactly two groups.\n\n")
}

cat("--- Procedure 2 Step 8 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
