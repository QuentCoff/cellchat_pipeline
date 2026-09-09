#!/usr/bin/env Rscript
# 04_circle_per_dataset.R
# CellChat Procedure 2 Step 7: Circle plots normalized across datasets
# Multi-group version: one circle plot per condition on a common scale.
#
# Usage:
#   Rscript 04_circle_per_dataset.R
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
                     merged_dir, "step7")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading object.list ===\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("Conditions: ", paste(names(object.list), collapse = ", "), "\n\n"))

# Keep only colors for the cell types present after subsetting
celltype_colors <- CELLTYPE_COLORS[levels(object.list[[1]]@idents)]
celltype_colors <- celltype_colors[!is.na(celltype_colors)]

# ============================================
# Step 7-i: Compute max weight for common scale
# ============================================

cat("=== Step 7-i: getMaxWeight ===\n")
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
cat(paste0("Max cells per group: ", weight.max[1], "\n"))
cat(paste0("Max interactions: ", weight.max[2], "\n\n"))

# ============================================
# Step 7-ii: Circle plot per dataset (count)
# ============================================

cat("=== Step 7-ii: Circle plots per dataset ===\n")

n <- length(object.list)
png(file.path(out_dir, "circle_per_dataset_count.png"),
    width = 700 * n, height = 700, res = 120)
par(mfrow = c(1, n), xpd = TRUE)

for (i in 1:n) {
  netVisual_circle(
    object.list[[i]]@net$count,
    weight.scale    = TRUE,
    edge.weight.max = weight.max[2],
    edge.width.max  = 8,
    title.name      = NULL,
    color.use       = celltype_colors,
    vertex.label.cex = 1.6
  )
  mtext(paste0("Number of interactions - ", names(object.list)[i]),
        side = 3, line = 0.5, cex = 1.2, font = 2)
}
dev.off()
cat(paste0("  Saved: circle_per_dataset_count.png\n"))

# ============================================
# Circle plots per dataset (weight/strength)
# ============================================

cat("\n=== Circle plots per dataset (strength) ===\n")

weight.max.w <- getMaxWeight(object.list, attribute = c("idents", "weight"))

png(file.path(out_dir, "circle_per_dataset_weight.png"),
    width = 700 * n, height = 700, res = 120)
par(mfrow = c(1, n), xpd = TRUE)

for (i in 1:n) {
  netVisual_circle(
    object.list[[i]]@net$weight,
    weight.scale    = TRUE,
    edge.weight.max = weight.max.w[2],
    edge.width.max  = 6,
    title.name      = NULL,
    color.use       = celltype_colors,
    vertex.label.cex = 1.6
  )
  mtext(paste0("Interaction strength - ", names(object.list)[i]),
        side = 3, line = 0.5, cex = 1.2, font = 2)
}
dev.off()
cat(paste0("  Saved: circle_per_dataset_weight.png\n\n"))

cat("--- Procedure 2 Step 7 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
