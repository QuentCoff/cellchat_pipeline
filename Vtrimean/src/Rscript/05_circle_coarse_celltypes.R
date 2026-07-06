#!/usr/bin/env Rscript
# 05_circle_coarse_celltypes.R
# CellChat Procedure 2 Step 8: Circle plots at coarse cell type level
# Aggregates cell populations into 5 coarse groups for simplified network view
#
# Cell type grouping (5 groups):
#   SSC        : SSC
#   Germ       : Spermatocyte, Spermatid
#   Sertoli    : Sertoli
#   Leydig     : Leydig
#   Somatic    : Myoid, Fibroblast, Endothelial
#
# Usage:
#   Rscript 05_circle_coarse_celltypes.R <prefix1> <prefix2>
#
# Example:
#   Rscript 05_circle_coarse_celltypes.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 05_circle_coarse_celltypes.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list <- file.path(base_dir, "Vtrimean", "Result", "data", "merged",
                        paste0(prefix1, "_vs_", prefix2),
                        paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "Vtrimean", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step8")

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

cat("\n=== Step 8-i: Define coarse cell types ===\n")

# Cell types must match the order of levels(object.list[[1]]@idents)
# Expected order: SSC, Spermatocyte, Spermatid, Sertoli, Leydig, Myoid, Fibroblast, Endothelial
cell_levels <- levels(object.list[[1]]@idents)
cat(paste0("Cell type order: ", paste(cell_levels, collapse = ", "), "\n"))

group.cellType <- ifelse(
  cell_levels == "SSC", "SSC",
  ifelse(cell_levels %in% c("Spermatocyte", "Spermatid"), "Germ",
  ifelse(cell_levels == "Sertoli", "Sertoli",
  ifelse(cell_levels == "Leydig", "Leydig",
         "Somatic")))
)
group.cellType <- factor(group.cellType, levels = c("SSC", "Germ", "Sertoli", "Leydig", "Somatic"))

cat(paste0("Grouping:\n"))
for (i in seq_along(cell_levels)) {
  cat(paste0("  ", cell_levels[i], " -> ", group.cellType[i], "\n"))
}

# ============================================
# Step 8-ii: Remerge based on coarse cell types
# ============================================

cat("\n=== Step 8-ii: mergeInteractions + mergeCellChat ===\n")

object.list <- lapply(object.list, function(x) {
  mergeInteractions(x, group.cellType)
})

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# ============================================
# Step 8-iii: Circle plots per dataset (count.merged)
# ============================================

cat("\n=== Step 8-iii: Circle plots per dataset (coarse) ===\n")

weight.max <- getMaxWeight(
  object.list,
  slot.name = c("idents", "net", "net"),
  attribute = c("idents", "count", "count.merged")
)

n <- length(object.list)

png(file.path(out_dir, "circle_coarse_per_dataset_count.png"),
    width = 600 * n, height = 600, res = 120)
par(mfrow = c(1, n), xpd = TRUE)

for (i in 1:n) {
  netVisual_circle(
    object.list[[i]]@net$count.merged,
    weight.scale    = TRUE,
    label.edge      = TRUE,
    edge.weight.max = weight.max[3],
    edge.width.max  = 12,
    title.name      = paste0("Number of interactions - ", names(object.list)[i])
  )
}
dev.off()
cat(paste0("  Saved: circle_coarse_per_dataset_count.png\n"))

# ============================================
# Step 8-iv: Differential count between coarse cell types
# ============================================

cat("\n=== Step 8-iv: Differential number of interactions (coarse) ===\n")

png(file.path(out_dir, "diff_coarse_count.png"), width = 700, height = 600, res = 120)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "count.merged", label.edge = TRUE)
dev.off()
cat(paste0("  Saved: diff_coarse_count.png\n"))

# ============================================
# Step 8-v: Differential strength between coarse cell types
# ============================================

cat("\n=== Step 8-v: Differential interaction strength (coarse) ===\n")

png(file.path(out_dir, "diff_coarse_weight.png"), width = 700, height = 600, res = 120)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight.merged", label.edge = TRUE)
dev.off()
cat(paste0("  Saved: diff_coarse_weight.png\n\n"))

cat("--- Procedure 2 Step 8 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
