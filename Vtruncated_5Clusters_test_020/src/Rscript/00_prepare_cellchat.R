#!/usr/bin/env Rscript
# 00_prepare_cellchat.R
# CellChat Procedure 1: Data input, Preprocessing & Inference
# Restreint aux 5 clusters : Sertoli, SSC, Leydig, Spermatid, Spermatocyte
# Isolated pipeline step: loads Seurat, filters by cell_type + detailed_group, runs inference, saves RDS
#
# Usage:
#   Rscript 01_prepare_cellchat.R <detailed_group>
#
# Example:
#   Rscript 01_prepare_cellchat.R Healthy
#   Rscript 01_prepare_cellchat.R Crypto       # regroupe Uni-Crypto + Bi-Crypto
#   Rscript 01_prepare_cellchat.R Immuno

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(future)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

# ============================================
# ARGUMENTS
# ============================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 01_prepare_cellchat.R <detailed_group>\n")
}

filter_group <- args[1]       # e.g. "Healthy", "Crypto", "Immuno"

# Crypto regroupe Uni-Crypto + Bi-Crypto
if (tolower(filter_group) == "crypto") {
  filter_values <- c("Uni-Crypto", "Bi-Crypto")
  out_prefix    <- "crypto"
} else {
  filter_values <- filter_group
  out_prefix    <- tolower(filter_group)
}

# Clusters d'intérêt (5)
target_celltypes <- c("Sertoli", "SSC", "Leydig", "Spermatid", "Spermatocyte")

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"
out_dir  <- file.path(base_dir, "Vtruncated_5Clusters_test_020", "Result", "data", "cellchat_prep", out_prefix)
out_rds  <- file.path(out_dir, paste0("cellchat_", out_prefix, ".rds"))

# ============================================
# 1. LOAD AND FILTER SEURAT
# ============================================

cat("=== Loading Seurat object ===\n")
seurat_obj <- readRDS(file.path(base_dir, "data_input", "Final_Annotated_Object_HUMAN.rds"))

cat(paste0("Total cells: ", ncol(seurat_obj), "\n"))
cat("Available cell types:\n")
print(table(seurat_obj$cell_type))
cat("Available detailed_group values:\n")
print(table(seurat_obj$detailed_group))

cat(paste0("\n=== Filtering for cell types: ", paste(target_celltypes, collapse = ", "), " ===\n"))
cells_celltype <- colnames(seurat_obj)[seurat_obj$cell_type %in% target_celltypes]
if (length(cells_celltype) == 0) {
  stop(paste0("No cells found for cell_type in '", paste(target_celltypes, collapse = ", "), "'\n"))
}
seurat_obj <- subset(seurat_obj, cells = cells_celltype)
seurat_obj$cell_type <- droplevels(seurat_obj$cell_type)
cat(paste0("Cells after cell type filtering: ", ncol(seurat_obj), "\n"))
cat("Cell types in subset:\n")
print(table(seurat_obj$cell_type))

cat(paste0("\n=== Filtering detailed_group for '", filter_group, "' ===\n"))
if (length(filter_values) > 1) {
  cat(paste0("  Groups included: ", paste(filter_values, collapse = ", "), "\n"))
}
cells_use <- colnames(seurat_obj)[seurat_obj$detailed_group %in% filter_values]

if (length(cells_use) == 0) {
  stop(paste0("No cells found for detailed_group in '", paste(filter_values, collapse = ", "), "'\n"))
}

seurat_filtered <- subset(seurat_obj, cells = cells_use)
seurat_filtered$cell_type <- droplevels(seurat_filtered$cell_type)
cat(paste0("Cells after group filtering: ", ncol(seurat_filtered), "\n"))
cat("Cell types in filtered data:\n")
print(table(seurat_filtered$cell_type))

# ============================================
# 2. CREATE CELLCHAT OBJECT (Option B: from Seurat)
# ============================================

cat("\n=== Creating CellChat object (Seurat Option B) ===\n")
cellchat <- createCellChat(object = seurat_filtered, group.by = "cell_type", assay = "RNA")

cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
cat(paste0("Groups: ", paste(levels(cellchat@idents), collapse = ", "), "\n"))

# ============================================
# 3. SELECT LIGAND-RECEPTOR DATABASE
# ============================================

cat("\n=== Configuring CellChatDB ===\n")
CellChatDB <- CellChatDB.human  # contient v1 + v2 (3233 interactions au total)

# Step 4 protocol: inspect database
showDatabaseCategory(CellChatDB)
glimpse(CellChatDB$interaction)

# Step 5 Option C: use entire database
cellchat@DB <- CellChatDB

# ============================================
# 4. PRE-PROCESSING
# ============================================

cat("\n=== Pre-processing ===\n")
cellchat <- subsetData(cellchat)

# Step 7: sequential execution (future parallélism blocks on large Seurat objects)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ============================================
# 5. INFERENCE OF CELL-CELL COMMUNICATION
# ============================================

cat("\n=== Inference ===\n")

# Step 11: computeCommunProb with truncatedMean (trim = 0.04)
cellchat <- computeCommunProb(
  cellchat,
  type            = "truncatedMean",
  trim            = 0.20,
  raw.use         = TRUE,
  population.size = TRUE
)

# Step 12: filter by minimum cell count per group
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Step 13: pathway-level inference
cellchat <- computeCommunProbPathway(cellchat)

# Step 14A: aggregate across all groups
cellchat <- aggregateNet(cellchat)

# ============================================
# 6. SAVE RDS (Pause Point)
# ============================================

cat(paste0("\n=== Saving CellChat object ===\n"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(cellchat, file = out_rds)
cat(paste0("Saved: ", out_rds, "\n"))
cat("\n--- Procedure 1 complete ---\n")
