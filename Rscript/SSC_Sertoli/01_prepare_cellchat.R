#!/usr/bin/env Rscript
# 01_prepare_cellchat.R
# Procedure 1 adaptee : subset SSC + Sertoli uniquement, puis CellChat par condition
#
# Usage:
#   Rscript 01_prepare_cellchat.R <detailed_group>
#
# Example:
#   Rscript 01_prepare_cellchat.R Healthy    # -> prefix "healthy"
#   Rscript 01_prepare_cellchat.R Crypto     # -> prefix "crypto" (Uni + Bi)

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(future)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 01_prepare_cellchat.R <detailed_group>\n")
}

filter_group <- args[1]

if (tolower(filter_group) == "crypto") {
  filter_values <- c("Uni-Crypto", "Bi-Crypto")
  out_prefix    <- "crypto"
} else {
  filter_values <- filter_group
  out_prefix    <- tolower(filter_group)
}

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"
out_dir  <- file.path(base_dir, "results", "SSC_Sertoli", "Procedure_1", out_prefix)
out_rds  <- file.path(out_dir, paste0("cellchat_", out_prefix, ".rds"))

# ============================================
# 1. LOAD AND FILTER SEURAT
# ============================================

cat("=== Loading Seurat object ===\n")
seurat_obj <- readRDS(file.path(base_dir, "data_input", "Final_Annotated_Object_HUMAN.rds"))

cat(paste0("Total cells: ", ncol(seurat_obj), "\n"))
cat("Available cell types:\n")
print(table(seurat_obj$cell_type))

cat("\n=== Filtering for cell types: SSC + Sertoli ===\n")
cells_celltype <- colnames(seurat_obj)[seurat_obj$cell_type %in% c("SSC", "Sertoli")]
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
# 2. CREATE CELLCHAT OBJECT
# ============================================

cat("\n=== Creating CellChat object ===\n")
cellchat <- createCellChat(object = seurat_filtered, group.by = "cell_type", assay = "RNA")

cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
cat(paste0("Groups: ", paste(levels(cellchat@idents), collapse = ", "), "\n"))

# ============================================
# 3. SELECT DATABASE
# ============================================

cat("\n=== Configuring CellChatDB ===\n")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

# ============================================
# 4. PRE-PROCESSING
# ============================================

cat("\n=== Pre-processing ===\n")
cellchat <- subsetData(cellchat)

options(future.globals.maxSize = 8 * 1024^3)
future::plan("multisession", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# ============================================
# 5. INFERENCE
# ============================================

cat("\n=== Inference ===\n")
cellchat <- computeCommunProb(
  cellchat,
  type            = "triMean",
  trim            = NULL,
  raw.use         = TRUE,
  population.size = TRUE
)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ============================================
# 6. SAVE RDS
# ============================================

cat("\n=== Saving CellChat object ===\n")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(cellchat, file = out_rds)
cat(paste0("Saved: ", out_rds, "\n"))
cat("\n--- Procedure 1 complete ---\n")
