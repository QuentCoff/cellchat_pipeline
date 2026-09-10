#!/usr/bin/env Rscript
# 00_prepare_cellchat.R
# CellChat Procedure 1: Data input, Preprocessing & Inference
#
# What it does:
#   Loads the Seurat object, filters cells by condition (detailed_group), creates
#   one CellChat object per condition, runs preprocessing and communication inference.
#
# Inputs:
#   - Seurat RDS file: data_input/<PREP_INPUT_RDS>
#   - Configuration file: config.R
#
# Outputs:
#   - One RDS per condition:
#     Result/data/cellchat_prep/<condition>/cellchat_<condition>.rds
#
# Previous step: none (this is the first script)
# Next step: 01_merge_cellchat.R
#
# Usage:
#   Rscript 00_prepare_cellchat.R
#
# Configuration is read from config.R.

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(future)
  library(dplyr)
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

# Helper: map a condition name to the metadata values used to filter the Seurat object.
# Crypto is a special case: it combines Uni-Crypto and Bi-Crypto samples.
get_condition_meta <- function(cond) {
  if (tolower(cond) == "crypto") {
    list(filter_values = c("Uni-Crypto", "Bi-Crypto"), out_prefix = "crypto")
  } else {
    list(filter_values = cond, out_prefix = tolower(cond))
  }
}

base_dir <- BASE_DIR  # project root from config.R

# ============================================
# 1. LOAD SEURAT AND DATABASE (once)
# ============================================

cat("=== Loading Seurat object ===\n")
# Read the annotated Seurat object provided by the user.
seurat_obj <- readRDS(file.path(base_dir, "data_input", PREP_INPUT_RDS))

cat(paste0("Total cells: ", ncol(seurat_obj), "\n"))
cat(paste0("Available ", SEURAT_CONDITION_COL, " values:\n"))
# Print the distribution of samples/conditions in the metadata.
print(table(seurat_obj@meta.data[[SEURAT_CONDITION_COL]]))

cat("\n=== Configuring CellChatDB ===\n")
# Select the appropriate ligand-receptor database (human or mouse).
CellChatDB <- switch(PREP_CELLCHAT_DB,
                       human = CellChatDB.human,
                       mouse = CellChatDB.mouse,
                       stop("Unknown PREP_CELLCHAT_DB: ", PREP_CELLCHAT_DB))
showDatabaseCategory(CellChatDB)
glimpse(CellChatDB$interaction)

# ============================================
# 2. PROCESS EACH CONDITION
# ============================================

# Loop over each biological condition and create one CellChat object.
for (filter_group in CONDITIONS) {
  meta <- get_condition_meta(filter_group)
  filter_values <- meta$filter_values   # metadata values to keep
  out_prefix <- meta$out_prefix       # lower-case prefix for output names

  # Output paths: one folder per condition.
  out_dir  <- file.path(base_dir, PROJECT_NAME, "Result", "data", "cellchat_prep", out_prefix)
  out_rds  <- file.path(out_dir, paste0("cellchat_", out_prefix, ".rds"))

  cat(paste0("\n##############################################\n"))
  cat(paste0("# Processing condition: ", filter_group, "\n"))
  cat(paste0("##############################################\n"))

  if (length(filter_values) > 1) {
    cat(paste0("Groups included: ", paste(filter_values, collapse = ", "), "\n"))
  }

  # Select cells whose condition column matches the current condition.
  cells_use <- colnames(seurat_obj)[seurat_obj@meta.data[[SEURAT_CONDITION_COL]] %in% filter_values]

  if (length(cells_use) == 0) {
    stop(paste0("No cells found for ", SEURAT_CONDITION_COL, " in '", paste(filter_values, collapse = ", "), "'\n"))
  }

  # Subset the Seurat object to keep only those cells.
  seurat_filtered <- subset(seurat_obj, cells = cells_use)
  cat(paste0("Cells after condition filtering: ", ncol(seurat_filtered), "\n"))
  cat("Cell types in filtered data:\n")
  print(table(seurat_filtered@meta.data[[SEURAT_CELL_TYPE_COL]]))

  cat("\n=== Creating CellChat object (all cell types) ===\n")
  # Build the CellChat object using the cell-type column defined in config.R.
  cellchat <- createCellChat(object = seurat_filtered, group.by = SEURAT_CELL_TYPE_COL, assay = "RNA")

  cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
  cat(paste0("Groups: ", paste(levels(cellchat@idents), collapse = ", "), "\n"))

  cat("\n=== Using CellChatDB ===\n")
  # Attach the selected ligand-receptor database.
  cellchat@DB <- CellChatDB

  cat("\n=== Pre-processing ===\n")
  # Identify over-expressed genes and interactions before inference.
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  cat("\n=== Inference ===\n")
  # Compute communication probabilities using the method defined in config.R.
  cellchat <- computeCommunProb(
    cellchat,
    type            = PREP_PROB_TYPE,
    trim            = PREP_PROB_TRIM,
    raw.use         = PREP_RAW_USE,
    population.size = PREP_POPULATION_SIZE
  )

  # Filter low-confidence communications and aggregate to pathway level.
  cellchat <- filterCommunication(cellchat, min.cells = PREP_MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cat(paste0("\n=== Saving CellChat object ===\n"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # Save the per-condition CellChat object for the merge step.
  saveRDS(cellchat, file = out_rds)
  cat(paste0("Saved: ", out_rds, "\n"))
}

cat("\n--- Procedure 1 complete ---\n")
