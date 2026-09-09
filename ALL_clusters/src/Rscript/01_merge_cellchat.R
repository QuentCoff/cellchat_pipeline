#!/usr/bin/env Rscript
# 01_merge_cellchat.R
# CellChat Procedure 2 Step 1-3: Load and merge CellChat objects for comparison
#
# Usage:
#   Rscript 01_merge_cellchat.R
#
# Configuration is read from config.R. It creates a multi-group merged object
# plus all pairwise merged objects defined in PAIRWISE.

suppressPackageStartupMessages({
  library(CellChat)
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

# ============================================
# Helper: merge a subset of conditions
# ============================================

merge_conditions <- function(conditions) {
  prefixes <- tolower(conditions)
  labels <- conditions

  object.list <- list()
  for (i in seq_along(conditions)) {
    rds <- file.path(base_dir, PROJECT_NAME, "Result", "data", "cellchat_prep", prefixes[i],
                     paste0("cellchat_", prefixes[i], ".rds"))
    cat(paste0("[", i, "/", length(conditions), "] Loading: ", rds, "\n"))
    if (!file.exists(rds)) {
      stop(paste0("File not found: ", rds, "\nRun 00_prepare_cellchat.R first.\n"))
    }
    obj <- readRDS(rds)
    cat(paste0("  Groups: ", length(levels(obj@idents)), " | ",
                paste(levels(obj@idents), collapse = ", "), "\n"))

    if (DO_SUBSET) {
      cat("  Subsetting to TARGET_CELLTYPES...\n")
      obj <- subsetCellChat(obj, idents.use = TARGET_CELLTYPES)
      cat(paste0("  Groups after subset: ", length(levels(obj@idents)), " | ",
                  paste(levels(obj@idents), collapse = ", "), "\n"))
    } else {
      cat("  Keeping all cell types (DO_SUBSET = FALSE)\n")
    }

    object.list[[labels[i]]] <- obj
  }

  cat(paste0("Merging: ", paste(labels, collapse = " + "), "\n"))
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))

  cat(paste0("Merged object groups: ", length(levels(cellchat@idents)), "\n"))
  cat(paste0("Merged object datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

  return(list(object.list = object.list, cellchat = cellchat))
}

# ============================================
# Save helper
# ============================================

save_merge <- function(conditions, merged) {
  prefixes <- tolower(conditions)
  merge_prefix <- paste(prefixes, collapse = "_vs_")

  out_dir <- file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                       MERGED_DIR_NAME, merge_prefix)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  rdata_list <- file.path(out_dir, paste0("cellchat_object.list_", merge_prefix, ".RData"))
  rdata_merged <- file.path(out_dir, paste0("cellchat_merged_", merge_prefix, ".RData"))

  cat("=== Saving merged object ===\n")
  e <- new.env()
  e$object.list <- merged$object.list
  e$cellchat <- merged$cellchat

  save(list = "object.list", file = rdata_list, envir = e)
  cat(paste0("RData saved: ", rdata_list, "\n"))

  save(list = "cellchat", file = rdata_merged, envir = e)
  cat(paste0("RData saved: ", rdata_merged, "\n\n"))

  return(invisible(NULL))
}

# ============================================
# 1. Multi-group merge (all conditions)
# ============================================

cat("=== Multi-group merge ===\n")
merged_all <- merge_conditions(CONDITIONS)
save_merge(CONDITIONS, merged_all)

# ============================================
# 2. Pairwise merges (for scripts that only support two groups)
# ============================================

if (length(PAIRWISE) > 0) {
  cat("=== Pairwise merges ===\n")
  for (pair in PAIRWISE) {
    merged_pair <- merge_conditions(pair)
    save_merge(pair, merged_pair)
  }
}

cat("--- Procedure 2 Steps 1-3 complete ---\n")
