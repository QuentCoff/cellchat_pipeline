#!/usr/bin/env Rscript
# 01_merge_cellchat.R
# CellChat Procedure 2 Step 1-3: Load and merge CellChat objects for comparison
#
# Usage:
#   Rscript 01_merge_cellchat.R <prefix1> <prefix2> [label1] [label2]
#
# Example:
#   Rscript 01_merge_cellchat.R healthy crypto
#   Rscript 01_merge_cellchat.R healthy crypto Healthy Crypto

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 01_merge_cellchat.R <prefix1> <prefix2> [label1] [label2]\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

# Optional labels for the conditions (default: same as prefix)
label1 <- ifelse(length(args) >= 3, args[3], args[1])
label2 <- ifelse(length(args) >= 4, args[4], args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# Read from Vtruncated_test/Result/data/cellchat_prep
rds1 <- file.path(base_dir, "Vtruncated_test", "Result", "data", "cellchat_prep", prefix1, paste0("cellchat_", prefix1, ".rds"))
rds2 <- file.path(base_dir, "Vtruncated_test", "Result", "data", "cellchat_prep", prefix2, paste0("cellchat_", prefix2, ".rds"))

# Write merged object to Vtruncated_test/Result/data/merged
out_dir <- file.path(base_dir, "Vtruncated_test", "Result", "data", "merged", paste0(prefix1, "_vs_", prefix2))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================
# Step 1: Load each CellChat object
# ============================================

cat("=== Loading CellChat objects ===\n")

cat(paste0("[1/2] Loading: ", rds1, "\n"))
if (!file.exists(rds1)) {
  stop(paste0("File not found: ", rds1, "\nRun 01_prepare_cellchat.R first.\n"))
}
cellchat.A <- readRDS(rds1)
cat(paste0("  Groups: ", length(levels(cellchat.A@idents)), " | ",
           paste(levels(cellchat.A@idents), collapse = ", "), "\n"))

cat(paste0("[2/2] Loading: ", rds2, "\n"))
if (!file.exists(rds2)) {
  stop(paste0("File not found: ", rds2, "\nRun 01_prepare_cellchat.R first.\n"))
}
cellchat.B <- readRDS(rds2)
cat(paste0("  Groups: ", length(levels(cellchat.B@idents)), " | ",
           paste(levels(cellchat.B@idents), collapse = ", "), "\n\n"))

# ============================================
# Step 2: (Optional) updateCellChat
# ============================================

# Skip — we use CellChat >= 1.6.0

cat("=== Step 2 (Optional): updateCellChat ===\n")
cat("  Skipped — CellChat >= 1.6.0 assumed.\n\n")

# ============================================
# Step 3: Merge CellChat objects
# ============================================

cat("=== Step 3: mergeCellChat ===\n")

object.list <- list()
object.list[[label1]] <- cellchat.A
object.list[[label2]] <- cellchat.B

cat(paste0("Merging: ", label1, " + ", label2, "\n"))
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cat(paste0("Merged object groups: ", length(levels(cellchat@idents)), "\n"))
cat(paste0("Merged object datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n\n"))

# ============================================
# Step 4: Save merged object
# ============================================

cat("=== Saving merged object ===\n")

rdata_list <- file.path(out_dir, paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))
rdata_merged <- file.path(out_dir, paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

save(object.list, file = rdata_list)
cat(paste0("RData saved: ", rdata_list, "\n"))

save(cellchat, file = rdata_merged)
cat(paste0("RData saved: ", rdata_merged, "\n\n"))

cat("--- Procedure 2 Steps 1-3 complete ---\n")
