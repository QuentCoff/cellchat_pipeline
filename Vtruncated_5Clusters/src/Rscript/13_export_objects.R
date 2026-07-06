#!/usr/bin/env Rscript
# 13_export_objects.R
# CellChat Procedure 2 Step 16: Export merged CellChat object and object.list
#
# Usage:
#   Rscript 13_export_objects.R <prefix1> <prefix2>
#
# Example:
#   Rscript 13_export_objects.R healthy crypto

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 13_export_objects.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list   <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

rdata_merged <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

rdata_deg    <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "plot",
                          paste0(prefix1, "_vs_", prefix2), "step12",
                          "cellchat_deg.RData")

out_dir <- file.path(base_dir, "Vtruncated_5Clusters", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step16")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Step 16: Export CellChat objects ===\n")

# ============================================
# Load object.list
# ============================================

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("  Loading object.list...\n")
load(rdata_list)  # loads 'object.list'
cat(paste0("    Datasets: ", paste(names(object.list), collapse = ", "), "\n"))

# Save to step16 directory
cat("  Saving object.list...\n")
save(object.list,
     file = file.path(out_dir, paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData")))

# ============================================
# Load merged cellchat (prefer DEG version)
# ============================================

if (file.exists(rdata_deg)) {
  cat("  Loading cellchat (with DEG results)...\n")
  load(rdata_deg)  # loads 'cellchat'
} else if (file.exists(rdata_merged)) {
  cat("  Loading merged cellchat...\n")
  load(rdata_merged)  # loads 'cellchat'
} else {
  stop(paste0("Merged cellchat not found. Run previous steps first.\n"))
}

cat("  Saving cellchat...\n")
save(cellchat,
     file = file.path(out_dir, paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData")))

cat("\n--- Procedure 2 Step 16 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
cat(paste0("  ", paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"), "\n"))
cat(paste0("  ", paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"), "\n"))
