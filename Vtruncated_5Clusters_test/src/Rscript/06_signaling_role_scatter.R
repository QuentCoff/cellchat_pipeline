#!/usr/bin/env Rscript
# 06_signaling_role_scatter.R
# CellChat Procedure 2 Step 9: Compare major sources and targets in 2D space
# Option A: identify notable changes in sending/receiving signals across datasets
# Option B: identify signaling changes of specific cell populations
#
# Usage:
#   Rscript 06_signaling_role_scatter.R <prefix1> <prefix2>
#
# Example:
#   Rscript 06_signaling_role_scatter.R healthy crypto

suppressPackageStartupMessages({
  library(CellChat)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 06_signaling_role_scatter.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list   <- file.path(base_dir, "Vtruncated_5Clusters_test", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))
rdata_merged <- file.path(base_dir, "Vtruncated_5Clusters_test", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "Vtruncated_5Clusters_test", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step9")

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
    title = names(object.list)[i], weight.MinMax = weight.MinMax)
}

gg_combined <- patchwork::wrap_plots(plots = gg)
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

# object.list already has centrality from Step 9A — rebuild merged from it
cellchat_b <- mergeCellChat(object.list, add.names = names(object.list))

cell_types <- c("Sertoli", "Leydig", "SSC")

for (ct in cell_types) {
  cat(paste0("  Processing: ", ct, "\n"))
  gg_b <- netAnalysis_signalingChanges_scatter(cellchat_b, idents.use = ct)
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("signaling_changes_", gsub(" ", "_", ct), ".png")),
    plot     = gg_b,
    width    = 6,
    height   = 5,
    dpi      = 150
  )
  cat(paste0("    Saved: signaling_changes_", gsub(" ", "_", ct), ".png\n"))
}

cat("\n")

cat("--- Procedure 2 Step 9 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
