#!/usr/bin/env Rscript
# 05_circle_per_dataset.R
# Step 7: Circle plots normalized across datasets
#
# Usage:
#   Rscript 05_circle_per_dataset.R <prefix1> <prefix2>

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 05_circle_per_dataset.R <prefix1> <prefix2>\n")
}

prefix1 <- tolower(args[1])
prefix2 <- tolower(args[2])

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_list <- file.path(base_dir, "results", "SSC_Sertoli", "Procedure_2",
                        paste0(prefix1, "_vs_", prefix2),
                        paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "results", "SSC_Sertoli", "Procedure_2",
                     paste0(prefix1, "_vs_", prefix2), "comparison", "step7")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\nRun 02_merge_cellchat.R first.\n"))
}

cat("=== Loading object.list ===\n")
load(rdata_list)
cat(paste0("Conditions: ", paste(names(object.list), collapse = ", "), "\n\n"))

cat("=== Step 7-i: getMaxWeight ===\n")
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
cat(paste0("Max cells per group: ", weight.max[1], "\n"))
cat(paste0("Max interactions: ", weight.max[2], "\n\n"))

# Custom colors: SSC = green, Sertoli = yellow
color.use <- c("SSC" = "#417505", "Sertoli" = "#FFD700")

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
    edge.width.max  = 12,
    title.name      = paste0("Number of interactions - ", names(object.list)[i]),
    color.use       = color.use
  )
}
dev.off()
cat("  Saved: circle_per_dataset_count.png\n")

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
    edge.width.max  = 12,
    title.name      = paste0("Interaction strength - ", names(object.list)[i]),
    color.use       = color.use
  )
}
dev.off()
cat("  Saved: circle_per_dataset_weight.png\n\n")

cat("--- Procedure 2 Step 7 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
