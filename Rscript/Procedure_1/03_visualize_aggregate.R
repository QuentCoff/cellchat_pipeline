#!/usr/bin/env Rscript
# 03_visualize_aggregate.R
# Visualisation du réseau agrégé global (circle plots) + séparation par cluster
# Inspiré de descended.R - réseau global "tous pathways confondus"
#
# Usage:
#   Rscript 03_visualize_aggregate.R <prefix>
#
# Example:
#   Rscript 03_visualize_aggregate.R healthy

suppressPackageStartupMessages({
  library(CellChat)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 03_visualize_aggregate.R <prefix>\n")
}

prefix <- tolower(args[1])
base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

in_rds  <- file.path(base_dir, "results", "Procedure_1", prefix, paste0("cellchat_", prefix, ".rds"))
out_dir <- file.path(base_dir, "results", "Procedure_1", prefix, "visualization", "aggregate")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Loading CellChat object ===\n")
cellchat <- readRDS(in_rds)

cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
cat(paste0("Groups: ", paste(levels(cellchat@idents), collapse = ", "), "\n\n"))

# ============================================
# RESEAU AGREGÉ GLOBAL (tous pathways)
# ============================================

cat("=== Generating aggregate network plots ===\n")

# Couleurs par cell type (même palette que descended.R)
celltype_colors <- c(
  "SSC"          = "#417505",
  "Spermatocyte" = "#4A90E2",
  "Spermatid"    = "#D0021B",
  "Sertoli"      = "#F5A623",
  "Leydig"       = "#880FF3",
  "Myoid"        = "#7ED321",
  "Fibroblast"   = "#5A5A5A",
  "Endothelial"  = "#F16B54"
)

# Filtrer les couleurs pour ne garder que les groupes présents
cell_groups <- levels(cellchat@idents)
color_use <- celltype_colors[cell_groups]

# --- 1. Réseau global: Nombre d'interactions ---
png(file.path(out_dir, "network_count.png"), width = 1000, height = 800, res = 120)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Number of interactions",
  edge.width.max = 10,
  color.use     = color_use
)
dev.off()
cat(paste0("  Saved: network_count.png\n"))

# --- 2. Réseau global: Poids / force des interactions ---
png(file.path(out_dir, "network_weight.png"), width = 1000, height = 800, res = 120)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = as.numeric(table(cellchat@idents)),
  weight.scale  = TRUE,
  label.edge    = FALSE,
  title.name    = "Interaction weights/strength",
  edge.width.max = 10,
  color.use     = color_use
)
dev.off()
cat(paste0("  Saved: network_weight.png\n"))

# --- 3. Réseau séparé par cluster (sender) ---
cat("\n=== Generating per-cluster sender plots ===\n")

mat <- cellchat@net$count
n_groups <- nrow(mat)

# Layout dynamique: 2 lignes, ceil(n/2) colonnes
cols <- ceiling(n_groups / 2)
png(file.path(out_dir, "network_per_cluster.png"),
    width = 400 * cols, height = 800, res = 120)
par(mfrow = c(2, cols), xpd = TRUE)

for (i in 1:n_groups) {
  mat2 <- matrix(0, nrow = n_groups, ncol = n_groups, dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight   = as.numeric(table(cellchat@idents)),
    weight.scale    = TRUE,
    edge.weight.max = max(mat),
    edge.width.max  = 12,
    title.name      = rownames(mat)[i],
    color.use       = color_use
  )
}
dev.off()
cat(paste0("  Saved: network_per_cluster.png (", n_groups, " clusters)\n"))

cat(paste0("\n--- Aggregate visualization complete ---\n"))
cat(paste0("Output: ", out_dir, "\n"))
