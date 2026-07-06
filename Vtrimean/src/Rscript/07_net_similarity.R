#!/usr/bin/env Rscript
# 07_net_similarity.R
# CellChat Procedure 2 Step 10: Signaling pathway similarity across conditions
# Functional or structural similarity + joint manifold learning + ranking
#
# Usage:
#   Rscript 07_net_similarity.R <prefix1> <prefix2> [type]
#
# Arguments:
#   type: "functional" (default) or "structural"
#
# Example:
#   Rscript 07_net_similarity.R healthy crypto
#   Rscript 07_net_similarity.R healthy crypto structural

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(6)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 07_net_similarity.R <prefix1> <prefix2> [type]\n")
}

prefix1  <- tolower(args[1])
prefix2  <- tolower(args[2])
sim_type <- if (length(args) >= 3) args[3] else "functional"

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

rdata_merged <- file.path(base_dir, "Vtrimean", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

out_dir <- file.path(base_dir, "Vtrimean", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step10")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged, "\nRun 01_merge_cellchat.R first.\n"))
}

cat("=== Loading merged CellChat object ===\n")
load(rdata_merged)  # loads 'cellchat'
cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n"))
cat(paste0("Similarity type: ", sim_type, "\n\n"))

# ============================================
# Step 10-i: Compute network similarity
# ============================================

cat("=== Step 10-i: computeNetSimilarityPairwise ===\n")
cellchat <- computeNetSimilarityPairwise(cellchat, type = sim_type)
cat("  Done.\n\n")

# ============================================
# Step 10-ii: Joint manifold learning
# ============================================

cat("=== Step 10-ii: netEmbedding ===\n")
cellchat <- netEmbedding(cellchat, type = sim_type, umap.method = "uwot")
cat("  Done.\n\n")

# ============================================
# Step 10-iii: Joint clustering
# ============================================

cat("=== Step 10-iii: netClustering ===\n")
cellchat <- netClustering(cellchat, type = sim_type)
cat("  Done.\n\n")

# ============================================
# Step 10-iv: netVisual_embeddingPairwise (version merged objects)
# ============================================
# netVisual_embedding hardcode 'comparison = "single"' donc inutilisable
# sur merged object — utiliser la version Pairwise.

cat("=== Step 10-iv: netVisual_embeddingPairwise ===\n")

tryCatch({
  gg_emb <- netVisual_embeddingPairwise(cellchat, type = sim_type, label.size = 3.5)
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("embedding_", sim_type, ".png")),
    plot     = gg_emb,
    width    = 10,
    height   = 8,
    dpi      = 150
  )
  cat(paste0("  Saved: embedding_", sim_type, ".png\n\n"))
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n\n"))
})

# ============================================
# Step 10-v: (Optional) Zoom in each group — version Pairwise
# ============================================

cat("=== Step 10-v: netVisual_embeddingPairwiseZoomIn ===\n")

tryCatch({
  gg_zoom <- netVisual_embeddingPairwiseZoomIn(cellchat, type = sim_type, nCol = 2)
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0("embedding_zoomin_", sim_type, ".png")),
    plot     = gg_zoom,
    width    = 10,
    height   = 8,
    dpi      = 150
  )
  cat(paste0("  Saved: embedding_zoomin_", sim_type, ".png\n\n"))
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n\n"))
})

# ============================================
# Step 10-vi: rankSimilarity — ACTIVÉ
# ============================================

cat("=== Step 10-vi: rankSimilarity ===\n")

# --- Extract numeric ranking data (same logic as rankSimilarity) ---
comparison.name <- paste(c(1, 2), collapse = "-")
Y <- cellchat@netP$similarity[[sim_type]]$dr[[comparison.name]]
group <- sub(".*--", "", rownames(Y))
comparison2.name <- unique(group)
data1 <- Y[group %in% comparison2.name[1], , drop = FALSE]
data2 <- Y[group %in% comparison2.name[2], , drop = FALSE]
rownames(data1) <- sub("--.*", "", rownames(data1))
rownames(data2) <- sub("--.*", "", rownames(data2))
pathway.show <- as.character(intersect(rownames(data1), rownames(data2)))
data1 <- data1[pathway.show, , drop = FALSE]
data2 <- data2[pathway.show, , drop = FALSE]

dist <- sapply(seq_len(nrow(data1)), function(i) {
  sqrt(sum((data1[i, ] - data2[i, ])^2))
})

rank_df <- data.frame(
  pathway = pathway.show,
  distance = dist,
  stringsAsFactors = FALSE
)
rank_df <- rank_df[order(rank_df$distance, decreasing = TRUE), ]
rownames(rank_df) <- NULL
rank_df$rank <- seq_len(nrow(rank_df))

write.csv(rank_df, file.path(out_dir, paste0("rank_similarity_data_", sim_type, ".csv")), row.names = FALSE)
cat(paste0("  Saved: rank_similarity_data_", sim_type, ".csv\n"))

# --- Plot ---
gg_rank <- rankSimilarity(cellchat, slot.name = "netP", type = sim_type,
                          comparison1 = c(1, 2), comparison2 = c(1, 2))
ggplot2::ggsave(
  filename = file.path(out_dir, paste0("rank_similarity_", sim_type, ".png")),
  plot     = gg_rank,
  width    = 8,
  height   = 6,
  dpi      = 150
)
cat(paste0("  Saved: rank_similarity_", sim_type, ".png\n\n"))

# Save cellchat with similarity/clustering for downstream river plot
save(cellchat, file = file.path(out_dir, paste0("cellchat_clustered_", sim_type, ".RData")))
cat(paste0("  Saved: cellchat_clustered_", sim_type, ".RData\n\n"))

cat("--- Procedure 2 Step 10 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
