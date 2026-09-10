#!/usr/bin/env Rscript
# 07_net_similarity.R
# CellChat Procedure 2 Step 10: Signaling pathway similarity across conditions
# Pairwise version: functional/structural similarity + joint manifold learning + ranking.
# Loops through all PAIRWISE comparisons defined in config.R.
#
# Usage:
#   Rscript 07_net_similarity.R
#
# Configuration is read from config.R. Set NET_SIM_TYPE to "structural" for structural similarity.

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
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
sim_type <- NET_SIM_TYPE
set.seed(NET_SIM_SEED)

pair_paths <- function(pair) {
  pair_prefix <- paste(tolower(pair), collapse = "_vs_")
  sub_dir <- paste0(pair[1], "vs", pair[2])  # e.g. HealthyvsCrypto
  list(
    prefix = pair_prefix,
    rdata_merged = file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                             MERGED_DIR_NAME, pair_prefix,
                             paste0("cellchat_merged_", pair_prefix, ".RData")),
    out_dir = file.path(base_dir, PROJECT_NAME, "Result", "plot",
                        MERGED_DIR_NAME, STEP10_DIR, sub_dir)
  )
}

run_pair <- function(pair) {
  paths <- pair_paths(pair)
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(paths$rdata_merged)) {
    cat(paste0("  SKIPPED: merged RData not found: ", paths$rdata_merged, "\n"))
    return(invisible(NULL))
  }

  cat("=== Loading merged CellChat object ===\n")
  load(paths$rdata_merged)  # loads 'cellchat'
  cat(paste0("Datasets: ", paste(unique(cellchat@meta$datasets), collapse = ", "), "\n"))
  cat(paste0("Similarity type: ", sim_type, "\n\n"))

  pair_title <- paste(pair, collapse = " vs ")

  # Step 10-i: Compute network similarity
  cat("=== Step 10-i: computeNetSimilarityPairwise ===\n")
  cellchat <- computeNetSimilarityPairwise(cellchat, type = sim_type)
  cat("  Done.\n\n")

  # Step 10-ii: Joint manifold learning
  cat("=== Step 10-ii: netEmbedding ===\n")
  cellchat <- netEmbedding(cellchat, type = sim_type, umap.method = NET_SIM_UMAP_METHOD)
  cat("  Done.\n\n")

  # Step 10-iii: Joint clustering
  cat("=== Step 10-iii: netClustering ===\n")
  cellchat <- netClustering(cellchat, type = sim_type)
  cat("  Done.\n\n")

  # Step 10-iv: netVisual_embeddingPairwise
  cat("=== Step 10-iv: netVisual_embeddingPairwise ===\n")
  tryCatch({
    gg_emb <- netVisual_embeddingPairwise(
      cellchat, type = sim_type,
      label.size = 5,
      dot.size   = c(1.5, 5),
      dot.alpha  = 0.8
    ) +
      ggtitle(paste0("Signaling pathway similarity embedding\n", pair_title)) +
      theme(
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title   = element_text(size = 16),
        axis.text    = element_text(size = 13),
        legend.text  = element_text(size = 13),
        legend.title = element_text(size = 14)
      )

    # Force ggrepel to display more labels (avoid removal due to overlap)
    text_idx <- which(sapply(gg_emb$layers, function(l) inherits(l$geom, "GeomTextRepel")))
    if (length(text_idx) > 0) {
      for (idx in text_idx) {
        gg_emb$layers[[idx]]$params$max.overlaps <- Inf
        gg_emb$layers[[idx]]$params$force       <- 0.5
        gg_emb$layers[[idx]]$params$box.padding <- unit(0.1, "lines")
      }
    }

    ggplot2::ggsave(
      filename = file.path(paths$out_dir, paste0("embedding_", sim_type, ".png")),
      plot     = gg_emb,
      width    = NET_SIM_WIDTH,
      height   = NET_SIM_HEIGHT,
      dpi      = NET_SIM_DPI
    )
    cat(paste0("  Saved: embedding_", sim_type, ".png\n\n"))
  }, error = function(e) {
    cat(paste0("  FAILED: ", e$message, "\n\n"))
  })

  # Step 10-v: netVisual_embeddingPairwiseZoomIn
  cat("=== Step 10-v: netVisual_embeddingPairwiseZoomIn ===\n")
  tryCatch({
    gg_zoom <- netVisual_embeddingPairwiseZoomIn(cellchat, type = sim_type, nCol = 2) +
      patchwork::plot_annotation(
        title = paste0("Signaling pathway similarity zoom-in\n", pair_title),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )
    ggplot2::ggsave(
      filename = file.path(paths$out_dir, paste0("embedding_zoomin_", sim_type, ".png")),
      plot     = gg_zoom,
      width    = NET_SIM_WIDTH,
      height   = NET_SIM_HEIGHT,
      dpi      = NET_SIM_DPI
    )
    cat(paste0("  Saved: embedding_zoomin_", sim_type, ".png\n\n"))
  }, error = function(e) {
    cat(paste0("  FAILED: ", e$message, "\n\n"))
  })

  # Step 10-vi: rankSimilarity
  cat("=== Step 10-vi: rankSimilarity ===\n")

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

  write.csv(rank_df, file.path(paths$out_dir, paste0("rank_similarity_data_", sim_type, ".csv")), row.names = FALSE)
  cat(paste0("  Saved: rank_similarity_data_", sim_type, ".csv\n"))

  gg_rank <- rankSimilarity(cellchat, slot.name = "netP", type = sim_type,
                            comparison1 = c(1, 2), comparison2 = c(1, 2)) +
    ggtitle(paste0("Pathway similarity ranking\n", pair_title)) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 12),
      axis.title.x = element_text(size = 14)
    )
  ggplot2::ggsave(
    filename = file.path(paths$out_dir, paste0("rank_similarity_", sim_type, ".png")),
    plot     = gg_rank,
    width    = NET_SIM_WIDTH,
    height   = NET_SIM_HEIGHT,
    dpi      = NET_SIM_DPI
  )
  cat(paste0("  Saved: rank_similarity_", sim_type, ".png\n\n"))

  # Save cellchat with similarity/clustering for downstream river plot
  save(cellchat, file = file.path(paths$out_dir, paste0("cellchat_clustered_", sim_type, ".RData")))
  cat(paste0("  Saved: cellchat_clustered_", sim_type, ".RData\n\n"))

  cat("--- Procedure 2 Step 10 complete ---\n")
  cat(paste0("Output: ", paths$out_dir, "\n"))
}

for (pair in PAIRWISE) {
  cat(paste0("\n##############################################\n"))
  cat(paste0("# Pairwise comparison: ", paste(pair, collapse = " vs "), "\n"))
  cat(paste0("##############################################\n"))
  run_pair(pair)
}
