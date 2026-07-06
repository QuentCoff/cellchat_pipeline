#!/usr/bin/env Rscript
# 07_robustness.R
# Robustness analysis of CellChat pathway similarity ranking across UMAP runs.
#
# Re-runs only the stochastic part (netEmbedding) N times on the same merged
# object, recomputes per-pathway distances and ranks, then produces:
#   - long table: (iter, pathway, distance, rank)
#   - per-pathway summary: mean/sd/cv/median/IQR of distance, rank stats,
#     frequency of being in top-K, wins/losses/ties from pairwise tests
#   - pairwise Wilcoxon signed-rank tests (paired, BH-corrected) between ALL
#     pathway pairs on the per-iteration distance differences
#
# Output directory is the directory of this script.
#
# Usage:
#   Rscript 07_robustness.R [type] [n_iter]
# Defaults: type=functional, n_iter=100

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
sim_type <- if (length(args) >= 1) args[1] else "functional"
n_iter   <- if (length(args) >= 2) as.integer(args[2]) else 100L

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"
prefix1  <- "healthy"
prefix2  <- "crypto"

# Read merged object from Vtruncated_5Clusters_test_020/Result/data/merged
rdata_merged <- file.path(base_dir, "Vtruncated_5Clusters_test_020", "Result", "data", "merged",
                          paste0(prefix1, "_vs_", prefix2),
                          paste0("cellchat_merged_", prefix1, "_", prefix2, ".RData"))

# Write to Vtruncated_5Clusters_test_020/Result/plot
out_dir <- file.path(base_dir, "Vtruncated_5Clusters_test_020", "Result", "plot",
                     paste0(prefix1, "_vs_", prefix2), "step10bis")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_merged)) {
  stop(paste0("Merged RData not found: ", rdata_merged))
}

cat("=== Loading merged CellChat object ===\n")
load(rdata_merged)  # loads 'cellchat'
cat(paste0("Similarity type: ", sim_type, "\n"))
cat(paste0("Iterations: ", n_iter, "\n\n"))

# ============================================
# Deterministic step: compute pairwise similarity ONCE
# ============================================
cat("=== computeNetSimilarityPairwise (deterministic, once) ===\n")
cellchat <- computeNetSimilarityPairwise(cellchat, type = sim_type)
cat("  Done.\n\n")

comparison.name <- paste(c(1, 2), collapse = "-")

# Helper: extract per-pathway distances from a netEmbedded cellchat
extract_dist <- function(cc) {
  Y <- cc@netP$similarity[[sim_type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  groups <- unique(group)
  data1 <- Y[group == groups[1], , drop = FALSE]
  data2 <- Y[group == groups[2], , drop = FALSE]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))
  pw <- intersect(rownames(data1), rownames(data2))
  d <- sapply(pw, function(p) sqrt(sum((data1[p, ] - data2[p, ])^2)))
  data.frame(pathway = pw, distance = as.numeric(d), stringsAsFactors = FALSE)
}

# ============================================
# Iterate netEmbedding without seed
# ============================================
cat(paste0("=== Looping ", n_iter, " UMAP runs ===\n"))

all_iters <- vector("list", n_iter)
t0 <- Sys.time()
for (i in seq_len(n_iter)) {
  # No set.seed -> each iteration uses a fresh random state
  cc_i <- tryCatch(
    netEmbedding(cellchat, type = sim_type, umap.method = "uwot"),
    error = function(e) { cat(sprintf("  iter %d FAILED: %s\n", i, e$message)); NULL }
  )
  if (is.null(cc_i)) next
  df_i <- extract_dist(cc_i)
  df_i$rank <- rank(-df_i$distance, ties.method = "min")
  df_i$iter <- i
  all_iters[[i]] <- df_i
  if (i %% 10 == 0 || i == n_iter) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  iter %d/%d  (%.1fs elapsed, %.1fs/iter)\n",
                i, n_iter, el, el / i))
  }
}

long_df <- do.call(rbind, all_iters)
write.csv(long_df,
          file.path(out_dir, paste0("iterations_long_", sim_type, ".csv")),
          row.names = FALSE)
cat(paste0("\nSaved long table: iterations_long_", sim_type, ".csv\n"))

# ============================================
# Per-pathway summary
# ============================================
cat("\n=== Per-pathway summary ===\n")
pathways <- sort(unique(long_df$pathway))
K <- length(pathways)
cat(paste0("Number of pathways (intersection of conditions): ", K, "\n"))

# Wide matrix [iter x pathway] of distances and ranks
dist_mat <- matrix(NA_real_, nrow = n_iter, ncol = K,
                   dimnames = list(seq_len(n_iter), pathways))
rank_mat <- dist_mat
for (it in unique(long_df$iter)) {
  sub <- long_df[long_df$iter == it, ]
  dist_mat[it, sub$pathway] <- sub$distance
  rank_mat[it, sub$pathway] <- sub$rank
}

summ <- data.frame(
  pathway       = pathways,
  mean_dist     = colMeans(dist_mat, na.rm = TRUE),
  sd_dist       = apply(dist_mat, 2, sd, na.rm = TRUE),
  median_dist   = apply(dist_mat, 2, median, na.rm = TRUE),
  q25_dist      = apply(dist_mat, 2, quantile, probs = 0.25, na.rm = TRUE),
  q75_dist      = apply(dist_mat, 2, quantile, probs = 0.75, na.rm = TRUE),
  median_rank   = apply(rank_mat, 2, median, na.rm = TRUE),
  q25_rank      = apply(rank_mat, 2, quantile, probs = 0.25, na.rm = TRUE),
  q75_rank      = apply(rank_mat, 2, quantile, probs = 0.75, na.rm = TRUE),
  min_rank      = apply(rank_mat, 2, min, na.rm = TRUE),
  max_rank      = apply(rank_mat, 2, max, na.rm = TRUE),
  freq_top1     = colMeans(rank_mat == 1, na.rm = TRUE),
  freq_top3     = colMeans(rank_mat <= 3, na.rm = TRUE),
  freq_top5     = colMeans(rank_mat <= 5, na.rm = TRUE),
  freq_top10    = colMeans(rank_mat <= 10, na.rm = TRUE),
  stringsAsFactors = FALSE
)
summ$cv_dist <- summ$sd_dist / summ$mean_dist

# Rank frequency table: pathway x rank -> P(rang == k)
rank_freq <- t(apply(rank_mat, 2, function(r) {
  tab <- tabulate(r[!is.na(r)], nbins = K)
  tab / sum(tab)
}))
colnames(rank_freq) <- paste0("rank_", seq_len(K))
rank_freq_df <- data.frame(pathway = rownames(rank_freq), rank_freq,
                           stringsAsFactors = FALSE, check.names = FALSE)
write.csv(rank_freq_df,
          file.path(out_dir, paste0("rank_frequency_", sim_type, ".csv")),
          row.names = FALSE)
cat(paste0("Saved: rank_frequency_", sim_type, ".csv\n"))

# ============================================
# Pairwise paired Wilcoxon tests (signed-rank)
# H0: median( dist[A,i] - dist[B,i] ) == 0
# Paired on iteration i (same UMAP run for both A and B)
# ============================================
cat("\n=== Pairwise Wilcoxon signed-rank tests ===\n")
pmat    <- matrix(NA_real_, K, K, dimnames = list(pathways, pathways))
diffmat <- matrix(NA_real_, K, K, dimnames = list(pathways, pathways))

for (a in seq_len(K - 1)) {
  for (b in (a + 1):K) {
    da <- dist_mat[, a]
    db <- dist_mat[, b]
    ok <- !is.na(da) & !is.na(db)
    if (sum(ok) < 5) next
    diffmat[a, b] <-  median(da[ok] - db[ok])
    diffmat[b, a] <- -diffmat[a, b]
    pv <- tryCatch(
      wilcox.test(da[ok], db[ok], paired = TRUE, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    pmat[a, b] <- pv
    pmat[b, a] <- pv
  }
}

# BH correction over the K*(K-1)/2 unique tests
upper_idx <- which(upper.tri(pmat), arr.ind = TRUE)
pvals <- pmat[upper_idx]
padj  <- p.adjust(pvals, method = "BH")
padj_mat <- matrix(NA_real_, K, K, dimnames = dimnames(pmat))
padj_mat[upper_idx] <- padj
lower_idx <- upper_idx[, c(2, 1)]
padj_mat[lower_idx] <- padj

write.csv(data.frame(pmat,     check.names = FALSE),
          file.path(out_dir, paste0("wilcoxon_pvalues_", sim_type, ".csv")))
write.csv(data.frame(padj_mat, check.names = FALSE),
          file.path(out_dir, paste0("wilcoxon_padj_BH_", sim_type, ".csv")))
write.csv(data.frame(diffmat,  check.names = FALSE),
          file.path(out_dir, paste0("median_diff_", sim_type, ".csv")))
cat("Saved: wilcoxon_pvalues / wilcoxon_padj_BH / median_diff\n")

# "Wins" score per pathway: number of OTHER pathways A is significantly
# greater than (padj < 0.05 AND median(A - B) > 0)
alpha <- 0.05
wins <- sapply(pathways, function(a) {
  sum(padj_mat[a, ] < alpha & diffmat[a, ] > 0, na.rm = TRUE)
})
losses <- sapply(pathways, function(a) {
  sum(padj_mat[a, ] < alpha & diffmat[a, ] < 0, na.rm = TRUE)
})
ties <- (K - 1) - wins - losses

summ$wins   <- wins[summ$pathway]
summ$losses <- losses[summ$pathway]
summ$ties   <- ties[summ$pathway]
summ <- summ[order(-summ$wins, summ$median_rank), ]
write.csv(summ,
          file.path(out_dir, paste0("per_pathway_summary_", sim_type, ".csv")),
          row.names = FALSE)
cat(paste0("Saved: per_pathway_summary_", sim_type, ".csv\n"))

# ============================================
# Plots
# ============================================
cat("\n=== Plots ===\n")

# Order pathways by median rank for plotting (top = best rank)
ord <- summ$pathway[order(summ$median_rank)]
long_df$pathway <- factor(long_df$pathway, levels = rev(ord))

gg_box <- ggplot(long_df, aes(x = pathway, y = distance)) +
  geom_boxplot(outlier.size = 0.5, fill = "#69b3a2", alpha = 0.7) +
  coord_flip() +
  theme_bw(base_size = 10) +
  labs(title = paste0("Distance distribution across ", n_iter,
                      " UMAP runs (", sim_type, ")"),
       x = NULL, y = "Euclidean distance in 2D embedding")
ggsave(file.path(out_dir, paste0("boxplot_distance_", sim_type, ".png")),
       gg_box, width = 8, height = max(4, K * 0.25), dpi = 150, limitsize = FALSE)

gg_rank <- ggplot(long_df, aes(x = pathway, y = rank)) +
  geom_boxplot(outlier.size = 0.5, fill = "#e07b91", alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(breaks = pretty(seq_len(K))) +
  theme_bw(base_size = 10) +
  labs(title = paste0("Rank distribution across ", n_iter,
                      " UMAP runs (", sim_type, ")"),
       x = NULL, y = "Rank (1 = largest distance)")
ggsave(file.path(out_dir, paste0("boxplot_rank_", sim_type, ".png")),
       gg_rank, width = 8, height = max(4, K * 0.25), dpi = 150, limitsize = FALSE)

cat("Saved: boxplot_distance / boxplot_rank\n")
cat("\n--- Done ---\n")
cat(paste0("Output: ", out_dir, "\n"))
