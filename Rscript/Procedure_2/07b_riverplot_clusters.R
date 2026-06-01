#!/usr/bin/env Rscript
# 07b_riverplot_clusters.R
# CellChat Procedure 2 Step 10 (extension):
# Alluvial plot: sender cell type -> receiver cell type -> cluster (from step 10)
#
# For each pathway in each dataset, we look at all (sender, receiver) edges with
# pathway prob > quantile threshold, and link them to the cluster of that pathway.
#
# Usage:
#   Rscript 07b_riverplot_clusters.R <prefix1> <prefix2> [type] [quantile] [top_pathways]
#
# Arguments:
#   type:         "functional" (default) or "structural"
#   quantile:     probability quantile threshold per pathway (default: 0.75)
#   top_pathways: keep only the top N pathways per cluster by total prob (default: all = 0)
#
# Example:
#   Rscript 07b_riverplot_clusters.R healthy crypto functional 0.75 0

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 07b_riverplot_clusters.R <prefix1> <prefix2> [type] [role] [top_n]\n")
}

prefix1      <- tolower(args[1])
prefix2      <- tolower(args[2])
sim_type     <- if (length(args) >= 3) args[3] else "functional"
prob_quant   <- if (length(args) >= 4) as.numeric(args[4]) else 0.75
top_pathways <- if (length(args) >= 5) as.integer(args[5]) else 0

base_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

step10_dir <- file.path(base_dir, "results", "Procedure_2",
                        paste0(prefix1, "_vs_", prefix2), "comparison", "step10")

rdata_clustered <- file.path(step10_dir, paste0("cellchat_clustered_", sim_type, ".RData"))
rdata_list      <- file.path(base_dir, "results", "Procedure_2",
                             paste0(prefix1, "_vs_", prefix2),
                             paste0("cellchat_object.list_", prefix1, "_", prefix2, ".RData"))

out_dir <- step10_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rdata_clustered)) {
  stop(paste0("Clustered RData not found: ", rdata_clustered,
              "\nRun 07_net_similarity.R first.\n"))
}
if (!file.exists(rdata_list)) {
  stop(paste0("object.list RData not found: ", rdata_list, "\n"))
}

cat("=== Loading clustered CellChat object ===\n")
load(rdata_clustered)  # 'cellchat'
load(rdata_list)       # 'object.list'

dataset_names <- names(object.list)
cat(paste0("Datasets: ", paste(dataset_names, collapse = ", "), "\n"))
cat(paste0("Similarity type: ", sim_type, " | quantile: ", prob_quant,
           " | top_pathways: ", top_pathways, "\n\n"))

# ============================================
# 1. Extract pathway -> cluster mapping
# ============================================

comparison_name <- paste(seq_along(object.list), collapse = "-")
group_vec <- cellchat@netP$similarity[[sim_type]]$group[[comparison_name]]

if (is.null(group_vec)) {
  stop(paste0("No clustering found at @netP$similarity$", sim_type,
              "$group$", comparison_name, ". Was netClustering() run?\n"))
}

# group_vec names look like: "PathwayName--DatasetName"
pw_cluster_df <- data.frame(
  key     = names(group_vec),
  pathway = sub("--.*", "",  names(group_vec)),
  dataset = sub(".*--", "",  names(group_vec)),
  cluster = paste0("C", as.integer(group_vec)),
  stringsAsFactors = FALSE
)

cat(paste0("  Pathways with cluster assignment: ", nrow(pw_cluster_df), "\n"))
cat(paste0("  Number of clusters: ", length(unique(pw_cluster_df$cluster)), "\n\n"))

# ============================================
# 2. Build (sender, receiver, cluster, dataset) edges from pathway prob arrays
# ============================================

# Optionally restrict to top N pathways per cluster (by total prob across datasets)
if (top_pathways > 0) {
  pw_prob_total <- sapply(seq_len(nrow(pw_cluster_df)), function(k) {
    pw <- pw_cluster_df$pathway[k]; ds <- pw_cluster_df$dataset[k]
    obj <- object.list[[ds]]
    if (!(pw %in% obj@netP$pathways)) return(0)
    sum(obj@netP$prob[, , which(obj@netP$pathways == pw)])
  })
  pw_cluster_df$total_prob <- pw_prob_total
  pw_cluster_df <- pw_cluster_df %>%
    group_by(cluster) %>%
    arrange(desc(total_prob)) %>%
    slice_head(n = top_pathways) %>%
    ungroup() %>%
    as.data.frame()
  cat(paste0("  Restricted to top ", top_pathways,
             " pathways per cluster -> ", nrow(pw_cluster_df), " rows\n"))
}

rows <- list()
for (k in seq_len(nrow(pw_cluster_df))) {
  pw  <- pw_cluster_df$pathway[k]
  ds  <- pw_cluster_df$dataset[k]
  cl  <- pw_cluster_df$cluster[k]
  obj <- object.list[[ds]]
  if (!(pw %in% obj@netP$pathways)) next
  mat <- obj@netP$prob[, , which(obj@netP$pathways == pw)]
  if (all(mat == 0)) next

  thresh <- max(quantile(mat[mat > 0], prob_quant), .Machine$double.eps)
  idx <- which(mat >= thresh, arr.ind = TRUE)
  if (nrow(idx) == 0) next

  rows[[length(rows) + 1]] <- data.frame(
    sender   = rownames(mat)[idx[, 1]],
    receiver = colnames(mat)[idx[, 2]],
    pathway  = pw,
    cluster  = cl,
    dataset  = ds,
    weight   = mat[idx],
    stringsAsFactors = FALSE
  )
}

if (length(rows) == 0) {
  stop("No (sender, receiver, cluster) edges could be built.\n")
}

df_full <- do.call(rbind, rows)
cat(paste0("  Total edges (with pathway): ", nrow(df_full), "\n"))

# Reshape: melt sender/receiver into a single 'celltype' column
# (each edge contributes its weight twice: once for sender, once for receiver)
df_send <- df_full %>%
  transmute(celltype = sender,   pathway, cluster, dataset, weight)
df_recv <- df_full %>%
  transmute(celltype = receiver, pathway, cluster, dataset, weight)

df <- rbind(df_send, df_recv) %>%
  group_by(celltype, pathway, cluster, dataset) %>%
  summarise(weight = sum(weight), .groups = "drop") %>%
  as.data.frame()

cat(paste0("  Aggregated (celltype, pathway, cluster, dataset) rows: ", nrow(df), "\n"))

write.csv(df_full, file.path(out_dir, paste0("riverplot_data_", sim_type, "_full.csv")),
          row.names = FALSE)
write.csv(df, file.path(out_dir, paste0("riverplot_data_", sim_type, "_aggregated.csv")),
          row.names = FALSE)
cat(paste0("  Saved: riverplot_data_", sim_type, "_full.csv (per-edge)\n"))
cat(paste0("  Saved: riverplot_data_", sim_type, "_aggregated.csv (celltype->pathway->cluster)\n\n"))

# ============================================
# 3. Build alluvial plot — one per dataset + combined
# ============================================

make_plot <- function(d, title_str) {
  # Floor weights so every stratum has a minimum visible height
  # (20% of max weight in the panel)
  min_w <- 0.20 * max(d$weight)
  d$weight_plot <- pmax(d$weight, min_w)

  ggplot(d,
         aes(axis1 = celltype, axis2 = pathway, axis3 = cluster,
             y = weight_plot)) +
    geom_alluvium(aes(fill = cluster), alpha = 0.7, width = 1/4,
                  show.legend = FALSE) +
    geom_stratum(width = 1/4, fill = "grey90", color = "grey40",
                 show.legend = FALSE) +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 3.2) +
    scale_x_discrete(limits = c("Cell type", "Pathway", "Cluster"),
                     expand = c(0.08, 0.08)) +
    labs(title = title_str, y = "Aggregated communication probability",
         fill = "Cluster") +
    theme_minimal(base_size = 12) +
    theme(panel.grid     = element_blank(),
          axis.text.y    = element_blank(),
          axis.ticks.y   = element_blank(),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

# Per dataset
for (ds in dataset_names) {
  d_sub <- df[df$dataset == ds, ]
  if (nrow(d_sub) == 0) next
  gg <- make_plot(d_sub, paste0("River plot — ", ds, " (", sim_type, ", q=", prob_quant, ")"))
  out_png <- file.path(out_dir, paste0("riverplot_", sim_type, "_", ds, ".png"))
  ggsave(out_png, gg, width = 12, height = max(6, nrow(d_sub) * 0.08), dpi = 150,
         bg = "white", limitsize = FALSE)
  cat(paste0("  Saved: ", basename(out_png), "\n"))
}

# Combined (faceted)
gg_all <- make_plot(df, paste0("River plot — all datasets (", sim_type, ", q=", prob_quant, ")")) +
  facet_wrap(~ dataset, scales = "free_y", ncol = 1)
out_png <- file.path(out_dir, paste0("riverplot_", sim_type, "_combined.png"))
ggsave(out_png, gg_all, width = 12, height = max(8, nrow(df) * 0.06), dpi = 150,
       bg = "white", limitsize = FALSE)
cat(paste0("  Saved: ", basename(out_png), "\n\n"))

cat("--- River plot complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
