#!/usr/bin/env Rscript
# 09_bubble_dysfunctional.R
# CellChat Procedure 2 Step 12: Identify dysfunctional signaling
# Pairwise version: loops through all PAIRWISE comparisons defined in config.R.
#
# Usage:
#   Rscript 09_bubble_dysfunctional.R

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# Helper: dynamic plot size based on bubble plot content
# y-axis = L-R pairs, x-axis = source -> target combinations
bubble_size <- function(gg) {
  d <- gg$data
  n_lr     <- length(unique(d$interaction_name_2))
  n_groups <- length(unique(paste(d$source, d$target)))
  list(
    width  = max(8,  0.5  * n_groups + 4),
    height = max(6,  0.22 * n_lr     + 3)
  )
}

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

pair_paths <- function(pair) {
  pair_prefix <- paste(tolower(pair), collapse = "_vs_")
  list(
    prefix = pair_prefix,
    rdata_merged = file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                             MERGED_DIR_NAME, pair_prefix,
                             paste0("cellchat_merged_", pair_prefix, ".RData")),
    rdata_list = file.path(base_dir, PROJECT_NAME, "Result", "data", "merged",
                           MERGED_DIR_NAME, pair_prefix,
                           paste0("cellchat_object.list_", pair_prefix, ".RData")),
    out_dir = file.path(base_dir, PROJECT_NAME, "Result", "plot",
                        MERGED_DIR_NAME, pair_prefix, "step12")
  )
}

run_pair <- function(pair) {
  paths <- pair_paths(pair)
  dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(paths$rdata_merged)) {
    cat(paste0("  SKIPPED: merged RData not found: ", paths$rdata_merged, "\n"))
    return(invisible(NULL))
  }
  if (!file.exists(paths$rdata_list)) {
    cat(paste0("  SKIPPED: object.list RData not found: ", paths$rdata_list, "\n"))
    return(invisible(NULL))
  }

  out_dir <- paths$out_dir

  cat("=== Loading CellChat objects ===\n")
  load(paths$rdata_merged)  # loads 'cellchat'
  load(paths$rdata_list)    # loads 'object.list'

  dataset_names <- names(object.list)
  cat(paste0("Datasets: ", paste(dataset_names, collapse = ", "), "\n"))
  cat(paste0("pos.dataset (ref): ", dataset_names[2], "\n\n"))

# Cell-type indices for Sertoli → Leydig / SSC analysis
idents <- levels(cellchat@idents$joint)
sertoli_idx <- match("Sertoli", idents)
leydig_idx  <- match("Leydig",  idents)
ssc_idx     <- match("SSC",     idents)

if (any(is.na(c(sertoli_idx, leydig_idx, ssc_idx)))) {
  stop("One or more cell types (Sertoli, Leydig, SSC) not found in idents\n")
}

cat(paste0("Cell-type indices — Sertoli: ", sertoli_idx,
           ", Leydig: ", leydig_idx, ", SSC: ", ssc_idx, "\n\n"))

# ============================================
# Step 12A: Bubble plots — communication probabilities
# ============================================

# --- 12A-i: All L-R pairs ---
cat("=== Step 12A-i: Bubble plot — all communications ===\n")

tryCatch({
  gg <- netVisual_bubble(cellchat,
                       sources.use = c(sertoli_idx, ssc_idx),
                       targets.use = c(sertoli_idx, ssc_idx),
                       comparison = c(1, 2), angle.x = 45)
  sz <- bubble_size(gg)
  cat(paste0("  Dynamic size: ", round(sz$width, 1), " x ", round(sz$height, 1), " in\n"))
  ggplot2::ggsave(
    filename = file.path(out_dir, "bubble_all.png"),
    plot     = gg,
    width    = sz$width,
    height   = sz$height,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: bubble_all.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# --- 12A-ii: Up-regulated in dataset 2 ---
cat("\n=== Step 12A-ii: Bubble plot — up-regulated in dataset 2 ===\n")

tryCatch({
  gg <- netVisual_bubble(cellchat,
                       sources.use = c(sertoli_idx, ssc_idx),
                       targets.use = c(sertoli_idx, ssc_idx),
                       comparison = c(1, 2), max.dataset = 2,
                       title.name = paste0("Increased signaling in ", dataset_names[2]),
                       angle.x = 45, remove.isolate = TRUE)
  sz <- bubble_size(gg)
  cat(paste0("  Dynamic size: ", round(sz$width, 1), " x ", round(sz$height, 1), " in\n"))
  ggplot2::ggsave(
    filename = file.path(out_dir, "bubble_up.png"),
    plot     = gg,
    width    = sz$width,
    height   = sz$height,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: bubble_up.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# --- 12A-iii: Down-regulated in dataset 2 (up in dataset 1) ---
cat("\n=== Step 12A-iii: Bubble plot — down-regulated in dataset 2 ===\n")

tryCatch({
  gg <- netVisual_bubble(cellchat,
                       sources.use = c(sertoli_idx, ssc_idx),
                       targets.use = c(sertoli_idx, ssc_idx),
                       comparison = c(1, 2), max.dataset = 1,
                       title.name = paste0("Decreased signaling in ", dataset_names[2]),
                       angle.x = 45, remove.isolate = TRUE)
  sz <- bubble_size(gg)
  cat(paste0("  Dynamic size: ", round(sz$width, 1), " x ", round(sz$height, 1), " in\n"))
  ggplot2::ggsave(
    filename = file.path(out_dir, "bubble_down.png"),
    plot     = gg,
    width    = sz$width,
    height   = sz$height,
    dpi      = 150,
    limitsize = FALSE
  )
  cat("  Saved: bubble_down.png\n")
}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

# ============================================
# Step 12B: Differential expression analysis
# ============================================

cat("\n=== Step 12B: Differential expression analysis ===\n")

pos.dataset   <- dataset_names[2]  # dataset 2 = condition test
features.name <- pos.dataset

cat(paste0("  pos.dataset: ", pos.dataset, "\n"))
cat("  Running identifyOverExpressedGenes (fast mode with presto)...\n")

tryCatch({
  cellchat <- identifyOverExpressedGenes(
    cellchat,
    group.dataset    = "datasets",
    pos.dataset      = pos.dataset,
    features.name    = features.name,
    only.pos         = FALSE,
    thresh.pc        = 0.1,
    thresh.fc        = 0.05,
    do.fast          = TRUE
  )
  cat("  Done.\n")

  # --- Map DEG results onto inferred communications ---
  cat("\n  Mapping DEG onto cell-cell communications...\n")
  net <- netMappingDEG(cellchat, features.name = features.name)
  cat("  Done.\n")

  # --- Up-regulated L-R pairs in pos.dataset ---
  cat(paste0("\n  Extracting up-regulated L-R pairs in ", pos.dataset, "...\n"))
  net.up <- subsetCommunication(
    cellchat, net = net, datasets = pos.dataset,
    sources.use = c("Sertoli", "SSC"),
    targets.use = c("Sertoli", "SSC"),
    ligand.logFC = 0.05, receptor.logFC = NULL
  )
  write.csv(net.up, file.path(out_dir, "net_up.csv"), row.names = FALSE)
  cat(paste0("  Saved: net_up.csv (n=", nrow(net.up), ")\n"))

  # --- Down-regulated L-R pairs in pos.dataset ---
  cat(paste0("\n  Extracting down-regulated L-R pairs in ", pos.dataset, "...\n"))
  net.down <- subsetCommunication(
    cellchat, net = net, datasets = dataset_names[1],
    sources.use = c("Sertoli", "SSC"),
    targets.use = c("Sertoli", "SSC"),
    ligand.logFC = -0.05, receptor.logFC = NULL
  )
  write.csv(net.down, file.path(out_dir, "net_down.csv"), row.names = FALSE)
  cat(paste0("  Saved: net_down.csv (n=", nrow(net.down), ")\n"))

  # --- Optional: extract individual signaling genes ---
  cat("\n  Extracting individual signaling genes (optional)...\n")
  gene.up   <- extractGeneSubsetFromPair(net.up,   cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

  write.csv(gene.up,   file.path(out_dir, "gene_up.csv"),   row.names = FALSE)
  write.csv(gene.down, file.path(out_dir, "gene_down.csv"), row.names = FALSE)
  cat("  Saved: gene_up.csv, gene_down.csv\n")

  # --- Optional: find enriched signaling ---
  cat("\n  Optional: findEnrichedSignaling...\n")
  if (length(gene.up) > 0) {
    df.up <- findEnrichedSignaling(
      object.list[[2]],
      features = gene.up[1:min(10, length(gene.up))],
      idents   = c("Sertoli", "SSC"),
      pattern  = "outgoing"
    )
    write.csv(df.up, file.path(out_dir, "enriched_signaling_up.csv"), row.names = FALSE)
    cat("  Saved: enriched_signaling_up.csv\n")
  }

  # Save DEG cellchat object
  save(cellchat, file = file.path(out_dir, "cellchat_deg.RData"))
  cat("  Saved: cellchat_deg.RData\n")

}, error = function(e) {
  cat(paste0("  FAILED: ", e$message, "\n"))
})

cat("\n--- Procedure 2 Step 12 complete ---\n")
cat(paste0("Output: ", out_dir, "\n"))
}

for (pair in PAIRWISE) {
  cat(paste0("\n##############################################\n"))
  cat(paste0("# Pairwise comparison: ", paste(pair, collapse = " vs "), "\n"))
  cat(paste0("##############################################\n"))
  run_pair(pair)
}
