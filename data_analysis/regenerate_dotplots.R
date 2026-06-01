#!/usr/bin/env Rscript
# Regenere les dotplots pathways x cell types via la fonction native CellChat netAnalysis_dot

suppressPackageStartupMessages({
  library(CellChat)
  library(NMF)
  library(ggalluvial)
})

celltype_colors <- c("SSC"="#417505", "Spermatocyte"="#4A90E2", "Spermatid"="#D0021B",
                     "Sertoli"="#F5A623", "Leydig"="#880FF3", "Myoid"="#7ED321",
                     "Fibroblast"="#5A5A5A", "Endothelial"="#F16B54")

make_dotplot <- function(rds_path, out_dir, label) {
  cat(">>>", label, ":", rds_path, "\n")
  if (!file.exists(rds_path)) { cat("  SKIP - introuvable\n"); return(invisible(NULL)) }
  cellchat <- readRDS(rds_path)

  # NMF patterns
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 6)

  # Bubble plot custom: X = cell types, Y = pathways, taille = outgoing strength
  prob <- cellchat@netP$prob                # [source, target, pathway]
  out_strength <- apply(prob, c(1, 3), sum) # [source, pathway]
  df <- as.data.frame.table(out_strength, responseName = "strength")
  colnames(df)[1:2] <- c("celltype", "pathway")
  df <- df[df$strength > 0, ]
  df$celltype <- factor(df$celltype,
                        levels = names(celltype_colors)[names(celltype_colors) %in% df$celltype])
  p_dot <- ggplot2::ggplot(df, ggplot2::aes(x = celltype, y = pathway,
                                            size = strength, color = celltype)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = celltype_colors) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, size = "Outgoing\nstrength",
                  title = paste0("Outgoing signaling - ", label)) +
    ggplot2::guides(color = "none")
  n_path <- length(unique(df$pathway))
  ggsave(file.path(out_dir, "bubble_pathway_celltype.png"),
         p_dot, width = 7, height = max(5, 0.22 * n_path + 2),
         dpi = 150, limitsize = FALSE)
  cat("  Bubble:", file.path(out_dir, "bubble_pathway_celltype.png"), "\n")

  # River plot - palette complète: cell types (nos couleurs) + patterns (couleurs neutres)
  pattern_colors <- setNames(
    scales::hue_pal()(6),
    paste0("Pattern ", 1:6)
  )
  river_colors <- c(celltype_colors, pattern_colors)
  p_river <- netAnalysis_river(cellchat, slot.name = "netP", pattern = "outgoing", cutoff = 0.5,
                               color.use = river_colors)
  ggsave(file.path(out_dir, "river_patterns_outgoing.png"),
         p_river, width = 10, height = 8, dpi = 150)
  cat("  River:", file.path(out_dir, "river_patterns_outgoing.png"), "\n\n")
}

base <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results"

make_dotplot(file.path(base, "descended/cellchat_descended.rds"),
             file.path(base, "descended/pathway"), "Descended")
make_dotplot(file.path(base, "crypto/cellchat_crypto.rds"),
             file.path(base, "crypto/pathway"), "Crypto")
make_dotplot(file.path(base, "immuno/cellchat_immuno.rds"),
             file.path(base, "immuno/pathway"), "Immuno")

cat("=== DONE ===\n")
