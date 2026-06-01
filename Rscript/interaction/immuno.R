#!/usr/bin/env Rscript
# CellChat analysis - Immuno samples only

library(CellChat)
library(ggplot2)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

# ============================================
# 1. CHARGEMENT ET FILTRAGE DES DONNEES
# ============================================

cat("=== Chargement des données ===\n")
seurat_obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object_HUMAN.rds")

cat(paste0("Dimensions totales: ", ncol(seurat_obj), " cellules\n"))
cat("Conditions disponibles:\n")
print(table(seurat_obj$condition))
cat("\nDetailed groups disponibles:\n")
print(table(seurat_obj$detailed_group))

# Filtrer sur les Immuno uniquement
cat("\n=== Filtrage sur Immuno ===\n")
cells_immuno <- colnames(seurat_obj)[seurat_obj$detailed_group == "Immuno"]
seurat_immuno <- subset(seurat_obj, cells = cells_immuno)

cat(paste0("Cellules Immuno: ", ncol(seurat_immuno), "\n"))
cat("Cell types dans Immuno:\n")
print(table(seurat_immuno$cell_type))

# ============================================
# 2. CREATION DE L'OBJET CELLCHAT (approche manuelle)
# ============================================

cat("\n=== Extraction des données pour CellChat ===\n")

# Extraire les données d'expression normalisées
data.input <- GetAssayData(seurat_immuno, assay = "RNA", layer = "data")
cat(paste0("Dimensions données: ", nrow(data.input), " gènes x ", ncol(data.input), " cellules\n"))

# Récupérer les labels des cellules (cell_type)
labels <- seurat_immuno$cell_type
cat(paste0("Nombre de groupes cellulaires: ", length(unique(labels)), "\n"))

# Créer le dataframe de métadonnées
meta <- data.frame(group = labels, row.names = names(labels))

# Créer l'objet CellChat manuellement
cat("\n=== Création CellChat ===\n")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

celltype_colors <- c("SSC"="#417505", "Spermatocyte"="#4A90E2", "Spermatid"="#D0021B",
                     "Sertoli"="#F5A623", "Leydig"="#880FF3", "Myoid"="#7ED321",
                     "Fibroblast"="#5A5A5A", "Endothelial"="#F16B54")

cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
cat("Groupes: ")
cat(paste(levels(cellchat@idents), collapse = ", "))
cat("\n")

# ============================================
# 3. BASE DE DONNEES LIGAND-RECEPTEUR
# ============================================

cat("\n=== Configuration CellChatDB ===\n")
CellChatDB <- CellChatDB.human  # Données déjà humanisées
cellchat@DB <- CellChatDB

# ============================================
# 4. PRE-TRAITEMENT
# ============================================

cat("\n=== Pré-traitement ===\n")
cellchat <- subsetData(cellchat)  # Subset expression data
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)  # do.fast=FALSE car pas de presto
cellchat <- identifyOverExpressedInteractions(cellchat)

# ============================================
# 5. INFERENCE DU RESEAU DE COMMUNICATION
# ============================================

cat("\n=== Inférence du réseau ===\n")
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Inférence au niveau des pathways
cellchat <- computeCommunProbPathway(cellchat)

# Réseau agrégé
cellchat <- aggregateNet(cellchat)

# ============================================
# 6. VISUALISATIONS
# ============================================

cat("\n=== Génération des visualisations ===\n")
out_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/immuno/interactions"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Réseau agrégé
png(file.path(out_dir, "network_count_immuno.png"), width = 1000, height = 800)
netVisual_circle(cellchat@net$count, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions",
                 edge.width.max = 10, color.use = celltype_colors)
dev.off()

png(file.path(out_dir, "network_weight_immuno.png"), width = 1000, height = 800)
netVisual_circle(cellchat@net$weight, vertex.weight = as.numeric(table(cellchat@idents)),
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength",
                 edge.width.max = 10, color.use = celltype_colors)
dev.off()


png(file.path(out_dir, "network_weight_separated_immuno.png"), width = 1500, height = 1200)
mat <- cellchat@net$count
par(mfrow = c(2,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = as.numeric(table(cellchat@idents)), weight.scale = T,
                   edge.weight.max = max(mat), edge.width.max = 12, title.name = rownames(mat)[i],
                   color.use = celltype_colors)
}
dev.off()

# ============================================
# 6b. ANALYSE DES PATHWAYS ET CONTRIBUTION LR
# ============================================

cat("\n=== Analyse des pathways ===\n")
pathway_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/immuno/pathway"
dir.create(pathway_dir, showWarnings = FALSE, recursive = TRUE)

# Sauvegarder la liste des pathways
pathways <- cellchat@netP$pathways
write.csv(data.frame(pathway = pathways), file = file.path(pathway_dir, "pathways_list.csv"), row.names = FALSE)
cat(paste0("Pathways détectés: ", length(pathways), "\n"))
cat(paste0("Liste exportée: ", file.path(pathway_dir, "pathways_list.csv"), "\n"))

# Exporter les informations sur les pathways depuis CellChatDB
cat("\n=== Extraction des informations des pathways (CellChatDB) ===\n")

# Extraire les infos depuis la base de données interaction
cellchatDB <- CellChatDB$interaction
pathway_annotations <- list()
all_communications <- list()

for (pathway in pathways) {
  # Chercher les interactions pour ce pathway
  pathway_rows <- cellchatDB[cellchatDB$pathway_name == pathway, ]
  
  if (nrow(pathway_rows) > 0) {
    # Prendre la première interaction pour les infos générales
    first_interaction <- pathway_rows[1, ]
    
    # Extraire les valeurs (gérer les NA)
    evidence_val <- ifelse(is.na(first_interaction$evidence) || first_interaction$evidence == "", "Not available", as.character(first_interaction$evidence))
    annotation_val <- ifelse(is.na(first_interaction$annotation) || first_interaction$annotation == "", "Not available", as.character(first_interaction$annotation))
    interaction_val <- ifelse(is.na(first_interaction$interaction_name) || first_interaction$interaction_name == "", "Not available", as.character(first_interaction$interaction_name))
    
    # Compter les couples LR significatifs détectés dans l'analyse
    lr_detected <- subsetCommunication(cellchat, signaling = pathway)
    
    # Compter les couples LR uniques (pas toutes les communications source-target)
    if (!is.null(lr_detected) && nrow(lr_detected) > 0) {
      lr_count_detected <- length(unique(lr_detected$interaction_name))
      # Ajouter à la liste des communications globales
      all_communications[[pathway]] <- lr_detected
    } else {
      lr_count_detected <- 0
    }
    
    pathway_annotations[[pathway]] <- data.frame(
      pathway_name = pathway,
      category = annotation_val,
      evidence = evidence_val,
      lr_count_CellChatDB = nrow(pathway_rows),
      lr_count_detected = lr_count_detected,
      lr_pair_example = interaction_val,
      stringsAsFactors = FALSE
    )
  } else {
    pathway_annotations[[pathway]] <- data.frame(
      pathway_name = pathway,
      category = "Not found in CellChatDB",
      evidence = "Not available",
      lr_count_CellChatDB = 0,
      lr_count_detected = 0,
      lr_pair_example = "Not available",
      stringsAsFactors = FALSE
    )
  }
}

# Sauvegarder toutes les communications dans un seul fichier
if (length(all_communications) > 0) {
  all_comm_df <- do.call(rbind, all_communications)
  write.csv(all_comm_df, file = file.path(pathway_dir, "all_communications.csv"), row.names = FALSE)
  cat(paste0("Toutes les communications exportées: ", file.path(pathway_dir, "all_communications.csv"), "\n"))
}

# Combiner en un seul dataframe
pathways_detected_info <- do.call(rbind, pathway_annotations)
write.csv(pathways_detected_info, file = file.path(pathway_dir, "pathways_info_CellChatDB.csv"), row.names = FALSE)
cat(paste0("Informations CellChatDB exportées: ", file.path(pathway_dir, "pathways_info_CellChatDB.csv"), "\n"))

# Extraire les couples LR enrichis pour tous les pathways
cat("\n=== Extraction des couples LR enrichis ===\n")
lr_enriched <- extractEnrichedLR(cellchat, signaling = pathways, geneLR.return = FALSE)

# Créer une liste des LR par pathway et exporter
lr_list <- list()
for (pathway in pathways) {
  lr_data <- extractEnrichedLR(cellchat, signaling = pathway, geneLR.return = TRUE)
  lr_list[[pathway]] <- lr_data
}

# Sauvegarder en format RDS pour conserver la structure
saveRDS(lr_list, file = file.path(pathway_dir, "LR_enriched_by_pathway.rds"))
cat(paste0("Couples LR enrichis exportés (RDS): ", file.path(pathway_dir, "LR_enriched_by_pathway.rds"), "\n"))

# Créer aussi une version texte résumée
sink(file.path(pathway_dir, "LR_enriched_summary.txt"))
cat("=== COUPLES LIGAND-RECEPTEUR ENRICHIS PAR PATHWAY ===\n\n")
for (pathway in names(lr_list)) {
  cat(paste0("Pathway: ", pathway, "\n"))
  cat(paste0("  Nombre de couples LR: ", length(lr_list[[pathway]]$interaction_name), "\n"))
  cat(paste0("  Couples: ", paste(head(lr_list[[pathway]]$interaction_name, 5), collapse = ", ")))
  if (length(lr_list[[pathway]]$interaction_name) > 5) {
    cat(paste0(" ... et ", length(lr_list[[pathway]]$interaction_name) - 5, " autres"))
  }
  cat("\n\n")
}
sink()
cat(paste0("Résumé LR exporté (TXT): ", file.path(pathway_dir, "LR_enriched_summary.txt"), "\n"))

# Visualisation de la contribution de chaque LR pair (taille adaptative)
cat("\n=== Génération des visualisations de contribution ===\n")

# Calculer la hauteur dynamiquement (~50 pixels par pathway, minimum 600)
plot_height <- max(50 * length(pathways), 600)
plot_width <- 1200

png(file.path(pathway_dir, "contribution_all_pathways.png"), width = plot_width, height = plot_height)
netAnalysis_contribution(cellchat, signaling = pathways, 
                         title = paste0("Contribution of each LR pairs - All ", length(pathways), " pathways"))
dev.off()
cat(paste0("Visualisation contribution exportée (", length(pathways), " pathways, ", plot_height, "px de haut): ", file.path(pathway_dir, "contribution_all_pathways.png"), "\n"))

# Visualisation pour les top 5 pathways (si plus de 5 pathways)
if (length(pathways) >= 5) {
  png(file.path(pathway_dir, "contribution_top5_pathways.png"), width = 1200, height = 800)
  netAnalysis_contribution(cellchat, signaling = pathways[1:5], 
                           title = paste("Top 5 pathways:", paste(pathways[1:5], collapse = ", ")))
  dev.off()
  cat(paste0("Visualisation top 5 pathways exportée: ", file.path(pathway_dir, "contribution_top5_pathways.png"), "\n"))
}

# ============================================
# 6c. DOTPLOT: PATHWAYS x CELL TYPES (fonction native CellChat)
# ============================================
cat("\n=== Génération du dotplot pathways x cell types (CellChat) ===\n")

suppressPackageStartupMessages(library(NMF))

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 6)

suppressPackageStartupMessages(library(ggalluvial))

# Bubble plot custom: X = cell types, Y = pathways, taille = outgoing strength
prob <- cellchat@netP$prob
out_strength <- apply(prob, c(1, 3), sum)
df_bub <- as.data.frame.table(out_strength, responseName = "strength")
colnames(df_bub)[1:2] <- c("celltype", "pathway")
df_bub <- df_bub[df_bub$strength > 0, ]
df_bub$celltype <- factor(df_bub$celltype,
                          levels = names(celltype_colors)[names(celltype_colors) %in% df_bub$celltype])
p_bubble <- ggplot2::ggplot(df_bub, ggplot2::aes(x = celltype, y = pathway,
                                                 size = strength, color = celltype)) +
  ggplot2::geom_point() +
  ggplot2::scale_color_manual(values = celltype_colors) +
  ggplot2::scale_size_continuous(range = c(1, 8)) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(x = NULL, y = NULL, size = "Outgoing\nstrength",
                title = "Outgoing signaling - Immuno") +
  ggplot2::guides(color = "none")
n_path <- length(unique(df_bub$pathway))
ggsave(file.path(pathway_dir, "bubble_pathway_celltype.png"),
       p_bubble, width = 7, height = max(5, 0.22 * n_path + 2),
       dpi = 150, limitsize = FALSE)
cat(paste0("Bubble plot exporté: ", file.path(pathway_dir, "bubble_pathway_celltype.png"), "\n"))

p_river <- netAnalysis_river(cellchat, slot.name = "netP", pattern = "outgoing", cutoff = 0.5)
ggsave(file.path(pathway_dir, "river_patterns_outgoing.png"),
       p_river, width = 10, height = 8, dpi = 150)
cat(paste0("River plot outgoing exporté: ", file.path(pathway_dir, "river_patterns_outgoing.png"), "\n"))

# ============================================
# 7. EXPORT DES GENES ET COUPLES LIGAND-RECEPTEUR
# ============================================

cat("\n=== Export des gènes utilisés par CellChat ===\n")

# Exporter les couples ligand-récepteur significatifs
lr_pairs <- cellchat@LR$LRsig
write.csv(lr_pairs, file = file.path(out_dir, "LR_couples_significatifs.csv"), row.names = FALSE)
cat(paste0("Couples LR exportés: ", file.path(out_dir, "LR_couples_significatifs.csv"), "\n"))

# Exporter la liste unique des gènes utilisés
genes_ligands <- unique(lr_pairs$ligand)
genes_receptors <- unique(lr_pairs$receptor)
all_genes <- unique(c(genes_ligands, genes_receptors))

genes_df <- data.frame(
  gene = all_genes,
  type = ifelse(all_genes %in% genes_ligands & all_genes %in% genes_receptors, "both",
                ifelse(all_genes %in% genes_ligands, "ligand", "receptor"))
)
write.csv(genes_df, file = file.path(out_dir, "genes_utilises.csv"), row.names = FALSE)
cat(paste0("Gènes utilisés exportés: ", file.path(out_dir, "genes_utilises.csv"), "\n"))
cat(paste0("  - Nombre total de gènes uniques: ", length(all_genes), "\n"))
cat(paste0("  - Ligands: ", length(genes_ligands), "\n"))
cat(paste0("  - Récepteurs: ", length(genes_receptors), "\n"))

# ============================================
# 8. EXPORT DES STATISTIQUES
# ============================================

cat("\n=== Export des statistiques ===\n")

# Créer le fichier de rapport
report_file <- file.path(out_dir, "statistiques_immuno.txt")
sink(report_file)

cat("========================================\n")
cat("RAPPORT STATISTIQUES - IMMUNO\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("========================================\n\n")

# 1. Dimensions totales
cat("1. DIMENSIONS GLOBALES\n")
cat("----------------------\n")
cat(paste0("Dimensions totales de l'objet Seurat: ", ncol(seurat_obj), " cellules\n"))
cat(paste0("Nombre de gènes: ", nrow(seurat_obj), "\n\n"))

cat("Conditions disponibles:\n")
print(table(seurat_obj$condition))
cat("\n")

# 2. Statistiques Immuno
cat("2. STATISTIQUES IMMUNO\n")
cat("----------------------\n")
cat(paste0("Cellules Immuno: ", ncol(seurat_immuno), "\n\n"))

cat("Cell types dans Immuno:\n")
cell_counts <- table(seurat_immuno$cell_type)
print(cell_counts)
cat("\n")

# 3. Gènes par cellule en moyenne par cluster
cat("3. NOMBRE DE GENES DETECTES PAR CELLULE (moyenne par cluster)\n")
cat("-------------------------------------------------------------\n")

# Calculer le nombre de gènes détectés par cellule
seurat_immuno$nGene <- Matrix::colSums(GetAssayData(seurat_immuno, assay = "RNA", layer = "counts") > 0)

# Moyenne par cluster
avg_genes_per_cluster <- tapply(seurat_immuno$nGene, seurat_immuno$cell_type, mean)
cat("Moyenne de gènes détectés par cellule, par cluster:\n")
print(round(avg_genes_per_cluster, 2))
cat("\n")

# 4. Couples récepteur-ligand par cluster
cat("4. COUPLES RECEPTEUR-LIGAND PAR CLUSTER\n")
cat("---------------------------------------\n")
cat("Nombre total de couples LR identifiés: ", length(cellchat@LR$LRsig$interaction_name), "\n")
cat("Fichier des couples LR: ", file.path(out_dir, "LR_couples_significatifs.csv"), "\n")
cat("Fichier des gènes utilisés: ", file.path(out_dir, "genes_utilises.csv"), "\n\n")

# Nombre d'interactions (count) par cluster (sender)
cat("Nombre d'interactions sortantes (count) par cluster:\n")
for (i in 1:nrow(cellchat@net$count)) {
  cluster_name <- rownames(cellchat@net$count)[i]
  total_interactions <- sum(cellchat@net$count[i, ])
  cat(paste0("  ", cluster_name, ": ", total_interactions, " interactions\n"))
}
cat("\n")

cat("Résumé des interactions:\n")
cat(paste0("  - Nombre total d'interactions: ", sum(cellchat@net$count), "\n"))
cat(paste0("  - Nombre moyen d'interactions par cluster: ", round(mean(rowSums(cellchat@net$count)), 2), "\n"))
cat("\n")

# 5. Résumé CellChat
cat("5. RESUME CELLCHAT\n")
cat("------------------\n")
cat(paste0("Cell groups: ", length(levels(cellchat@idents)), "\n"))
cat(paste0("Groupes: ", paste(levels(cellchat@idents), collapse = ", "), "\n\n"))

# Fermer le fichier
sink()
cat(paste0("Statistiques exportées dans: ", report_file, "\n"))

# Sauvegarde de l'objet CellChat
rds_file <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/immuno/cellchat_immuno.rds"
saveRDS(cellchat, file = rds_file)
cat(paste0("Objet CellChat sauvegardé dans: ", rds_file, "\n"))