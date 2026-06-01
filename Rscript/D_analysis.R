# ==============================================================================
# ANALYSE CELLCHAT — FOCUS DESCENDED
# Comparaisons individuelles entre échantillons Descended
# Sains : D-CV19, D-CV25, D-CV30, D-CV42
# Uni-Descended : D-CV22, D-CV44
# ==============================================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. CHARGEMENT
# ==============================================================================
print("--- Chargement de l'objet Seurat humanisé ---")
obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds")
Idents(obj) <- "cell_type"

# Subset Descended uniquement
obj_desc <- subset(obj, subset = condition == "Descended")

# Classification des échantillons
sains <- c("D-CV19", "D-CV25", "D-CV30", "D-CV42")
uni   <- c("D-CV22", "D-CV44")

cat("Échantillons Descended disponibles:\n")
print(table(obj_desc$Horse))
cat("\nSains:", paste(sains, collapse=", "), "\n")
cat("Uni-Descended:", paste(uni, collapse=", "), "\n")

# Dossier de sortie
outdir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/descended_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

CellChatDB <- CellChatDB.human

# ==============================================================================
# 2. PIPELINE CELLCHAT PAR ÉCHANTILLON
# ==============================================================================

run_cellchat_sample <- function(seurat_obj, horse_id) {
  cat("\n=== CellChat pour:", horse_id, "===\n")
  
  obj_sub <- subset(seurat_obj, subset = Horse == horse_id)
  cat("  Cellules:", ncol(obj_sub), "\n")
  cat("  Types cellulaires:", paste(sort(unique(as.character(obj_sub$cell_type))), collapse=", "), "\n")
  
  # Vérifier qu'il y a assez de types cellulaires
  if (length(unique(obj_sub$cell_type)) < 2) {
    cat("  SKIP: pas assez de types cellulaires\n")
    return(NULL)
  }
  
  cc <- createCellChat(object = obj_sub, group.by = "cell_type", assay = "RNA")
  cc@DB <- CellChatDB
  
  cc <- subsetData(cc)
  cat("  Identification des gènes surexprimés...\n")
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  cc <- identifyOverExpressedInteractions(cc)
  
  cat("  Calcul des probabilités de communication...\n")
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  
  cat("  Calcul des voies de signalisation...\n")
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  
  cat("  Calcul de la centralité...\n")
  cc <- netAnalysis_computeCentrality(cc)
  
  cat("  Done!\n")
  return(cc)
}

# Lancer pour chaque échantillon Descended
all_horses <- c(sains, uni)
cellchat_list <- list()

for (h in all_horses) {
  cellchat_list[[h]] <- run_cellchat_sample(obj_desc, h)
  if (!is.null(cellchat_list[[h]])) {
    saveRDS(cellchat_list[[h]], file.path(outdir, paste0("cellchat_", h, ".rds")))
  }
}

# Retirer les NULL (échantillons sans assez de cellules)
cellchat_list <- cellchat_list[!sapply(cellchat_list, is.null)]

# ==============================================================================
# 3. VISUALISATIONS INDIVIDUELLES
# ==============================================================================
print("--- Visualisations individuelles par échantillon ---")

for (h in names(cellchat_list)) {
  cc <- cellchat_list[[h]]
  
  # Réseau nb interactions
  png(file.path(outdir, paste0("network_count_", h, ".png")), width=800, height=800)
  netVisual_circle(cc@net$count, vertex.weight = table(cc@idents),
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste("Nb interactions -", h))
  dev.off()
  
  # Réseau force interactions
  png(file.path(outdir, paste0("network_weight_", h, ".png")), width=800, height=800)
  netVisual_circle(cc@net$weight, vertex.weight = table(cc@idents),
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste("Force interactions -", h))
  dev.off()
}

# ==============================================================================
# 4. COMPARAISONS SAINS vs SAINS
# ==============================================================================
print("--- Comparaisons entre échantillons sains ---")

dir_sain_vs_sain <- file.path(outdir, "sains_vs_sains")
dir.create(dir_sain_vs_sain, showWarnings = FALSE)

sains_available <- intersect(sains, names(cellchat_list))

if (length(sains_available) >= 2) {
  combis_sains <- combn(sains_available, 2, simplify = FALSE)
  
  for (pair in combis_sains) {
    h1 <- pair[1]; h2 <- pair[2]
    cat("  Comparaison:", h1, "vs", h2, "\n")
    
    obj_merge <- tryCatch({
      merged <- mergeCellChat(list(cellchat_list[[h1]], cellchat_list[[h2]]),
                              add.names = c(h1, h2))
      merged
    }, error = function(e) { cat("    Erreur merge:", conditionMessage(e), "\n"); NULL })
    
    if (!is.null(obj_merge)) {
      # Barplot nb + force
      png(file.path(dir_sain_vs_sain, paste0("compare_", h1, "_vs_", h2, ".png")), width=1200, height=600)
      gg1 <- compareInteractions(obj_merge, show.legend = FALSE, group = c(1,2))
      gg2 <- compareInteractions(obj_merge, show.legend = FALSE, group = c(1,2), measure = "weight")
      print(gg1 + gg2)
      dev.off()
      
      # Réseau différentiel
      png(file.path(dir_sain_vs_sain, paste0("diff_net_", h1, "_vs_", h2, ".png")), width=1200, height=600)
      par(mfrow = c(1,2))
      netVisual_diffInteraction(obj_merge, weight.scale = TRUE)
      netVisual_diffInteraction(obj_merge, weight.scale = TRUE, measure = "weight")
      dev.off()
      
      # Heatmap différentielle
      png(file.path(dir_sain_vs_sain, paste0("diff_heatmap_", h1, "_vs_", h2, ".png")), width=1200, height=600)
      gg1 <- netVisual_heatmap(obj_merge)
      gg2 <- netVisual_heatmap(obj_merge, measure = "weight")
      print(gg1 + gg2)
      dev.off()
      
      saveRDS(obj_merge, file.path(dir_sain_vs_sain, paste0("merged_", h1, "_vs_", h2, ".rds")))
    }
  }
}

# ==============================================================================
# 5. COMPARAISONS UNI vs CHAQUE SAIN
# ==============================================================================
print("--- Comparaisons Uni-Descended vs Sains ---")

dir_uni_vs_sain <- file.path(outdir, "uni_vs_sains")
dir.create(dir_uni_vs_sain, showWarnings = FALSE)

uni_available <- intersect(uni, names(cellchat_list))

for (u in uni_available) {
  for (s in sains_available) {
    cat("  Comparaison:", u, "(Uni) vs", s, "(Sain)\n")
    
    obj_merge <- tryCatch({
      merged <- mergeCellChat(list(cellchat_list[[u]], cellchat_list[[s]]),
                              add.names = c(paste0(u, "_Uni"), paste0(s, "_Sain")))
      merged
    }, error = function(e) { cat("    Erreur merge:", conditionMessage(e), "\n"); NULL })
    
    if (!is.null(obj_merge)) {
      # Barplot nb + force
      png(file.path(dir_uni_vs_sain, paste0("compare_", u, "_vs_", s, ".png")), width=1200, height=600)
      gg1 <- compareInteractions(obj_merge, show.legend = FALSE, group = c(1,2))
      gg2 <- compareInteractions(obj_merge, show.legend = FALSE, group = c(1,2), measure = "weight")
      print(gg1 + gg2)
      dev.off()
      
      # Réseau différentiel
      png(file.path(dir_uni_vs_sain, paste0("diff_net_", u, "_vs_", s, ".png")), width=1200, height=600)
      par(mfrow = c(1,2))
      netVisual_diffInteraction(obj_merge, weight.scale = TRUE)
      netVisual_diffInteraction(obj_merge, weight.scale = TRUE, measure = "weight")
      dev.off()
      
      # Heatmap différentielle
      png(file.path(dir_uni_vs_sain, paste0("diff_heatmap_", u, "_vs_", s, ".png")), width=1200, height=600)
      gg1 <- netVisual_heatmap(obj_merge)
      gg2 <- netVisual_heatmap(obj_merge, measure = "weight")
      print(gg1 + gg2)
      dev.off()
      
      saveRDS(obj_merge, file.path(dir_uni_vs_sain, paste0("merged_", u, "_vs_", s, ".rds")))
    }
  }
}

# ==============================================================================
# 6. RÉSUMÉ
# ==============================================================================
print("--- Résumé des interactions par échantillon ---")

summary_df <- data.frame(
  Horse = character(),
  Group = character(),
  Nb_Cells = integer(),
  Nb_Interactions = integer(),
  Interaction_Strength = numeric(),
  stringsAsFactors = FALSE
)

for (h in names(cellchat_list)) {
  cc <- cellchat_list[[h]]
  grp <- ifelse(h %in% sains, "Sain", "Uni-Descended")
  summary_df <- rbind(summary_df, data.frame(
    Horse = h,
    Group = grp,
    Nb_Cells = length(cc@idents),
    Nb_Interactions = sum(cc@net$count),
    Interaction_Strength = round(sum(cc@net$weight), 3)
  ))
}

print(summary_df)
write.csv(summary_df, file.path(outdir, "summary_interactions.csv"), row.names = FALSE)

# Barplot résumé
p_summary <- ggplot(summary_df, aes(x = reorder(Horse, -Nb_Interactions), y = Nb_Interactions, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Sain" = "#4E79A7", "Uni-Descended" = "#F28E2B")) +
  theme_minimal() +
  labs(title = "Nombre d'interactions par échantillon Descended",
       x = "Échantillon", y = "Nb interactions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "summary_barplot_interactions.png"), p_summary, width = 8, height = 6)

p_strength <- ggplot(summary_df, aes(x = reorder(Horse, -Interaction_Strength), y = Interaction_Strength, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Sain" = "#4E79A7", "Uni-Descended" = "#F28E2B")) +
  theme_minimal() +
  labs(title = "Force des interactions par échantillon Descended",
       x = "Échantillon", y = "Force des interactions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "summary_barplot_strength.png"), p_strength, width = 8, height = 6)

print("SUCCÈS : Analyse Descended terminée ! Résultats dans /results/descended_analysis/")
