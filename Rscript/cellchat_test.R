# ==============================================================================
# CELLCHAT : ANALYSE DES COMMUNICATIONS INTERCELLULAIRES
# Testicule de cheval (gènes humanisés) — Descended vs Crypto vs Immuno
# ==============================================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# 1. CHARGEMENT DE L'OBJET HUMANISÉ
# ==============================================================================
print("--- Chargement de l'objet Seurat humanisé ---")
obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds")

# ==============================================================================
# 2. VÉRIFICATION DES MÉTADONNÉES (déjà présentes dans l'objet annoté)
# ==============================================================================
print("--- Vérification des métadonnées ---")
cat("Colonnes disponibles:\n")
print(colnames(obj@meta.data))

Idents(obj) <- "cell_type"

cat("\nDistribution des types cellulaires:\n")
print(table(obj$cell_type))
cat("\nDistribution par condition:\n")
print(table(obj$condition))

# ==============================================================================
# 3. PIPELINE CELLCHAT PAR CONDITION
# ==============================================================================
print("--- Lancement du pipeline CellChat ---")

outdir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

CellChatDB <- CellChatDB.human

run_cellchat <- function(seurat_obj, condition_name) {
  cat("\n=== CellChat pour:", condition_name, "===\n")
  
  # Subset par condition
  obj_sub <- subset(seurat_obj, subset = condition == condition_name)
  cat("  Cellules:", ncol(obj_sub), "\n")
  cat("  Types cellulaires:", paste(unique(obj_sub$cell_type), collapse=", "), "\n")
  
  # Créer l'objet CellChat
  cc <- createCellChat(object = obj_sub, group.by = "cell_type", assay = "RNA")
  cc@DB <- CellChatDB
  
  # Pipeline standard
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

# Lancer pour chaque condition
cellchat.descended <- run_cellchat(obj, "Descended")
cellchat.crypto    <- run_cellchat(obj, "Crypto")
cellchat.immuno    <- run_cellchat(obj, "Immuno")

# ==============================================================================
# 4. SAUVEGARDE DES OBJETS INDIVIDUELS
# ==============================================================================
print("--- Sauvegarde des objets CellChat ---")
saveRDS(cellchat.descended, file.path(outdir, "cellchat_descended.rds"))
saveRDS(cellchat.crypto,    file.path(outdir, "cellchat_crypto.rds"))
saveRDS(cellchat.immuno,    file.path(outdir, "cellchat_immuno.rds"))

# ==============================================================================
# 5. COMPARAISON DESCENDED vs CRYPTO
# ==============================================================================
print("--- Comparaison Descended vs Crypto ---")

object.list <- list(Descended = cellchat.descended, Crypto = cellchat.crypto)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))

# Nombre et force des interactions
png(file.path(outdir, "comparison_nb_interactions.png"), width=1200, height=600)
gg1 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat.merged, show.legend = FALSE, group = c(1,2), measure = "weight")
print(gg1 + gg2)
dev.off()

# Différentiel des interactions par type cellulaire
png(file.path(outdir, "diff_interactions_network.png"), width=1200, height=600)
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE)
netVisual_diffInteraction(cellchat.merged, weight.scale = TRUE, measure = "weight")
dev.off()

# Heatmap du différentiel
png(file.path(outdir, "diff_interactions_heatmap.png"), width=1200, height=600)
gg1 <- netVisual_heatmap(cellchat.merged)
gg2 <- netVisual_heatmap(cellchat.merged, measure = "weight")
print(gg1 + gg2)
dev.off()

# ==============================================================================
# 6. COMPARAISON DESCENDED vs IMMUNO
# ==============================================================================
print("--- Comparaison Descended vs Immuno ---")

object.list2 <- list(Descended = cellchat.descended, Immuno = cellchat.immuno)
cellchat.merged2 <- mergeCellChat(object.list2, add.names = names(object.list2))

png(file.path(outdir, "comparison_desc_vs_immuno.png"), width=1200, height=600)
gg1 <- compareInteractions(cellchat.merged2, show.legend = FALSE, group = c(1,2))
gg2 <- compareInteractions(cellchat.merged2, show.legend = FALSE, group = c(1,2), measure = "weight")
print(gg1 + gg2)
dev.off()

# ==============================================================================
# 7. VISUALISATIONS PAR CONDITION (réseaux individuels)
# ==============================================================================
print("--- Visualisations réseau par condition ---")

for (name in c("descended", "crypto", "immuno")) {
  cc <- get(paste0("cellchat.", name))
  
  png(file.path(outdir, paste0("network_count_", name, ".png")), width=800, height=800)
  netVisual_circle(cc@net$count, vertex.weight = table(cc@idents),
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste("Nb interactions -", name))
  dev.off()
  
  png(file.path(outdir, paste0("network_weight_", name, ".png")), width=800, height=800)
  netVisual_circle(cc@net$weight, vertex.weight = table(cc@idents),
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = paste("Force interactions -", name))
  dev.off()
}

# ==============================================================================
# 8. SAUVEGARDE FINALE
# ==============================================================================
saveRDS(cellchat.merged,  file.path(outdir, "cellchat_merged_desc_vs_crypto.rds"))
saveRDS(cellchat.merged2, file.path(outdir, "cellchat_merged_desc_vs_immuno.rds"))

print("SUCCÈS : Analyse CellChat terminée ! Résultats dans /results/")