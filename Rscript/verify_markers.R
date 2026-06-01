# ==============================================================================
# VÉRIFICATION DE LA PRÉSENCE DES GÈNES MARQUEURS DANS L'OBJET HUMANISÉ
# ==============================================================================

library(Seurat)

# 1. Liste des marqueurs d'intérêt
marqueurs <- c("DAZL", "DMRT1", "GFRA1", "SETX", "HORMAD2", "PIWIL1", 
               "CCDC168", "EFCAB3", "CREM", "ERBB4", "SYNE2", "FSHR", "SOX9", 
               "LHCGR", "STAR", "ACTA2", "MYH11", "DCN", "TCF21", "EBF1", 
               "VWF", "ETS1", "FLI1")

# 2. Chargement de l'objet humanisé
print("--- Chargement de l'objet Seurat Humanisé ---")
obj_path <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds"
obj <- readRDS(obj_path)

# 3. Extraction des gènes disponibles
genes_dispos <- rownames(obj)

# 4. Comparaison
marqueurs_trouves <- marqueurs[marqueurs %in% genes_dispos]
marqueurs_perdus <- setdiff(marqueurs, genes_dispos)

# 5. Affichage des résultats
cat("\n=================================================================\n")
cat(sprintf("RÉSULTATS DE LA VÉRIFICATION DES MARQUEURS (%d testés)\n", length(marqueurs)))
cat("=================================================================\n\n")

cat(sprintf("✅ Marqueurs PRÉSENTS dans l'objet (%d) :\n", length(marqueurs_trouves)))
if(length(marqueurs_trouves) > 0) {
  print(marqueurs_trouves)
} else {
  cat("Aucun marqueur trouvé.\n")
}

cat(sprintf("\n❌ Marqueurs ABSENTS de l'objet (%d) :\n", length(marqueurs_perdus)))
if(length(marqueurs_perdus) > 0) {
  print(marqueurs_perdus)
} else {
  cat("Tous les marqueurs sont présents !\n")
}

cat("\n=================================================================\n")
