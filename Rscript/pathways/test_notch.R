#!/usr/bin/env Rscript
# Test script pour verifier la presence du pathway NOTCH dans l'analyse Descended

library(CellChat)

# Charger l'objet CellChat
file_path <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/descended/cellchat_descended.rds"
cat(">>> Chargement de l'objet CellChat :", file_path, "\n")
cellchat <- readRDS(file_path)

# 1. Liste des pathways
cat("\n=== 1. PATHWAYS DETECTES ===\n")
pathways <- cellchat@netP$pathways
cat("Nombre total de pathways:", length(pathways), "\n")
cat("Liste complete:\n")
print(pathways)

# 2. Chercher NOTCH
cat("\n=== 2. RECHERCHE DE NOTCH ===\n")
notch_found <- "NOTCH" %in% pathways
cat("NOTCH est-il dans les pathways ?", notch_found, "\n")

# Chercher aussi des variantes (case-insensitive, partial match)
notch_matches <- pathways[grep("notch", pathways, ignore.case = TRUE)]
cat("Pathways contenant 'notch':\n")
print(notch_matches)

# 3. Test de contribution pour NOTCH
if (notch_found) {
  cat("\n=== 3. GENERATION DE LA FIGURE NOTCH ===\n")
  out_png <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/descended/pathway/contribution_NOTCH_only.png"
  png(out_png, width = 1000, height = 600)
  netAnalysis_contribution(cellchat, signaling = "NOTCH",
                           title = "Contribution LR pairs - NOTCH pathway")
  dev.off()
  cat("Figure exportee:", out_png, "\n")
} else {
  cat("\n=== 3. NOTCH NON TROUVE - PAS DE FIGURE GENEREE ===\n")
}

# 4. Comparaison avec all_communications.csv
cat("\n=== 4. VERIFICATION DANS all_communications.csv ===\n")
csv_file <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/descended/pathway/all_communications.csv"
if (file.exists(csv_file)) {
  all_comms <- read.csv(csv_file, stringsAsFactors = FALSE)
  if ("pathway_name" %in% names(all_comms)) {
    unique_pathways_csv <- unique(all_comms$pathway_name)
    cat("Pathways uniques dans all_communications.csv:", length(unique_pathways_csv), "\n")
    notch_in_csv <- "NOTCH" %in% unique_pathways_csv
    cat("NOTCH est-il dans le CSV ?", notch_in_csv, "\n")
    if (notch_in_csv) {
      notch_rows <- all_comms[all_comms$pathway_name == "NOTCH", ]
      cat("Nombre de communications NOTCH dans le CSV:", nrow(notch_rows), "\n")
      cat("Premiers couples LR NOTCH:\n")
      print(head(unique(notch_rows$interaction_name_2), 10))
    }
  } else {
    cat("Colonne 'pathway_name' non trouvee dans le CSV. Colonnes disponibles:\n")
    print(names(all_comms))
  }
} else {
  cat("Fichier CSV non trouve:", csv_file, "\n")
}

cat("\n=== DONE ===\n")
