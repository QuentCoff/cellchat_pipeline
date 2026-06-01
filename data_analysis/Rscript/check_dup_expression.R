#!/usr/bin/env Rscript
# Vérifie si les gènes cheval bien nommés sont plus exprimés
# que les gènes LOC uncharacterized dans les doublons humains

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# 1. Charger données
# ---------------------------------------------------------------------------
SEURAT_PATH <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds"
MAPPING_CSV <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/ortholog_mapping_with_horse_symbol.csv"

message(">>> Chargement du Seurat...")
seurat_obj <- readRDS(SEURAT_PATH)
mat <- GetAssayData(seurat_obj, layer = "counts")

message(">>> Chargement du mapping...")
mapping <- read.csv(MAPPING_CSV, stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# 2. Filtrer les doublons humains (status OK)
# ---------------------------------------------------------------------------
mapping_ok <- mapping %>%
  filter(status == "OK", !is.na(human_symbol), human_symbol != "") %>%
  mutate(
    horse_symbol = trimws(horse_symbol),
    human_symbol = trimws(human_symbol)
  )

# Compter les occurrences par human_symbol
counts_human <- table(mapping_ok$human_symbol)
dup_human <- names(counts_human)[counts_human > 1]

message(">>> ", length(dup_human), " symboles humains en double")

# ---------------------------------------------------------------------------
# 3. Pour chaque doublon, comparer l'expression
# ---------------------------------------------------------------------------
results <- list()

for (hsym in dup_human) {
  horse_genes <- mapping_ok$horse_symbol[mapping_ok$human_symbol == hsym]
  
  # Garder ceux présents dans le Seurat
  horse_genes_present <- horse_genes[horse_genes %in% rownames(mat)]
  
  if (length(horse_genes_present) == 0) next
  
  # Somme des counts par gène (expression totale dans tout l'objet)
  expr <- Matrix::rowSums(mat[horse_genes_present, , drop = FALSE])
  
  # Identifier le gène "bien nommé" : celui qui n'est pas LOC
  non_loc <- horse_genes_present[!grepl("^LOC", horse_genes_present)]
  
  # S'il y a plusieurs non-LOC, prendre celui qui match exactement le human_symbol
  if (length(non_loc) > 1) {
    exact_match <- non_loc[toupper(non_loc) == toupper(hsym)]
    if (length(exact_match) > 0) {
      non_loc <- exact_match[1]
    } else {
      non_loc <- non_loc[1]  # fallback
    }
  }
  
  best_gene <- names(which.max(expr))
  best_expr <- max(expr)
  
  is_non_loc_best <- ifelse(length(non_loc) > 0, best_gene == non_loc, NA)
  
  results[[hsym]] <- data.frame(
    human_symbol = hsym,
    n_horse_genes = length(horse_genes),
    n_present = length(horse_genes_present),
    best_gene = best_gene,
    best_expr = best_expr,
    non_loc_gene = ifelse(length(non_loc) > 0, non_loc, NA),
    non_loc_expr = ifelse(length(non_loc) > 0, expr[non_loc], NA),
    non_loc_is_best = is_non_loc_best,
    all_genes = paste(paste0(horse_genes_present, "=", expr[horse_genes_present]), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

results_df <- do.call(rbind, results)
rownames(results_df) <- NULL

# ---------------------------------------------------------------------------
# 4. Résumé
# ---------------------------------------------------------------------------
message("\n=== RÉSULTATS ===")
message("Total doublons analysés : ", nrow(results_df))

# Cas où il y a un non-LOC
has_nonloc <- results_df[!is.na(results_df$non_loc_gene), ]
message("Doublons avec un gène non-LOC : ", nrow(has_nonloc))
message("  -> non-LOC est le plus exprimé : ", sum(has_nonloc$non_loc_is_best, na.rm = TRUE))
message("  -> non-LOC N'EST PAS le plus exprimé : ", sum(!has_nonloc$non_loc_is_best, na.rm = TRUE))

# Exporter
out_csv <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/doublons_expression_check.csv"
write.csv(results_df, out_csv, row.names = FALSE)
message("\nExport : ", out_csv)

# Afficher les cas où non-LOC n'est pas le meilleur
bad <- has_nonloc[!has_nonloc$non_loc_is_best, ]
if (nrow(bad) > 0) {
  message("\n=== Cas où le non-LOC n'est PAS le plus exprimé ===")
  print(bad[, c("human_symbol", "best_gene", "non_loc_gene", "best_expr", "non_loc_expr", "all_genes")])
}
