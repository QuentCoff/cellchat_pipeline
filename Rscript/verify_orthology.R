# ==============================================================================
# VÉRIFICATION GLOBALE DE L'ORTHOLOGIE
# ==============================================================================

library(Seurat)
library(biomaRt)

# ==============================================================================
# 1. STATISTIQUES GLOBALES SUR LE MAPPING LOCAL
# ==============================================================================
print("--- Statistiques Globales de l'Orthologie ---")

# Chargement des gènes originaux
obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds")
genes_horse <- rownames(obj)
total_genes <- length(genes_horse)

# Chargement de la table téléchargée
local_csv <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/ortho_complete_horse_human.csv"
all_ortho <- read.csv(local_csv)
colnames(all_ortho) <- c("Horse_Gene", "Human_Gene")

# Nettoyage des lignes vides
all_ortho_clean <- all_ortho[all_ortho$Human_Gene != "" & all_ortho$Horse_Gene != "", ]

# Filtre sur nos gènes
our_ortho <- all_ortho_clean[all_ortho_clean$Horse_Gene %in% genes_horse, ]

# Calculs
genes_with_match <- length(unique(our_ortho$Horse_Gene))
genes_lost <- total_genes - genes_with_match
success_rate <- (genes_with_match / total_genes) * 100
loss_rate <- (genes_lost / total_genes) * 100

# Doublons
dup_horse <- sum(duplicated(our_ortho$Horse_Gene)) # 1 gène cheval -> plusieurs gènes humains
dup_human <- sum(duplicated(our_ortho$Human_Gene)) # Plusieurs gènes cheval -> 1 gène humain

cat(sprintf("Total des gènes de l'objet (Cheval) : %d\n", total_genes))
cat(sprintf("Gènes avec une correspondance Humaine : %d (%.2f%%)\n", genes_with_match, success_rate))
cat(sprintf("Gènes perdus (sans correspondance)  : %d (%.2f%%)\n", genes_lost, loss_rate))
cat("\n--- Analyse des doublons ---\n")
cat(sprintf("Doublons côté Cheval (1 gène cheval -> N humains) : %d\n", dup_horse))
cat(sprintf("Doublons côté Humain (N gènes cheval -> 1 humain) : %d\n", dup_human))

# ==============================================================================
# 2. VÉRIFICATION EN DIRECT (LIVE TEST) SUR L'ENSEMBLE DES GÈNES
# ==============================================================================
print("\n--- Vérification en direct via Ensembl (Tous les gènes) ---")
print("Cela peut prendre un peu de temps (requête par lots)...")

connect_horse <- function() {
  for (m in c("www", "useast", "asia")) {
    cat("  Tentative miroir:", m, "...\n")
    mart <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "ecaballus_gene_ensembl", mirror = m),
      error = function(e) { cat("    Échec:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(mart)) { cat("    OK!\n"); return(mart) }
  }
  for (v in c(113, 112, 111)) {
    cat("  Tentative archive version", v, "...\n")
    mart <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "ecaballus_gene_ensembl", version = v),
      error = function(e) { cat("    Échec:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(mart)) { cat("    OK!\n"); return(mart) }
  }
  stop("Impossible de se connecter à Ensembl")
}

horse <- connect_horse()

# Découper en lots de 500 gènes pour éviter les timeouts
batch_size <- 500
batches <- split(genes_horse, ceiling(seq_along(genes_horse) / batch_size))

live_results_list <- list()

for (i in seq_along(batches)) {
  cat(sprintf("  Requête lot %d/%d...\n", i, length(batches)))
  for(attempt in 1:3) {
    res <- tryCatch({
      getBM(
        attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
        filters = "external_gene_name",
        values = batches[[i]],
        mart = horse
      )
    }, error = function(e) NULL)
    
    if(!is.null(res)) {
      live_results_list[[i]] <- res
      break
    } else {
      cat(sprintf("    Échec lot %d (tentative %d). Pause de 5s...\n", i, attempt))
      Sys.sleep(5)
    }
  }
}

# Fusionner tous les résultats
live_results <- do.call(rbind, live_results_list)
colnames(live_results) <- c("Horse_Gene", "Ensembl_Human_Gene")

# Nettoyage
live_results <- live_results[live_results$Ensembl_Human_Gene != "", ]
live_results <- live_results[!duplicated(live_results$Horse_Gene), ]

comparison <- merge(live_results, our_ortho[!duplicated(our_ortho$Horse_Gene), c("Horse_Gene", "Human_Gene")], by = "Horse_Gene", all.x = TRUE)
colnames(comparison)[3] <- "Local_Human_Gene"
comparison$Match <- comparison$Ensembl_Human_Gene == comparison$Local_Human_Gene

# ==============================================================================
# 3. GÉNÉRATION DU RAPPORT RÉCAPITULATIF
# ==============================================================================
print("\n--- Génération du rapport Markdown ---")

report_path <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/RAPPORT_ORTHOLOGIE.md"

sink(report_path)
cat("# Rapport d'Orthologie : Cheval -> Humain\n\n")

cat("## 1. Statistiques Globales (Mapping Local)\n")
cat(sprintf("- **Total des gènes de l'objet Seurat (Cheval)** : %d\n", total_genes))
cat(sprintf("- **Gènes avec une correspondance Humaine** : %d (%.2f%%)\n", genes_with_match, success_rate))
cat(sprintf("- **Gènes perdus (sans correspondance)** : %d (%.2f%%)\n", genes_lost, loss_rate))

cat("\n### Détails des Doublons\n")
cat(sprintf("- **Doublons côté Cheval** : %d (Un même gène cheval pointe vers plusieurs gènes humains. Seul le premier est conservé lors du mapping.)\n", dup_horse))
cat(sprintf("- **Doublons côté Humain** : %d (Plusieurs gènes de cheval différents convergent vers le même gène humain. Cela peut entraîner une fusion des comptes de ces gènes lors de la création de l'objet Seurat.)\n", dup_human))

cat("\n## 2. Vérification en direct (Ensembl via biomaRt)\n")
cat("*Un test a été effectué en interrogeant directement les serveurs d'Ensembl aujourd'hui pour vérifier si la base locale est toujours à jour.*\n\n")

cat(sprintf("- **Total gènes cherchés dans Seurat** : %d\n", total_genes))
cat(sprintf("- **Total gènes trouvés sur Ensembl aujourd'hui** : %d\n", nrow(comparison)))
cat(sprintf("- **Correspondances parfaites avec le fichier local** : %d (%.2f%% des gènes trouvés)\n", 
            sum(comparison$Match, na.rm = TRUE),
            (sum(comparison$Match, na.rm = TRUE) / nrow(comparison)) * 100))

if (all(comparison$Match, na.rm = TRUE)) {
  cat("\n**Conclusion** : ✅ Le mapping local est **100% identique** aux données d'Ensembl actuelles.\n")
} else {
  diffs <- sum(!comparison$Match, na.rm = TRUE)
  cat(sprintf("\n**Conclusion** : ⚠️ **%d différences** ont été trouvées entre Ensembl aujourd'hui et le fichier local.\n", diffs))
  cat("\n### Aperçu des différences :\n")
  cat("```\n")
  print(head(comparison[!comparison$Match, ]))
  cat("```\n")
}
sink()

cat(sprintf("Rapport généré avec succès : %s\n", report_path))

