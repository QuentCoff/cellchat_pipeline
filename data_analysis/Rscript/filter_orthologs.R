#!/usr/bin/env Rscript

# Script pour filtrer les genes avec orthologues
# 1. Matcher les genes de gene_export.txt avec la colonne initial_alias
# 2. Filtrer ceux qui ont des orthologues (ortholog_name et ortholog_ensg != "N/A")

library(ggplot2)

# Chemin des fichiers
input_file <- "gProfiler_ecaballus_hsapiens_16-04-2026_09-57-05.csv"
gene_list_file <- "../gene_export.txt"
output_file <- "gProfiler_ecaballus_hsapiens_16-04-2026_09-57-05_filtered_orthologs.csv"

# Lecture du fichier CSV
data <- read.csv(input_file, stringsAsFactors = FALSE)

# DEDOUBLONNAGE: garder uniquement la premiere occurrence de chaque gene
data <- data[!duplicated(data$initial_alias), ]
cat("Nombre de genes uniques dans gProfiler (apres dedoublonnage):", nrow(data), "\n")

# Lecture de la liste des genes a matcher
gene_list <- readLines(gene_list_file)
gene_list <- trimws(gene_list)
gene_list <- gene_list[gene_list != ""]  # Supprimer les lignes vides

# Filtrage ETAPES 1: ne garder que les genes presents dans gene_export.txt
data_matched <- data[data$initial_alias %in% gene_list, ]

# Genes supprimes = genes du CSV qui ne sont PAS dans gene_export.txt
data_removed <- data[!(data$initial_alias %in% gene_list), ]
unique_removed_genes <- unique(data_removed$initial_alias)

# Calcul des genes
total_csv_genes <- length(unique(data$initial_alias))
num_unique_matched <- length(unique(data_matched$initial_alias))
num_unique_removed <- length(unique_removed_genes)

# Affichage des statistiques de matching
cat("\n=== ETAPE 1: FILTRAGE PAR gene_export.txt ===\n")
cat("Nombre total de genes uniques dans gProfiler CSV:", total_csv_genes, "\n")
cat("Nombre de genes dans gene_export.txt:", length(gene_list), "\n")
cat("Genes CONSERVES (dans gene_export.txt):", num_unique_matched, "\n")
cat("Genes SUPPRIMES (pas dans gene_export.txt):", num_unique_removed, "\n")

# Export de la liste des genes supprimes (du CSV)
if (num_unique_removed > 0) {
  writeLines(unique_removed_genes, "genes_supprimes_non_dans_gene_export.txt")
  cat("Liste des genes supprimes exportee: genes_supprimes_non_dans_gene_export.txt\n")
}

# Genes de gene_export.txt NON TROUVES dans le CSV
genes_from_list_not_in_csv <- setdiff(gene_list, data_matched$initial_alias)
if (length(genes_from_list_not_in_csv) > 0) {
  writeLines(genes_from_list_not_in_csv, "genes_non_trouves_dans_gProfiler.txt")
  cat("Genes de gene_export.txt non trouves dans gProfiler:", length(genes_from_list_not_in_csv), "\n")
  cat("Liste exportee: genes_non_trouves_dans_gProfiler.txt\n")
}

# Filtrage des genes avec orthologues (ortholog_name et ortholog_ensg != "N/A")
filtered_data <- data_matched[data_matched$ortholog_name != "N/A" & data_matched$ortholog_ensg != "N/A", ]

# Calcul des nombres pour les genes filtres
total_matched <- nrow(data_matched)
genes_with_ortholog <- nrow(filtered_data)
genes_without_ortholog <- total_matched - genes_with_ortholog

# Calcul des pourcentages
pct_with <- round(genes_with_ortholog / total_matched * 100, 2)
pct_without <- round(genes_without_ortholog / total_matched * 100, 2)

# Affichage des resultats
cat("\n=== ETAPE 2: FILTRAGE DES ORTHOLOGUES ===\n")
cat("Nombre total de genes matche's:", total_matched, "\n")
cat("Genes avec orthologues:", genes_with_ortholog, "(", pct_with, "%)\n")
cat("Genes sans orthologues:", genes_without_ortholog, "(", pct_without, "%)\n")

# Export des genes sans orthologues pour analyse biomaRt
genes_sans_ortho <- data_matched[data_matched$ortholog_name == "N/A" | data_matched$ortholog_ensg == "N/A", "initial_alias"]
writeLines(genes_sans_ortho, "genes_sans_orthologues.txt")
cat("\nListe des genes sans orthologues exportee: genes_sans_orthologues.txt\n")

# Ecriture du fichier filtre
write.csv(filtered_data, output_file, row.names = FALSE)
cat("\nFichier filtre sauvegarde:", output_file, "\n")

# Creation du camembert
df_pie <- data.frame(
  category = c(paste0("Avec orthologues\n", genes_with_ortholog, " genes (", pct_with, "%)"),
               paste0("Sans orthologues\n", genes_without_ortholog, " genes (", pct_without, "%)")),
  count = c(genes_with_ortholog, genes_without_ortholog),
  stringsAsFactors = FALSE
)

# Generation du graphique
colors <- c("#2E7D32", "#C62828")  # Vert et rouge

p <- ggplot(df_pie, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = colors) +
  labs(title = "Distribution des genes avec/sans orthologues (genes matche's)",
       subtitle = paste0("Total: ", total_matched, " genes"),
       fill = NULL) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12))

# Sauvegarde du PNG
ggsave("orthologs_distribution.png", p, width = 8, height = 6, dpi = 300)
cat("Graphique sauvegarde: orthologs_distribution.png\n")
