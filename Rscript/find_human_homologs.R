#!/usr/bin/env Rscript
# Script pour trouver les homologues humains des genes de gene_export.csv

# Lecture des fichiers
genes_export <- read.csv('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/gene_export.csv', 
                          header = FALSE, stringsAsFactors = FALSE)
colnames(genes_export)[1] <- "gene_symbol"

homology <- read.csv('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/homology.txt', 
                      stringsAsFactors = FALSE)

horse_ids <- read.csv('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/horse_full_ensembl_id.txt', 
                       stringsAsFactors = FALSE)

# Normalisation des noms de colonnes
colnames(horse_ids) <- c("Gene.stable.ID", "Gene.name", "Gene.Synonym")

# Etape 1: Trouver les Ensembl IDs des genes de gene_export.csv
# On cherche par Gene.name ou Gene.Synonym
matched_genes <- data.frame(
  gene_symbol = genes_export$gene_symbol,
  Gene.stable.ID = NA,
  stringsAsFactors = FALSE
)

for (i in 1:nrow(matched_genes)) {
  symbol <- matched_genes$gene_symbol[i]
  
  # Chercher par Gene.name
  match_name <- horse_ids$Gene.stable.ID[horse_ids$Gene.name == symbol]
  if (length(match_name) > 0) {
    matched_genes$Gene.stable.ID[i] <- match_name[1]  # Prendre le premier match
    next
  }
  
  # Chercher par Gene.Synonym
  match_syn <- horse_ids$Gene.stable.ID[horse_ids$Gene.Synonym == symbol]
  if (length(match_syn) > 0) {
    matched_genes$Gene.stable.ID[i] <- match_syn[1]
  }
}

# Statistiques
cat("=== MAPPING GENE SYMBOL -> ENSEMBL ID ===\n")
found <- sum(!is.na(matched_genes$Gene.stable.ID))
cat(paste0("Genes trouves avec Ensembl ID: ", found, "/", nrow(matched_genes), "\n"))
cat(paste0("Genes non trouves: ", sum(is.na(matched_genes$Gene.stable.ID)), "\n"))

# Genes non trouves
not_found <- matched_genes$gene_symbol[is.na(matched_genes$Gene.stable.ID)]
if (length(not_found) > 0) {
  cat("\nGenes non trouves:\n")
  print(not_found)
}

# Etape 2: Trouver les homologues humains
matched_genes$Human.gene.stable.ID <- NA
matched_genes$Human.gene.name <- NA
matched_genes$Orthology.confidence <- NA

for (i in 1:nrow(matched_genes)) {
  if (!is.na(matched_genes$Gene.stable.ID[i])) {
    # Chercher dans homology
    hom_match <- homology[homology$Gene.stable.ID == matched_genes$Gene.stable.ID[i], ]
    if (nrow(hom_match) > 0) {
      # Prendre le premier match avec la meilleure confiance
      best_match <- hom_match[order(hom_match$Human.orthology.confidence..0.low..1.high., decreasing = TRUE)[1], ]
      matched_genes$Human.gene.stable.ID[i] <- best_match$Human.gene.stable.ID
      matched_genes$Human.gene.name[i] <- best_match$Human.gene.name
      matched_genes$Orthology.confidence[i] <- best_match$Human.orthology.confidence..0.low..1.high.
    }
  }
}

# Statistiques homologie
found_human <- sum(!is.na(matched_genes$Human.gene.stable.ID))
cat(paste0("\n=== MAPPING ENSEMBL ID -> HOMOLOGUE HUMAIN ===\n"))
cat(paste0("Genes avec homologue humain: ", found_human, "/", found, "\n"))

# Export
out_dir <- '/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(matched_genes, 
          file.path(out_dir, "gene_export_with_human_homologs.csv"), 
          row.names = FALSE, na = "")

# Resume par confiance d'orthologie
cat("\n=== CONFIANCE DE L'ORTHOLOGIE ===\n")
conf_table <- table(matched_genes$Orthology.confidence, useNA = "ifany")
print(conf_table)

# Genes avec forte confiance (1)
high_conf <- matched_genes[!is.na(matched_genes$Orthology.confidence) & matched_genes$Orthology.confidence == 1, ]
cat(paste0("\nGenes avec haute confiance d'orthologie: ", nrow(high_conf), "\n"))

write.csv(high_conf, 
          file.path(out_dir, "gene_export_high_confidence_homologs.csv"), 
          row.names = FALSE)

cat(paste0("\n=== FICHIERS EXPORTES ===\n"))
cat("- gene_export_with_human_homologs.csv (tous les genes avec mapping)\n")
cat("- gene_export_high_confidence_homologs.csv (genes avec confiance = 1)\n")
