#!/usr/bin/env Rscript
# Script pour trouver les genes en commun entre homology.txt et horse_full_ensembl_id.txt

# Lecture des fichiers
homology <- read.csv('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/homology.txt', 
                      stringsAsFactors = FALSE)
horse_ids <- read.csv('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/horse_full_ensembl_id.txt', 
                       stringsAsFactors = FALSE)

# Extraction des Gene stable ID uniques
genes_homology <- unique(homology$Gene.stable.ID)
genes_horse <- unique(horse_ids$Gene.stable.ID)

# Genes en commun
common_genes <- intersect(genes_homology, genes_horse)

# Statistiques
cat("=== STATISTIQUES ===\n")
cat(paste0("Genes uniques dans homology.txt: ", length(genes_homology), "\n"))
cat(paste0("Genes uniques dans horse_full_ensembl_id.txt: ", length(genes_horse), "\n"))
cat(paste0("Genes en commun: ", length(common_genes), "\n"))
cat(paste0("Genes uniquement dans homology: ", length(setdiff(genes_homology, genes_horse)), "\n"))
cat(paste0("Genes uniquement dans horse: ", length(setdiff(genes_horse, genes_homology)), "\n"))

# Export des genes en commun
out_dir <- '/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(data.frame(Gene.stable.ID = common_genes), 
          file.path(out_dir, "common_genes_ids.csv"), 
          row.names = FALSE)

# Export avec infos supplementaires des deux fichiers
common_homology <- homology[homology$Gene.stable.ID %in% common_genes, ]
common_horse <- horse_ids[horse_ids$Gene.stable.ID %in% common_genes, ]

write.csv(common_homology, 
          file.path(out_dir, "common_genes_homology_details.csv"), 
          row.names = FALSE)

write.csv(common_horse, 
          file.path(out_dir, "common_genes_horse_details.csv"), 
          row.names = FALSE)

# Liste des genes uniques a chaque fichier
only_homology <- setdiff(genes_homology, genes_horse)
only_horse <- setdiff(genes_horse, genes_homology)

write.csv(data.frame(Gene.stable.ID = only_homology), 
          file.path(out_dir, "only_in_homology.csv"), 
          row.names = FALSE)

write.csv(data.frame(Gene.stable.ID = only_horse), 
          file.path(out_dir, "only_in_horse.csv"), 
          row.names = FALSE)

cat(paste0("\n=== FICHIERS EXPORTES DANS ", out_dir, " ===\n"))
cat("- common_genes_ids.csv (liste des IDs communs)\n")
cat("- common_genes_homology_details.csv (details depuis homology.txt)\n")
cat("- common_genes_horse_details.csv (details depuis horse_full_ensembl_id.txt)\n")
cat("- only_in_homology.csv (genes uniquement dans homology)\n")
cat("- only_in_horse.csv (genes uniquement dans horse_full_ensembl_id.txt)\n")
