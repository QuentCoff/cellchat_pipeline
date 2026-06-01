#!/usr/bin/env Rscript
# Script pour vérifier les informations géniques dans l'objet Seurat

library(Seurat)

cat("=== Chargement de l'objet Seurat ===\n")
seurat_obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds")

cat("\n=== STRUCTURE DE L'OBJET ===\n")
print(seurat_obj)

cat("\n=== ASSAYS DISPONIBLES ===\n")
print(Assays(seurat_obj))

cat("\n=== FEATURES (GÈNES) - 10 premiers ===\n")
features <- rownames(seurat_obj)
print(head(features, 20))
cat(paste0("\nTotal gènes: ", length(features), "\n"))

cat("\n=== VÉRIFICATION: Séquences géniques ? ===\n")
cat("Dans un objet Seurat standard, on a:\n")
cat("- rownames = IDs/noms des gènes (ex: ENS Equus caballus ou noms de gènes)\n")
cat("- Mais PAS les séquences nucléotidiques/protéiques\n")
cat("- Les séquences sont dans les bases de données externes (Ensembl, NCBI...)\n")

cat("\n=== INFOS DANS LES META.FEATURES ===\n")
for (assay in Assays(seurat_obj)) {
  cat(paste0("\nAssay: ", assay, "\n"))
  meta_features <- seurat_obj[[assay]]@meta.features
  if (ncol(meta_features) > 0) {
    cat("Colonnes disponibles:\n")
    print(colnames(meta_features))
    cat("Premières lignes:\n")
    print(head(meta_features, 3))
  } else {
    cat("Aucune meta.feature\n")
  }
}

cat("\n=== EXEMPLE DE GÈNES ===\n")
cat("Format des noms de gènes:\n")
print(head(features, 10))

# Vérifier si ce sont des IDs Ensembl ou des noms de gènes
if (any(grepl("^ENSECAG", features))) {
  cat("\n=> Ce sont des IDs Ensembl (ENSECAG... = Equus caballus)\n")
  cat("Pour avoir les séquences, il faut interroger Ensembl ou BioMart\n")
} else if (any(grepl("^[A-Z][a-z]{2}[0-9]{2,}$", features))) {
  cat("\n=> Ce sont des noms de gènes standard (ex: LAMA1, ITGB1...)\n")
}

cat("\n=== SLots disponibles ===\n")
print(slotNames(seurat_obj))
