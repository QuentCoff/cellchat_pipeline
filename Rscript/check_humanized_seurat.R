#!/usr/bin/env Rscript
# Script pour inspecter Horse_Humanized_Seurat.rds

library(Seurat)

# Chargement
obj <- readRDS('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds')

cat("=== CLASSE ET DIMENSIONS ===\n")
print(class(obj))
print(dim(obj))

cat("\n=== ASSAYS ===\n")
print(names(obj@assays))

cat("\n=== META DATA ===\n")
print(colnames(obj@meta.data))
print(head(obj@meta.data[, 1:min(5, ncol(obj@meta.data))]))

cat("\n=== CELL TYPES ===\n")
if ('cell_type' %in% colnames(obj@meta.data)) {
  print(table(obj@meta.data$cell_type))
} else if ('CellType' %in% colnames(obj@meta.data)) {
  print(table(obj@meta.data$CellType))
} else {
  print(table(Idents(obj)))
}

cat("\n=== SAMPLES ===\n")
if ('sample' %in% colnames(obj@meta.data)) {
  print(table(obj@meta.data$sample))
}

cat("\n=== CONDITIONS ===\n")
if ('condition' %in% colnames(obj@meta.data)) {
  print(table(obj@meta.data$condition))
}

# Verification des genes (humanises ?)
cat("\n=== APERCU DES GENES (premiers rownames) ===\n")
print(head(rownames(obj), 20))

cat("\n=== NOMS DES GENES ===\n")
# Verifier si les noms commencent par ENSG (human Ensembl)
human_ensg <- sum(grepl("^ENSG", rownames(obj)))
horse_enscag <- sum(grepl("^ENSECAG", rownames(obj)))
symbols <- sum(!grepl("^ENS", rownames(obj)))

cat(paste0("Genes avec prefixe ENSG (humain): ", human_ensg, "\n"))
cat(paste0("Genes avec prefixe ENSECAG (cheval): ", horse_enscag, "\n"))
cat(paste0("Autres (symboles): ", symbols, "\n"))
