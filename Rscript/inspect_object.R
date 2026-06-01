#!/usr/bin/env Rscript
# Script pour inspecter le contenu de Final_Annotated_Object.rds et exporter en CSV

library(Seurat)

# Chargement du fichier
obj <- readRDS('/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds')

# Dossier de sortie
out_dir <- '/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/inspect_output'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Classe et dimensions
cat("=== CLASSE ET DIMENSIONS ===\n")
print(class(obj))
print(dim(obj))
info <- data.frame(
  property = c("class", "n_cells", "n_features"),
  value = c(paste(class(obj), collapse = "; "), ncol(obj), nrow(obj))
)
write.csv(info, file.path(out_dir, "object_info.csv"), row.names = FALSE)

# Slots disponibles
cat("\n=== SLOTS DISPONIBLES ===\n")
slots <- slotNames(obj)
print(slots)
write.csv(data.frame(slots = slots), file.path(out_dir, "slots.csv"), row.names = FALSE)

# Assays
cat("\n=== ASSAYS ===\n")
assays <- names(obj@assays)
print(assays)
write.csv(data.frame(assays = assays), file.path(out_dir, "assays.csv"), row.names = FALSE)

# Meta data - colonnes
cat("\n=== META DATA (colonnes) ===\n")
meta_cols <- colnames(obj@meta.data)
print(meta_cols)
write.csv(data.frame(metadata_columns = meta_cols), file.path(out_dir, "metadata_columns.csv"), row.names = FALSE)

# Meta data complete
cat("\n=== EXPORT META DATA COMPLETE ===\n")
write.csv(obj@meta.data, file.path(out_dir, "metadata_all_cells.csv"), row.names = TRUE)

# Cell types (Idents)
cat("\n=== CELL TYPES (Idents) ===\n")
idents_table <- table(Idents(obj))
print(idents_table)
idents_df <- as.data.frame(idents_table)
colnames(idents_df) <- c("cell_type", "count")
write.csv(idents_df, file.path(out_dir, "cell_types_idents.csv"), row.names = FALSE)

# Samples
if ('sample' %in% colnames(obj@meta.data)) {
  cat("\n=== SAMPLES ===\n")
  sample_table <- table(obj@meta.data$sample)
  print(sample_table)
  sample_df <- as.data.frame(sample_table)
  colnames(sample_df) <- c("sample", "count")
  write.csv(sample_df, file.path(out_dir, "samples.csv"), row.names = FALSE)
}

# Conditions
if ('condition' %in% colnames(obj@meta.data)) {
  cat("\n=== CONDITIONS ===\n")
  cond_table <- table(obj@meta.data$condition)
  print(cond_table)
  cond_df <- as.data.frame(cond_table)
  colnames(cond_df) <- c("condition", "count")
  write.csv(cond_df, file.path(out_dir, "conditions.csv"), row.names = FALSE)
}

# Cell types (metadata)
if ('cell_type' %in% colnames(obj@meta.data)) {
  cat("\n=== CELL TYPES (cell_type) ===\n")
  ct_table <- table(obj@meta.data$cell_type)
  print(ct_table)
  ct_df <- as.data.frame(ct_table)
  colnames(ct_df) <- c("cell_type", "count")
  write.csv(ct_df, file.path(out_dir, "cell_types_metadata.csv"), row.names = FALSE)
}

# Crosstab sample x cell_type
cat("\n=== CROSSTAB SAMPLE x CELL TYPE ===\n")
if ('sample' %in% colnames(obj@meta.data) && 'cell_type' %in% colnames(obj@meta.data)) {
  crosstab <- table(obj@meta.data$sample, obj@meta.data$cell_type)
  print(crosstab)
  crosstab_df <- as.data.frame.matrix(crosstab)
  crosstab_df$sample <- rownames(crosstab_df)
  crosstab_df <- crosstab_df[, c("sample", setdiff(colnames(crosstab_df), "sample"))]
  write.csv(crosstab_df, file.path(out_dir, "sample_celltype_crosstab.csv"), row.names = FALSE)
}

cat(paste("\n=== FICHIERS CSV EXPORTES DANS", out_dir, "===\n"))
