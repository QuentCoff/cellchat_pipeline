#!/usr/bin/env Rscript
# Convertit un objet Seurat cheval -> humain via mapping d'orthologues NCBI
# L'objet original reste intact

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

# ---------------------------------------------------------------------------
# 1. Chemins
# ---------------------------------------------------------------------------
SEURAT_ORIGINAL <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds"
MAPPING_CSV    <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/ortholog_mapping_with_horse_symbol.csv"
OUTPUT_RDS     <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object_HUMAN.rds"

# ---------------------------------------------------------------------------
# 2. Charger l'objet original
# ---------------------------------------------------------------------------
message(">>> Chargement de l'objet Seurat original...")
seurat_horse <- readRDS(SEURAT_ORIGINAL)

message(">>> Objet original : ", ncol(seurat_horse), " cellules, ", nrow(seurat_horse), " gènes")

# ---------------------------------------------------------------------------
# 3. Lire le mapping
# ---------------------------------------------------------------------------
message(">>> Lecture du mapping orthologues...")
mapping <- read.csv(MAPPING_CSV, stringsAsFactors = FALSE)

# Filtrer : garder uniquement les orthologues trouvés (status == "OK")
mapping_ok <- mapping %>%
  filter(status == "OK", !is.na(human_symbol), human_symbol != "") %>%
  mutate(
    human_gene_id = as.integer(human_gene_id),
    human_symbol  = trimws(human_symbol)
  ) %>%
  select(horse_symbol, human_symbol)

message(">>> ", nrow(mapping_ok), " gènes avec orthologue humain")

# ---------------------------------------------------------------------------
# 4. Intersection mapping / Seurat
# ---------------------------------------------------------------------------
genes_in_seurat <- rownames(seurat_horse)

mapping_present <- mapping_ok %>% filter(horse_symbol %in% genes_in_seurat)

missing <- setdiff(mapping_ok$horse_symbol, genes_in_seurat)
if (length(missing) > 0) {
  message(">>> ", length(missing), " gènes du mapping absents de l'objet Seurat (ignorés)")
}

message(">>> ", nrow(mapping_present), " gènes cheval présents dans l'objet Seurat avec orthologue")

# ---------------------------------------------------------------------------
# 5. Préparer les nouveaux noms : tous les gènes, orthologués ou non
# ---------------------------------------------------------------------------
new_names <- genes_in_seurat
names(new_names) <- genes_in_seurat

# Renommer les gènes avec orthologue (premier mapping si doublon cheval)
mapping_first <- mapping_present %>%
  group_by(horse_symbol) %>%
  slice(1) %>%
  ungroup()

new_names[mapping_first$horse_symbol] <- mapping_first$human_symbol

# ---------------------------------------------------------------------------
# 6. Extraire, renommer et agréger toute la matrice
# ---------------------------------------------------------------------------
message(">>> Extraction de la matrice de comptage (counts)...")
mat <- GetAssayData(seurat_horse, layer = "counts")

# Renommer toutes les lignes
rownames(mat) <- new_names[rownames(mat)]

# Agréger les lignes avec le même nom (doublons humains + éventuels conflits)
all_names <- rownames(mat)
unique_names <- unique(all_names)

if (length(unique_names) < length(all_names)) {
  message(">>> Agrégation des lignes dupliquées (sparse)...")
  groups <- factor(all_names, levels = unique_names)
  agg_mat <- sparseMatrix(
    i = as.integer(groups),
    j = seq_along(groups),
    x = rep(1, length(groups)),
    dims = c(nlevels(groups), length(groups)),
    dimnames = list(levels(groups), NULL)
  )
  mat_human <- agg_mat %*% mat
} else {
  mat_human <- mat
}

message(">>> Matrice finale : ", nrow(mat_human), " gènes, ", ncol(mat_human), " cellules")

# ---------------------------------------------------------------------------
# 7. Créer le nouvel objet Seurat
# ---------------------------------------------------------------------------
message(">>> Création du nouvel objet Seurat humanisé...")
seurat_human <- CreateSeuratObject(
  counts    = mat_human,
  project   = seurat_horse@project.name,
  meta.data = seurat_horse@meta.data
)

# ---------------------------------------------------------------------------
# 8. Normalisation (sans recalcul PCA/UMAP)
# ---------------------------------------------------------------------------
message(">>> Normalisation...")
seurat_human <- NormalizeData(seurat_human)

message(">>> FindVariableFeatures...")
seurat_human <- FindVariableFeatures(seurat_human, selection.method = "vst", nfeatures = 2000)

message(">>> Scaling...")
seurat_human <- ScaleData(seurat_human)

# ---------------------------------------------------------------------------
# 9. Transférer les identités et réductions de l'original
# ---------------------------------------------------------------------------
Idents(seurat_human) <- Idents(seurat_horse)

# Copier les réductions (PCA, UMAP, etc.) de l'objet original
for (reduc_name in names(seurat_horse@reductions)) {
  seurat_human@reductions[[reduc_name]] <- seurat_horse@reductions[[reduc_name]]
  message(">>> Reduction transférée : ", reduc_name)
}

# Copier aussi les colonnes de clustering si elles existent dans meta.data
for (col in c("seurat_clusters", "cell_type", "condition", "orig.ident")) {
  if (col %in% colnames(seurat_horse@meta.data) && !(col %in% colnames(seurat_human@meta.data))) {
    seurat_human@meta.data[[col]] <- seurat_horse@meta.data[[col]]
  }
}

# ---------------------------------------------------------------------------
# 10. Sauvegarder
# ---------------------------------------------------------------------------
message(">>> Sauvegarde dans : ", OUTPUT_RDS)
saveRDS(seurat_human, file = OUTPUT_RDS)

message(">>> DONE !")
message("    Objet original    : ", ncol(seurat_horse),  " cellules, ", nrow(seurat_horse),  " gènes")
message("    Objet final       : ", ncol(seurat_human), " cellules, ", nrow(seurat_human), " gènes")
message("    Gènes renommés    : ", nrow(mapping_first))
message("    Gènes conservés    : ", nrow(seurat_horse) - nrow(mapping_first))
