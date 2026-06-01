#!/usr/bin/env Rscript
# Vérifie la couverture des gènes CellChatDB dans l'objet Seurat

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# 1. Charger l'objet Seurat (original ou humanisé)
# ---------------------------------------------------------------------------
SEURAT_PATH <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object_HUMAN.rds"

message(">>> Chargement du Seurat...")
seurat_obj <- readRDS(SEURAT_PATH)
seurat_genes <- rownames(seurat_obj)
message(">>> ", length(seurat_genes), " gènes dans l'objet Seurat")

# ---------------------------------------------------------------------------
# 2. Charger CellChatDB.human et extraire les gènes
# ---------------------------------------------------------------------------
message(">>> Chargement de CellChatDB.human...")
CellChatDB <- CellChatDB.human

# Extraire tous les gènes utilisés dans la DB
genes_db <- unique(c(
  CellChatDB$interaction$ligand,
  CellChatDB$interaction$receptor,
  CellChatDB$interaction$ agonist,    # si existe
  CellChatDB$interaction$ antagonist, # si existe
  CellChatDB$interaction$ co_A,       # co-agonist
  CellChatDB$interaction$ co_I      # co-inhibitor
))

# Nettoyer : séparer les complexes (gènes séparés par : ou +)
genes_db_split <- unlist(strsplit(genes_db, "[:+]_"))
genes_db_split <- unique(trimws(genes_db_split))
genes_db_split <- genes_db_split[genes_db_split != "" & !is.na(genes_db_split)]

message(">>> ", length(genes_db_split), " gènes uniques dans CellChatDB")

# ---------------------------------------------------------------------------
# 3. Comparer
# ---------------------------------------------------------------------------
present <- genes_db_split[genes_db_split %in% seurat_genes]
missing <- genes_db_split[!genes_db_split %in% seurat_genes]

message("\n=== RÉSULTATS ===")
message("Gènes CellChatDB présents dans le Seurat : ", length(present), " / ", length(genes_db_split))
message("Gènes CellChatDB MANQUANTS          : ", length(missing), " / ", length(genes_db_split))
message("Proportion présente                 : ", round(100 * length(present) / length(genes_db_split), 2), "%")

# ---------------------------------------------------------------------------
# 4. Export
# ---------------------------------------------------------------------------
out_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis"

write.csv(data.frame(gene = present, status = "PRESENT"),
          file.path(out_dir, "cellchatdb_genes_present.csv"), row.names = FALSE)
write.csv(data.frame(gene = missing, status = "MISSING"),
          file.path(out_dir, "cellchatdb_genes_missing.csv"), row.names = FALSE)

message("\nExports :")
message("  - cellchatdb_genes_present.csv")
message("  - cellchatdb_genes_missing.csv")

# Top 20 manquants
message("\n=== Top 20 gènes manquants ===")
print(head(missing, 20))
