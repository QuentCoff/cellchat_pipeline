#!/usr/bin/env Rscript
# Breakdown des gènes CellChatDB manquants : pas d'orthologue vs pas exprimé

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
})

SEURAT_PATH <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object_HUMAN.rds"
MAPPING_CSV <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/ortholog_mapping_with_horse_symbol.csv"

# ---------------------------------------------------------------------------
# 1. Charger
# ---------------------------------------------------------------------------
seurat_obj <- readRDS(SEURAT_PATH)
seurat_genes <- rownames(seurat_obj)

mapping <- read.csv(MAPPING_CSV, stringsAsFactors = FALSE)
mapping_ok <- mapping %>%
  filter(status == "OK", !is.na(human_symbol), human_symbol != "") %>%
  mutate(human_symbol = trimws(human_symbol),
         horse_symbol = trimws(horse_symbol))

# Gènes humains avec orthologue cheval dans le mapping
human_with_ortho <- unique(mapping_ok$human_symbol)

# ---------------------------------------------------------------------------
# 2. CellChatDB genes
# ---------------------------------------------------------------------------
CellChatDB <- CellChatDB.human
genes_db <- unique(c(
  CellChatDB$interaction$ligand,
  CellChatDB$interaction$receptor,
  CellChatDB$interaction$co_A,
  CellChatDB$interaction$co_I
))
genes_db_split <- unlist(strsplit(genes_db, "[:+]_"))
genes_db_split <- unique(trimws(genes_db_split))
genes_db_split <- genes_db_split[genes_db_split != "" & !is.na(genes_db_split)]

# ---------------------------------------------------------------------------
# 3. Classification des manquants
# ---------------------------------------------------------------------------
missing <- genes_db_split[!genes_db_split %in% seurat_genes]

# 3a. Pas d'orthologue dans le mapping
no_ortho <- missing[!missing %in% human_with_ortho]

# 3b. Orthologue existe mais pas dans le Seurat (donc pas exprimé)
with_ortho_but_missing <- missing[missing %in% human_with_ortho]

# 3c. Pour les avec orthologue : l'orthologue cheval est-il dans le Seurat ?
horse_counterparts <- mapping_ok %>%
  filter(human_symbol %in% with_ortho_but_missing) %>%
  group_by(human_symbol) %>%
  summarise(
    horse_genes = paste(unique(horse_symbol), collapse = "; "),
    in_seurat = any(horse_symbol %in% seurat_genes),
    .groups = 'drop'
  )

with_ortho_in_seurat <- horse_counterparts %>% filter(in_seurat == TRUE)
with_ortho_not_in_seurat <- horse_counterparts %>% filter(in_seurat == FALSE)

# ---------------------------------------------------------------------------
# 4. Résumé
# ---------------------------------------------------------------------------
message("=== BREAKDOWN DES GENES CELLCHATDB MANQUANTS ===")
message("Total gènes CellChatDB : ", length(genes_db_split))
message("Manquants dans le Seurat : ", length(missing))
message("")
message("1. PAS d'orthologue cheval identifié : ", length(no_ortho), " gènes")
message("   -> Manque dans le CSV de mapping")
message("")
message("2. Orthologue cheval identifié, mais horse_symbol PAS dans le Seurat : ", nrow(with_ortho_not_in_seurat), " gènes")
message("   -> Orthologue existe mais probablement pas exprimé dans tes données")
message("")
message("3. Orthologue cheval identifié ET horse_symbol présent dans le Seurat : ", nrow(with_ortho_in_seurat), " gènes")
message("   -> Devrait être là ! Nom différent ou filtré ?")
message("")

if (nrow(with_ortho_in_seurat) > 0) {
  message("=== Cas étranges (orthologue cheval présent dans Seurat mais human_symbol manquant) ===")
  print(head(with_ortho_in_seurat, 20))
}

# ---------------------------------------------------------------------------
# 5. Export
# ---------------------------------------------------------------------------
out_dir <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis"

write.csv(data.frame(gene = no_ortho, reason = "NO_ORTHOLOG_IN_MAPPING"),
          file.path(out_dir, "cellchatdb_missing_no_ortholog.csv"), row.names = FALSE)

write.csv(with_ortho_not_in_seurat,
          file.path(out_dir, "cellchatdb_missing_ortholog_not_expressed.csv"), row.names = FALSE)

write.csv(with_ortho_in_seurat,
          file.path(out_dir, "cellchatdb_missing_ortholog_present.csv"), row.names = FALSE)

message("\nExports:")
message("  - cellchatdb_missing_no_ortholog.csv")
message("  - cellchatdb_missing_ortholog_not_expressed.csv")
message("  - cellchatdb_missing_ortholog_present.csv")
