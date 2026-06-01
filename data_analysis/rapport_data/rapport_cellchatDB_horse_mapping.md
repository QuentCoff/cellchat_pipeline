=== BREAKDOWN DES GENES CELLCHATDB MANQUANTS ===
Total gènes CellChatDB : 1479
Manquants dans le Seurat : 505

1. PAS d'orthologue cheval identifié : 505 gènes
   -> Manque dans le CSV de mapping

2. Orthologue cheval identifié, mais horse_symbol PAS dans le Seurat : 0 gènes
   -> Orthologue existe mais probablement pas exprimé dans tes données

3. Orthologue cheval identifié ET horse_symbol présent dans le Seurat : 0 gènes
   -> Devrait être là ! Nom différent ou filtré ?


Exports:
  - cellchatdb_missing_no_ortholog.csv
  - cellchatdb_missing_ortholog_not_expressed.csv
  - cellchatdb_missing_ortholog_present.csv