# ==============================================================================
# SCRIPT D'ORTHOLOGIE : CHEVAL -> HUMAIN (Version Miroir Résilient)
# ==============================================================================

# 1. Charger les bibliothèques
library(biomaRt)
library(Seurat)

# 2. CHARGER L'OBJET
print("--- Chargement de l'objet Seurat original ---")
obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Final_Annotated_Object.rds")

# 3. CHARGEMENT DE LA TABLE D'ORTHOLOGIE (déjà téléchargée)
print("--- Chargement de la table d'orthologie locale ---")
genes_horse <- rownames(obj)

ortho_csv <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/ortho_complete_horse_human.csv"
all_ortho <- read.csv(ortho_csv)
colnames(all_ortho) <- c("Horse_Gene", "Human_Gene")
cat("  Table chargée:", nrow(all_ortho), "lignes\n")

# Mapping local avec nos gènes
ortho_table <- all_ortho[all_ortho$Horse_Gene %in% genes_horse, ]

cat("  Orthologues trouvés pour nos gènes:", nrow(ortho_table), "\n")
no_ortho <- setdiff(genes_horse, ortho_table$Horse_Gene)
cat("  Gènes cheval sans orthologue humain:", length(no_ortho), "\n")
print("--- Aperçu de la table de conversion : ---")
print(head(ortho_table))

# 5. NETTOYAGE ET GESTION DES DOUBLONS
print("--- Nettoyage de la table d'orthologie ---")

# Supprimer les gènes vides
ortho_table <- ortho_table[ortho_table$Human_Gene != "", ]

# Supprimer les doublons (crucial pour le renommage de la matrice)
ortho_table <- ortho_table[!duplicated(ortho_table$Horse_Gene), ]
ortho_table <- ortho_table[!duplicated(ortho_table$Human_Gene), ]

# Sauvegarde locale du mapping pour archive
write.csv(ortho_table, "ortho_mapping_final.csv", row.names = FALSE)

# 6. CRÉATION DE L'OBJET SEURAT HUMANISÉ
print("--- Conversion de l'objet Seurat vers noms Humains ---")

# On ne garde que les gènes présents dans notre table de mapping
obj_sub <- subset(obj, features = ortho_table$Horse_Gene)

# Extraction des données
counts <- GetAssayData(obj_sub, assay = "RNA", layer = "counts")
data <- GetAssayData(obj_sub, assay = "RNA", layer = "data")

# Création du vecteur de mapping
new_names <- ortho_table$Human_Gene
names(new_names) <- ortho_table$Horse_Gene

# Renommage des lignes
rownames(counts) <- new_names[rownames(counts)]
rownames(data) <- new_names[rownames(data)]

# Reconstruction de l'objet Seurat "Humain"
obj_human <- CreateSeuratObject(counts = counts, meta.data = obj_sub@meta.data)
obj_human <- SetAssayData(obj_human, layer = "data", new.data = data)

# Transfert des Idents (types cellulaires)
cat("Colonnes de métadonnées disponibles:\n")
print(colnames(obj_sub@meta.data))

# Chercher automatiquement la colonne de type cellulaire
cell_col <- intersect(c("cell_type", "celltype", "CellType", "cell.type", 
                        "seurat_clusters", "cluster", "ident"), 
                      colnames(obj_sub@meta.data))
if (length(cell_col) > 0) {
  cat("Colonne utilisée pour Idents:", cell_col[1], "\n")
  Idents(obj_human) <- obj_sub@meta.data[[cell_col[1]]]
} else {
  cat("ATTENTION: Aucune colonne de type cellulaire trouvée, Idents non transférés.\n")
  cat("Tu pourras les assigner manuellement après.\n")
}

# 7. VÉRIFICATION ET SAUVEGARDE FINALE
print("Vérification des 5 premiers gènes humanisés :")
print(head(rownames(obj_human)))

print("--- Sauvegarde de l'objet FINAL (Format RDS) ---")
saveRDS(obj_human, "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds")

print("SUCCÈS : Ton objet humanisé est prêt et sauvegardé !")