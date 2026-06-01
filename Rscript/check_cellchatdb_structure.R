#!/usr/bin/env Rscript
# Script pour explorer la structure de CellChatDB

library(CellChat)

# Charger la base de données
cellchatDB <- CellChatDB.human

cat("=== STRUCTURE DE CellChatDB ===\n\n")

# Voir les slots/noms disponibles
cat("Noms disponibles dans CellChatDB:\n")
print(names(cellchatDB))

cat("\n=== TABLE INTERACTION ===\n")
cat("Dimensions:", dim(cellchatDB$interaction), "\n")
cat("Colonnes:\n")
print(colnames(cellchatDB$interaction))

cat("\n=== PREMIERES LIGNES INTERACTION ===\n")
print(head(cellchatDB$interaction[, 1:min(10, ncol(cellchatDB$interaction))], 3))

cat("\n=== EXEMPLE PATHWAY NOTCH ===\n")
notch_data <- subsetDB(cellchatDB, search = "NOTCH")
if (!is.null(notch_data$interaction)) {
  cat("Colonnes dans interaction:\n")
  print(colnames(notch_data$interaction))
  cat("\nPremière interaction NOTCH:\n")
  print(notch_data$interaction[1, ])
}

cat("\n=== AUTRES TABLES DISPONIBLES ===\n")
for (name in names(cellchatDB)) {
  obj <- cellchatDB[[name]]
  if (is.data.frame(obj) || is.matrix(obj)) {
    cat(paste0(name, ": ", nrow(obj), " rows x ", ncol(obj), " cols\n"))
  } else if (is.list(obj)) {
    cat(paste0(name, ": liste avec ", length(obj), " éléments\n"))
  }
}
