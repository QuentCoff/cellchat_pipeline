#!/usr/bin/env Rscript
# ==============================================================================
# SCRIPT: Recherche d'orthologues pour genes perdus via biomaRt
# ==============================================================================
# Ce script prend les genes sans orthologues (de l'analyse gProfiler)
# et fait une requete biomaRt pour trouver d'eventuels orthologues alternatifs

library(biomaRt)

# ==============================================================================
# CONFIGURATION DES FICHIERS
# ==============================================================================
LOST_GENES_FILE <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/genes_sans_orthologues.txt"
OUTPUT_FILE <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/lost_genes_biomart_results.csv"
REPORT_FILE <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/rapport_biomart.txt"

# Initialisation du rapport
report_lines <- c(
  paste(rep("=", 80), collapse = ""),
  "RAPPORT DE RECHERCHE BIOMART - ORTHOLOGUES CHEVAL -> HUMAIN",
  paste(rep("=", 80), collapse = ""),
  paste0("Date d'execution: ", Sys.time()),
  ""
)

# ==============================================================================
# VERIFICATION INITIALE
# ==============================================================================
cat("\n[1/7] VERIFICATION DU FICHIER D'ENTREE...\n")
report_lines <- c(report_lines, "[1/7] VERIFICATION DU FICHIER D'ENTREE", "")

if (!file.exists(LOST_GENES_FILE)) {
  msg <- paste0("ERREUR: Fichier non trouve: ", LOST_GENES_FILE)
  cat(msg, "\n")
  cat("Assure-toi d'avoir lance filter_orthologs.R d'abord.\n")
  stop(msg)
}

cat("  ✓ Fichier trouve:", LOST_GENES_FILE, "\n")
report_lines <- c(report_lines, paste0("  Fichier d'entree: ", LOST_GENES_FILE))

# ==============================================================================
# LECTURE DES GENES PERDUS
# ==============================================================================
cat("\n[2/7] LECTURE DES GENES PERDUS...\n")
report_lines <- c(report_lines, "", "[2/7] LECTURE DES GENES PERDUS", "")

lost_genes <- readLines(LOST_GENES_FILE)
lost_genes <- trimws(lost_genes)
lost_genes <- lost_genes[lost_genes != ""]

nb_genes <- length(lost_genes)
cat("  ✓ Nombre de genes a analyser:", nb_genes, "\n")
report_lines <- c(report_lines, paste0("  Nombre de genes a analyser: ", nb_genes))

# Apercu des premiers genes
if (nb_genes > 0) {
  cat("  ✓ Apercu des 10 premiers genes:\n")
  for (i in 1:min(10, nb_genes)) {
    cat("     -", lost_genes[i], "\n")
  }
  if (nb_genes > 10) {
    cat("     ... et", nb_genes - 10, "autres\n")
  }
}

# ==============================================================================
# CONNEXION A ENSEMBL
# ==============================================================================
cat("\n[3/7] CONNEXION AUX BASES ENSEMBL...\n")
report_lines <- c(report_lines, "", "[3/7] CONNEXION AUX BASES ENSEMBL", "")

cat("  → Tentative de connexion au serveur principal...\n")
connection_success <- FALSE

tryCatch({
  ensembl_horse <- useMart("ensembl", dataset = "ecaballus_gene_ensembl")
  cat("  ✓ Connecte a ecaballus_gene_ensembl (serveur principal)\n")
  report_lines <- c(report_lines, "  Serveur: principal (www.ensembl.org)")
  report_lines <- c(report_lines, "  Dataset cheval: ecaballus_gene_ensembl ✓")
  
  ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  cat("  ✓ Connecte a hsapiens_gene_ensembl (serveur principal)\n")
  report_lines <- c(report_lines, "  Dataset humain: hsapiens_gene_ensembl ✓")
  connection_success <- TRUE
  
}, error = function(e) {
  cat("  ✗ Echec connexion principale:", e$message, "\n")
  cat("  → Tentative avec le serveur mirror (useast)...\n")
  
  tryCatch({
    ensembl_horse <<- useMart("ensembl", dataset = "ecaballus_gene_ensembl", 
                              host = "https://useast.ensembl.org")
    cat("  ✓ Connecte a ecaballus_gene_ensembl (mirror USEAST)\n")
    report_lines <<- c(report_lines, "  Serveur: mirror (useast.ensembl.org)")
    report_lines <<- c(report_lines, "  Dataset cheval: ecaballus_gene_ensembl ✓")
    
    ensembl_human <<- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                              host = "https://useast.ensembl.org")
    cat("  ✓ Connecte a hsapiens_gene_ensembl (mirror USEAST)\n")
    report_lines <<- c(report_lines, "  Dataset humain: hsapiens_gene_ensembl ✓")
    connection_success <<- TRUE
    
  }, error = function(e2) {
    cat("  ✗ Echec total de connexion:", e2$message, "\n")
    report_lines <<- c(report_lines, "  ERREUR: Impossible de se connecter a Ensembl")
  })
})

if (!connection_success) {
  stop("Impossible de se connecter aux serveurs Ensembl")
}

# ==============================================================================
# ETAPE 1: MAPPING GENE SYMBOL -> ENSEMBL ID CHEVAL
# ==============================================================================
cat("\n[4/7] RECHERCHE DES ENSEMBL IDs CHEVAL...\n")
report_lines <- c(report_lines, "", "[4/7] RECHERCHE DES ENSEMBL IDs CHEVAL", "")

cat("  → Requete biomaRt en cours...\n")
cat("  → Attributs demandes: external_gene_name, ensembl_gene_id, gene_biotype\n")

t1 <- Sys.time()
horse_ensembl <- tryCatch({
  getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", "gene_biotype"),
    filters = "external_gene_name",
    values = lost_genes,
    mart = ensembl_horse
  )
}, error = function(e) {
  cat("  ✗ ERREUR requete cheval:", conditionMessage(e), "\n")
  cat("  → Nouvelle tentative avec serveur USEAST...\n")
  
  # Reconnexion au mirror
  ensembl_horse_mirror <- useMart("ensembl", dataset = "ecaballus_gene_ensembl",
                                   host = "https://useast.ensembl.org")
  getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", "gene_biotype"),
    filters = "external_gene_name",
    values = lost_genes,
    mart = ensembl_horse_mirror
  )
})
t2 <- Sys.time()

cat("  ✓ Requete terminee en", round(difftime(t2, t1, units = "secs"), 2), "secondes\n")

found_horse <- length(unique(horse_ensembl$external_gene_name))
not_found <- nb_genes - found_horse

cat("  ✓ Resultats:\n")
cat("     - Genes trouves dans Ensembl cheval:", found_horse, "/", nb_genes, "\n")
cat("     - Genes NON TROUVES:", not_found, "\n")
cat("     - Taux de recuperation:", round(found_horse/nb_genes*100, 2), "%\n")

report_lines <- c(report_lines,
  paste0("  Genes trouves: ", found_horse, "/", nb_genes),
  paste0("  Genes non trouves: ", not_found),
  paste0("  Taux de recuperation: ", round(found_horse/nb_genes*100, 2), "%"),
  paste0("  Temps de requete: ", round(difftime(t2, t1, units = "secs"), 2), " secondes")
)

# Genes non trouves dans Ensembl cheval
not_found_horse <- setdiff(lost_genes, horse_ensembl$external_gene_name)
if (length(not_found_horse) > 0) {
  cat("  → Export des genes non trouves...\n")
  NOT_FOUND_FILE <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_analysis/genes_non_trouves_ensembl_horse.txt"
  writeLines(not_found_horse, NOT_FOUND_FILE)
  cat("  ✓ Liste exportee:", NOT_FOUND_FILE, "\n")
  report_lines <- c(report_lines, paste0("  Fichier genes non trouves: ", NOT_FOUND_FILE))
  
  # Apercu des genes non trouves
  cat("  → Apercu des genes non trouves (max 10):\n")
  for (i in 1:min(10, length(not_found_horse))) {
    cat("     -", not_found_horse[i], "\n")
  }
}

# ==============================================================================
# ETAPE 2: RECHERCHE DES ORTHOLOGUES HUMAINS
# ==============================================================================
cat("\n[5/7] RECHERCHE DES ORTHOLOGUES HUMAINS...\n")
report_lines <- c(report_lines, "", "[5/7] RECHERCHE DES ORTHOLOGUES HUMAINS", "")

horse_ids <- unique(horse_ensembl$ensembl_gene_id)
cat("  →", length(horse_ids), "Ensembl IDs cheval a traiter\n")
cat("  → Attributs demandes:\n")
cat("     - ensembl_gene_id (ID cheval)\n")
cat("     - external_gene_name (nom cheval)\n")
cat("     - hsapiens_homolog_ensembl_gene (ID humain)\n")
cat("     - hsapiens_homolog_associated_gene_name (nom humain)\n")
cat("     - hsapiens_homolog_orthology_type (type d'orthologie)\n")
cat("     - hsapiens_homolog_orthology_confidence (confiance 0-1)\n")
cat("     - hsapiens_homolog_perc_id (% identite proteique)\n")
cat("     - hsapiens_homolog_perc_id_r1 (% identite reciproque)\n")

cat("  → Requete biomaRt en cours...\n")
t1 <- Sys.time()
orthologs <- tryCatch({
  getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "hsapiens_homolog_ensembl_gene",
      "hsapiens_homolog_associated_gene_name",
      "hsapiens_homolog_orthology_type",
      "hsapiens_homolog_orthology_confidence",
      "hsapiens_homolog_perc_id",
      "hsapiens_homolog_perc_id_r1"
    ),
    filters = "ensembl_gene_id",
    values = horse_ids,
    mart = ensembl_horse
  )
}, error = function(e) {
  cat("  ✗ ERREUR requete orthologues:", conditionMessage(e), "\n")
  cat("  → Nouvelle tentative avec serveur USEAST...\n")
  ensembl_horse_mirror <- useMart("ensembl", dataset = "ecaballus_gene_ensembl",
                                   host = "https://useast.ensembl.org")
  getBM(
    attributes = c(
      "ensembl_gene_id",
      "external_gene_name",
      "hsapiens_homolog_ensembl_gene",
      "hsapiens_homolog_associated_gene_name",
      "hsapiens_homolog_orthology_type",
      "hsapiens_homolog_orthology_confidence",
      "hsapiens_homolog_perc_id",
      "hsapiens_homolog_perc_id_r1"
    ),
    filters = "ensembl_gene_id",
    values = horse_ids,
    mart = ensembl_horse_mirror
  )
})
t2 <- Sys.time()

cat("  ✓ Requete terminee en", round(difftime(t2, t1, units = "secs"), 2), "secondes\n")

total_ortho <- nrow(orthologs)
cat("  ✓ Total orthologues bruts trouves:", total_ortho, "\n")

# Filtrer pour garder uniquement les orthologues avec donnees humaines valides
orthologs_valid <- orthologs[
  orthologs$hsapiens_homolog_ensembl_gene != "" & 
  !is.na(orthologs$hsapiens_homolog_ensembl_gene),
]

valid_ortho <- nrow(orthologs_valid)
cat("  ✓ Orthologues avec donnees humaines valides:", valid_ortho, "\n")
cat("  ✓ Orthologues exclus (donnees vides):", total_ortho - valid_ortho, "\n")

report_lines <- c(report_lines,
  paste0("  Ensembl IDs cheval traites: ", length(horse_ids)),
  paste0("  Orthologues bruts trouves: ", total_ortho),
  paste0("  Orthologues valides: ", valid_ortho),
  paste0("  Orthologues exclus: ", total_ortho - valid_ortho),
  paste0("  Temps de requete: ", round(difftime(t2, t1, units = "secs"), 2), " secondes")
)

# ==============================================================================
# ETAPE 3: ENRICHISSEMENT AVEC LES DONNEES HUMAINES
# ==============================================================================
if (nrow(orthologs_valid) > 0) {
  cat("\n[6/7] ENRICHISSEMENT DES DONNEES HUMAINES...\n")
  report_lines <- c(report_lines, "", "[6/7] ENRICHISSEMENT DES DONNEES HUMAINES", "")
  
  human_ids <- unique(orthologs_valid$hsapiens_homolog_ensembl_gene)
  cat("  →", length(human_ids), "genes humains uniques a enrichir\n")
  cat("  → Attributs demandes:\n")
  cat("     - ensembl_gene_id\n")
  cat("     - external_gene_name\n")
  cat("     - description\n")
  cat("     - gene_biotype\n")
  cat("     - chromosome_name\n")
  
  cat("  → Requete biomaRt sur hsapiens...\n")
  t1 <- Sys.time()
  human_info <- tryCatch({
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "description",
        "gene_biotype",
        "chromosome_name"
      ),
      filters = "ensembl_gene_id",
      values = human_ids,
      mart = ensembl_human
    )
  }, error = function(e) {
    cat("  ✗ ERREUR requete humain:", conditionMessage(e), "\n")
    cat("  → Nouvelle tentative avec serveur USEAST...\n")
    ensembl_human_mirror <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                                     host = "https://useast.ensembl.org")
    getBM(
      attributes = c(
        "ensembl_gene_id",
        "external_gene_name",
        "description",
        "gene_biotype",
        "chromosome_name"
      ),
      filters = "ensembl_gene_id",
      values = human_ids,
      mart = ensembl_human_mirror
    )
  })
  t2 <- Sys.time()
  
  cat("  ✓ Requete terminee en", round(difftime(t2, t1, units = "secs"), 2), "secondes\n")
  cat("  ✓ Informations recuperees pour", nrow(human_info), "genes humains\n")
  
  report_lines <- c(report_lines,
    paste0("  Genes humains a enrichir: ", length(human_ids)),
    paste0("  Genes enrichis: ", nrow(human_info)),
    paste0("  Temps de requete: ", round(difftime(t2, t1, units = "secs"), 2), " secondes")
  )
  
  # Fusion des donnees
  results <- merge(
    orthologs_valid,
    human_info,
    by.x = "hsapiens_homolog_ensembl_gene",
    by.y = "ensembl_gene_id",
    all.x = TRUE
  )
  
  # Renommage des colonnes apres merge
  colnames(results) <- c(
    "Human_Ensembl_ID",
    "Horse_Ensembl_ID",
    "Horse_Gene_Symbol",
    "Human_Gene_Symbol_Homolog",
    "Orthology_Type",
    "Confidence",
    "Perc_Identity",
    "Perc_Identity_Reciprocal",
    "Human_Gene_Symbol",
    "Human_Description",
    "Human_Biotype",
    "Human_Chromosome"
  )
  
  # Reordonner les colonnes
  results <- results[, c(
    "Horse_Gene_Symbol",
    "Horse_Ensembl_ID",
    "Human_Ensembl_ID",
    "Human_Gene_Symbol",
    "Human_Gene_Symbol_Homolog",
    "Orthology_Type",
    "Confidence",
    "Perc_Identity",
    "Perc_Identity_Reciprocal",
    "Human_Description",
    "Human_Biotype",
    "Human_Chromosome"
  )]
  
  # Renommage final des colonnes
  colnames(results) <- c(
    "Horse_Gene_Symbol",
    "Horse_Ensembl_ID",
    "Human_Ensembl_ID",
    "Human_Gene_Symbol_BioMart",
    "Human_Gene_Symbol_Homolog",
    "Orthology_Type",
    "Confidence",
    "Perc_Identity",
    "Perc_Identity_Reciprocal",
    "Human_Description",
    "Human_Biotype",
    "Human_Chromosome"
  )
  
  # Export
  write.csv(results, OUTPUT_FILE, row.names = FALSE)
  cat("\n[7/7] EXPORT ET STATISTIQUES FINALES...\n")
  report_lines <- c(report_lines, "", "[7/7] EXPORT ET STATISTIQUES FINALES", "")
  
  cat("  ✓ Fichier CSV exporte:", OUTPUT_FILE, "\n")
  cat("  ✓ Nombre total d'orthologues:", nrow(results), "\n")
  report_lines <- c(report_lines,
    paste0("  Fichier CSV: ", OUTPUT_FILE),
    paste0("  Nombre total d'orthologues: ", nrow(results))
  )
  
  # ==============================================================================
  # STATISTIQUES DETAILLEES
  # ==============================================================================
  cat("\n  --- STATISTIQUES PAR CONFIANCE ---\n")
  conf_stats <- table(results$Confidence, useNA = "ifany")
  cat("  Confiance 0 (faible):", sum(conf_stats[names(conf_stats) == "0"], na.rm = TRUE), "\n")
  cat("  Confiance 1 (haute):", sum(conf_stats[names(conf_stats) == "1"], na.rm = TRUE), "\n")
  cat("  Confiance NA:", sum(is.na(results$Confidence)), "\n")
  
  report_lines <- c(report_lines, "",
    "STATISTIQUES PAR CONFIANCE",
    paste0("  Confiance 0 (faible): ", sum(conf_stats[names(conf_stats) == "0"], na.rm = TRUE)),
    paste0("  Confiance 1 (haute): ", sum(conf_stats[names(conf_stats) == "1"], na.rm = TRUE)),
    paste0("  Confiance NA: ", sum(is.na(results$Confidence)))
  )
  
  # Genes avec haute confiance
  high_conf <- results[results$Confidence == 1, ]
  cat("\n  --- GENES AVEC HAUTE CONFIANCE (1) ---\n")
  cat("  Nombre:", nrow(high_conf), "\n")
  report_lines <- c(report_lines, "",
    "GENES AVEC HAUTE CONFIANCE",
    paste0("  Nombre: ", nrow(high_conf))
  )
  
  # Distribution par type d'orthologie
  cat("\n  --- DISTRIBUTION PAR TYPE D'ORTHOLOGIE ---\n")
  ortho_types <- table(results$Orthology_Type)
  print(ortho_types)
  report_lines <- c(report_lines, "",
    "DISTRIBUTION PAR TYPE D'ORTHOLOGIE"
  )
  for (type in names(ortho_types)) {
    report_lines <- c(report_lines, paste0("  ", type, ": ", ortho_types[type]))
    cat("  ", type, ":", ortho_types[type], "\n")
  }
  
  # Genes avec meilleure identite reciproque
  best_reciprocal <- results[results$Perc_Identity_Reciprocal == 100, ]
  cat("\n  --- GENES AVEC 100% IDENTITE RECIPROQUE ---\n")
  cat("  Nombre:", nrow(best_reciprocal), "\n")
  report_lines <- c(report_lines, "",
    "GENES AVEC 100% IDENTITE RECIPROQUE",
    paste0("  Nombre: ", nrow(best_reciprocal))
  )
  
  # Genes recuperes
  genes_recovered <- length(unique(results$Horse_Gene_Symbol))
  recovery_rate <- round(genes_recovered / found_horse * 100, 2)
  cat("\n  --- TAUX DE RECUPERATION GLOBAL ---\n")
  cat("  Genes cheval avec orthologue trouve:", genes_recovered, "/", found_horse, "\n")
  cat("  Taux de recuperation:", recovery_rate, "%\n")
  report_lines <- c(report_lines, "",
    "TAUX DE RECUPERATION GLOBAL",
    paste0("  Genes avec orthologue: ", genes_recovered, "/", found_horse),
    paste0("  Taux: ", recovery_rate, "%")
  )
  
  # Apercu des meilleurs resultats
  if (nrow(high_conf) > 0) {
    cat("\n  --- APERCU DES 5 MEILLEURS ORTHOLOGUES (confiance=1) ---\n")
    for (i in 1:min(5, nrow(high_conf))) {
      cat("  ", i, ".", high_conf$Horse_Gene_Symbol[i], "→", 
          high_conf$Human_Gene_Symbol_BioMart[i], "\n")
    }
  }
  
} else {
  cat("\n[6-7] AUCUN ORTHOLOGUE TROUVE...\n")
  cat("  ✗ Aucun orthologue valide trouve via biomaRt pour ces genes.\n")
  results <- data.frame(
    Horse_Gene_Symbol = character(),
    Message = character(),
    stringsAsFactors = FALSE
  )
  write.csv(results, OUTPUT_FILE, row.names = FALSE)
  report_lines <- c(report_lines, "",
    "RESULTAT: AUCUN ORTHOLOGUE TROUVE",
    "  Aucun orthologue valide trouve via biomaRt"
  )
}

# ==============================================================================
# ECRITURE DU RAPPORT
# ==============================================================================
report_lines <- c(report_lines, "",
  paste(rep("=", 80), collapse = ""),
  "FIN DU RAPPORT",
  paste(rep("=", 80), collapse = ""),
  paste0("Heure de fin: ", Sys.time())
)

writeLines(report_lines, REPORT_FILE)
cat("\n" , paste(rep("=", 80), collapse = ""), "\n", sep = "")
cat("RAPPORT COMPLET EXPORTE:\n")
cat("  ", REPORT_FILE, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("\n=== TERMINE ===\n")
