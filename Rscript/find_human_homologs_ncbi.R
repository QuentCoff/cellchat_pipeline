#!/usr/bin/env Rscript
# Script pour trouver les homologues humains via NCBI Datasets API
# Map les symboles de genes cheval vers leurs orthologues humains

library(httr)
library(jsonlite)

# Lecture des genes
genes_file <- '/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/gene_export.csv'
genes_export <- read.csv(genes_file, header = FALSE, stringsAsFactors = FALSE)
gene_symbols <- unique(genes_export$V1)

cat("=== MAPPING VIA NCBI DATASETS ===\n")
cat(paste0("Nombre de genes a mapper: ", length(gene_symbols), "\n\n"))

# Fonction pour chercher un gene cheval et trouver son orthologue humain
find_ortholog <- function(horse_symbol) {
  result <- data.frame(
    horse_symbol = horse_symbol,
    horse_gene_id = NA,
    horse_description = NA,
    human_symbol = NA,
    human_gene_id = NA,
    human_description = NA,
    ortholog_type = NA,
    stringsAsFactors = FALSE
  )
  
  tryCatch({
    # 1. Chercher le gene cheval dans NCBI (taxon 9796 = horse)
    search_url <- paste0(
      "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/symbol/", 
      URLencode(horse_symbol), 
      "/taxon/9796"
    )
    
    response <- GET(search_url, timeout(30))
    
    if (status_code(response) != 200) {
      return(result)
    }
    
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))
    
    if (is.null(data$genes) || length(data$genes) == 0) {
      return(result)
    }
    
    # Extraire info gene cheval
    horse_gene <- data$genes[[1]]$gene
    result$horse_gene_id <- horse_gene$gene_id
    if (!is.null(horse_gene$description)) {
      result$horse_description <- horse_gene$description
    }
    
    # 2. Chercher les orthologues
    if (!is.null(data$genes[[1]]$orthologs) && length(data$genes[[1]]$orthologs) > 0) {
      orthologs <- data$genes[[1]]$orthologs
      
      # Chercher l'orthologue humain (taxon 9606)
      for (ortho in orthologs) {
        if (ortho$gene$taxname == "Homo sapiens" || ortho$gene$tax_id == 9606) {
          result$human_symbol <- ortho$gene$symbol
          result$human_gene_id <- ortho$gene$gene_id
          if (!is.null(ortho$gene$description)) {
            result$human_description <- ortho$gene$description
          }
          if (!is.null(ortho$method)) {
            result$ortholog_type <- paste(ortho$method, collapse = "; ")
          }
          break
        }
      }
    }
    
  }, error = function(e) {
    # Erreur silencieuse, retourne le resultat vide
  })
  
  return(result)
}

# Traitement avec delai pour respecter les limites NCBI (3 requetes/sec)
results_list <- list()
total <- length(gene_symbols)

for (i in seq_along(gene_symbols)) {
  symbol <- gene_symbols[i]
  
  if (i %% 10 == 1 || i == total) {
    cat(paste0("Traitement: ", i, "/", total, " (", symbol, ")\n"))
  }
  
  result <- find_ortholog(symbol)
  results_list[[i]] <- result
  
  # Delai pour respecter la limite de 3 requetes par seconde
  Sys.sleep(0.4)
}

# Combinaison des resultats
results <- do.call(rbind, results_list)

# Statistiques detaillees
cat("\n=== STATISTIQUES DETAILLEES ===\n")

total_input <- length(gene_symbols)
found_in_ncbi <- sum(!is.na(results$horse_gene_id))
with_human_ortholog <- sum(!is.na(results$human_symbol))
not_in_ncbi <- sum(is.na(results$horse_gene_id))
in_ncbi_no_ortholog <- sum(!is.na(results$horse_gene_id) & is.na(results$human_symbol))

# Calculs de pourcentages
pct_match <- round(with_human_ortholog / total_input * 100, 2)
pct_loss_total <- round((total_input - with_human_ortholog) / total_input * 100, 2)
pct_not_in_ncbi <- round(not_in_ncbi / total_input * 100, 2)
pct_no_ortholog <- round(in_ncbi_no_ortholog / total_input * 100, 2)

cat(paste0("=== DONNEES DE BASE ===\n"))
cat(paste0("Nombre de genes en entree (gene_export.csv): ", total_input, "\n"))

cat(paste0("\n=== MAPPING NCBI ===\n"))
cat(paste0("Genes trouves dans NCBI (cheval): ", found_in_ncbi, " (", round(found_in_ncbi/total_input*100, 2), "%)\n"))
cat(paste0("  -> Non trouves dans NCBI: ", not_in_ncbi, " (", pct_not_in_ncbi, "%)\n"))

cat(paste0("\n=== ORTHOLOGIE ===\n"))
cat(paste0("Genes avec orthologue humain: ", with_human_ortholog, " (", pct_match, "%)\n"))
cat(paste0("  -> Trouves dans NCBI mais SANS orthologue humain: ", in_ncbi_no_ortholog, " (", pct_no_ortholog, "%)\n"))

cat(paste0("\n=== PERTE TOTALE ===\n"))
cat(paste0("Genes PERDUS (sans orthologue humain): ", total_input - with_human_ortholog, " (", pct_loss_total, "%)\n"))
cat(paste0("Genes CONSERVES (avec orthologue): ", with_human_ortholog, " (", pct_match, "%)\n"))

# Creation d'un resume pour export
summary_stats <- data.frame(
  metric = c(
    "Total genes input",
    "Found in NCBI (horse)",
    "Not found in NCBI",
    "With human ortholog",
    "Found in NCBI but no human ortholog",
    "Total genes lost",
    "Total genes conserved"
  ),
  count = c(
    total_input,
    found_in_ncbi,
    not_in_ncbi,
    with_human_ortholog,
    in_ncbi_no_ortholog,
    total_input - with_human_ortholog,
    with_human_ortholog
  ),
  percentage = c(
    "100%",
    paste0(round(found_in_ncbi/total_input*100, 2), "%"),
    paste0(pct_not_in_ncbi, "%"),
    paste0(pct_match, "%"),
    paste0(pct_no_ortholog, "%"),
    paste0(pct_loss_total, "%"),
    paste0(pct_match, "%")
  )
)

# Export
cat("\n=== EXPORT ===\n")
out_dir <- '/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(out_dir, "gene_export_ncbi_orthologs.csv")
write.csv(results, output_file, row.names = FALSE, na = "")
cat(paste0("Fichier exporte: ", output_file, "\n"))

# Export du resume des statistiques
summary_file <- file.path(out_dir, "gene_export_ncbi_summary.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(paste0("Resume statistique exporte: ", summary_file, "\n"))

# Liste des genes non mappes
not_mapped <- results$horse_symbol[is.na(results$human_symbol)]
if (length(not_mapped) > 0) {
  missing_file <- file.path(out_dir, "gene_export_not_mapped_ncbi.csv")
  write.csv(data.frame(horse_symbol = not_mapped), missing_file, row.names = FALSE)
  cat(paste0("Genes non mappes exportes: ", missing_file, "\n"))
}

cat("\nTermine!\n")
