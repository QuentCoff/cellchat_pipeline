# config.R
# Central configuration for the ALL_clusters CellChat pipeline
# Edit the values below to change the conditions and reference.

PROJECT_NAME <- "ALL_clusters"
BASE_DIR <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# Conditions to process (must match values in seurat_obj$detailed_group)
CONDITIONS <- c("Healthy", "Crypto", "Immuno")

# Reference condition for pairwise comparisons
REFERENCE <- "Healthy"

# If TRUE, subset CellChat objects to TARGET_CELLTYPES before merging.
# If FALSE, keep all cell types.
DO_SUBSET <- FALSE

# Target cell types to retain for analysis (all clusters)
TARGET_CELLTYPES <- c("SSC", "Sertoli", "Leydig")

# Output folder name under Result/data/merged and Result/plot:
#   "ALL" when DO_SUBSET is FALSE, or the joined target cell-type names otherwise.
MERGED_DIR_NAME <- if (DO_SUBSET && length(TARGET_CELLTYPES) > 0) {
  paste(TARGET_CELLTYPES, collapse = "_")
} else {
  "ALL"
}

# Fixed colors for cell type plots
CELLTYPE_COLORS <- c(
  "SSC"          = "#417505",
  "Spermatocyte" = "#4A90E2",
  "Spermatid"    = "#D0021B",
  "Sertoli"      = "#F5A623",
  "Leydig"       = "#880FF3",
  "Myoid"        = "#7ED321",
  "Fibroblast"   = "#5A5A5A",
  "Endothelial"  = "#F16B54"
)

# Fixed colors for condition/dataset plots
CONDITION_COLORS <- c(
  "Healthy" = "#e40e03",  # red
  "Crypto"  = "#2692a5",  # blue
  "Immuno"  = "#b0d3b8"   # green
)

# Pairwise comparisons to run for scripts that only support two groups
PAIRWISE <- list(
  c("Healthy", "Crypto"),
  c("Healthy", "Immuno")
)

# Multi-group output prefix, e.g. "healthy_vs_crypto_vs_immuno"
MERGE_PREFIX <- paste(tolower(CONDITIONS), collapse = "_vs_")
