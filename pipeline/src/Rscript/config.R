# config.R
# Central configuration for the CellChat pipeline
#
# What it does:
#   Defines all project paths, biological conditions, cell-type subsetting rules,
#   color palettes, pairwise comparisons and script-specific parameters used by
#   every R script in this folder. A change here is picked up by all scripts on
#   the next run; individual R files do not need to be edited.
#
# Used by:
#   00_prepare_cellchat.R, 01_merge_cellchat.R, 02_compare_interactions.R,
#   03_diff_interactions.R, 04_circle_per_dataset.R, 05_circle_coarse_celltypes.R,
#   06_signaling_role_scatter.R, 07_net_similarity.R, 08_ranknet_heatmap.R,
#   09_bubble_dysfunctional.R, 10_dysfunctional_viz.R, 11_pathway_viz.R,
#   12_gene_expression_viz.R, 13_export_objects.R
#
# Edit the values below to change the conditions, subsetting, colors and
# script-specific parameters.

# ============================================================
# 1. PROJECT PATHS
# ============================================================

PROJECT_NAME <- "pipeline"
BASE_DIR <- "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# ============================================================
# 2. BIOLOGICAL CONDITIONS
# ============================================================

# Conditions to process (must match values in the condition metadata column)
CONDITIONS <- c("Healthy", "Crypto", "Immuno")

# Reference condition for pairwise comparisons
REFERENCE <- "Healthy"

# Seurat metadata column used for condition filtering and cell-type grouping
SEURAT_CONDITION_COL <- "detailed_group"  # column containing CONDITIONS values
SEURAT_CELL_TYPE_COL   <- "cell_type"        # column containing cell-type labels

# ============================================================
# 3. CELL-TYPE SUBSETTING
# ============================================================

# If TRUE, subset CellChat objects to TARGET_CELLTYPES before merging.
# If FALSE, keep all cell types.
DO_SUBSET <- FALSE

# Target cell types to retain for analysis (used only when DO_SUBSET = TRUE)
TARGET_CELLTYPES <- c("SSC", "Sertoli", "Leydig")

# ============================================================
# 4. OUTPUT NAMING
# ============================================================

# Output folder name under Result/data/merged and Result/plot:
#   "ALL" when DO_SUBSET is FALSE, or the joined target cell-type names otherwise.
MERGED_DIR_NAME <- if (DO_SUBSET && length(TARGET_CELLTYPES) > 0) {
  paste(TARGET_CELLTYPES, collapse = "_")
} else {
  "ALL"
}

# Multi-group output prefix, e.g. "healthy_vs_crypto_vs_immuno"
MERGE_PREFIX <- paste(tolower(CONDITIONS), collapse = "_vs_")

# ============================================================
# 5. PLOT COLORS
# ============================================================

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

# ============================================================
# 6. PAIRWISE COMPARISONS
# ============================================================

# Pairwise comparisons to run for scripts that only support two groups
PAIRWISE <- list(
  c("Healthy", "Crypto"),
  c("Healthy", "Immuno")
)

# ============================================================
# 7. SCRIPT-SPECIFIC PARAMETERS
# ============================================================

# --- 00_prepare_cellchat.R ---
PREP_INPUT_RDS          <- "Final_Annotated_Object_HUMAN.rds"
PREP_CELLCHAT_DB        <- "human"   # "human" or "mouse"
PREP_PROB_TYPE          <- "triMean"   # "triMean" or "truncatedMean"
PREP_PROB_TRIM          <- 0.30        # trim value when PREP_PROB_TYPE = "truncatedMean"
PREP_RAW_USE            <- TRUE
PREP_POPULATION_SIZE    <- FALSE
PREP_MIN_CELLS          <- 10          # min.cells for filterCommunication

# --- 01_merge_cellchat.R ---
# Controlled by DO_SUBSET and TARGET_CELLTYPES above.

# --- 02_compare_interactions.R ---
STEP5_DIR        <- "step5"
COMPARE_WIDTH    <- 600
COMPARE_HEIGHT   <- 500
COMPARE_RES      <- 120

# --- 03_diff_interactions.R ---
STEP6_DIR              <- "step6"
DIFF_CIRCLE_WIDTH      <- 1800
DIFF_CIRCLE_HEIGHT     <- 800
DIFF_CIRCLE_RES        <- 120

# --- 04_circle_per_dataset.R ---
STEP7_DIR                  <- "step7"
CIRCLE_PER_DS_WIDTH        <- 700     # width per panel
CIRCLE_PER_DS_HEIGHT       <- 700
CIRCLE_PER_DS_RES          <- 120
CIRCLE_EDGE_WIDTH_COUNT    <- 8
CIRCLE_EDGE_WIDTH_WEIGHT   <- 6

# --- 05_circle_coarse_celltypes.R ---
STEP8_DIR             <- "step8"
COARSE_WIDTH          <- 600          # width per panel
COARSE_HEIGHT         <- 600
COARSE_RES            <- 120
COARSE_EDGE_WIDTH     <- 12

# Fine cell type -> coarse group mapping.
# Any fine cell type not explicitly listed is assigned to COARSE_GROUP_DEFAULT.
COARSE_GROUP_DEFAULT <- "Somatic"
COARSE_GROUP_MAP <- c(
  "SSC"          = "SSC",
  "Spermatocyte" = "Germ",
  "Spermatid"    = "Germ",
  "Sertoli"      = "Sertoli",
  "Leydig"       = "Leydig",
  "Myoid"        = "Somatic",
  "Fibroblast"   = "Somatic",
  "Endothelial"  = "Somatic"
)
COARSE_GROUP_LEVELS <- c("SSC", "Germ", "Sertoli", "Leydig", "Somatic")

# --- 06_signaling_role_scatter.R ---
STEP9_DIR                      <- "step9"
SIGNALING_ROLE_CELL_TYPES      <- c("Sertoli", "Leydig", "SSC")
SIGNALING_ROLE_WIDTH_PER_PANEL <- 5
SIGNALING_ROLE_HEIGHT          <- 5
SIGNALING_ROLE_DPI             <- 150

# --- 07_net_similarity.R ---
STEP10_DIR            <- "step10"
NET_SIM_TYPE          <- "functional"  # "functional" or "structural"
NET_SIM_UMAP_METHOD   <- "uwot"
NET_SIM_SEED          <- 6
NET_SIM_WIDTH         <- 10
NET_SIM_HEIGHT        <- 8
NET_SIM_DPI           <- 150

# --- 08_ranknet_heatmap.R ---
STEP11_DIR               <- "step11"
RANKNET_WIDTH            <- 10
RANKNET_HEIGHT_PER_PW    <- 0.25        # inches per pathway for rankNet bar plots
HEATMAP_HEIGHT_PER_PW    <- 0.30        # inches per pathway for ComplexHeatmap
HEATMAP_HEIGHT_PER_LR    <- 0.18        # inches per L-R pair for L-R heatmaps
HEATMAP_PATTERNS         <- c("outgoing", "incoming", "all")

# --- 09_bubble_dysfunctional.R ---
STEP12_DIR                  <- "step12"
DYSFUNCTIONAL_SOURCES     <- c("Sertoli", "SSC")
DYSFUNCTIONAL_TARGETS     <- c("Sertoli", "SSC")
DYSFUNCTIONAL_LIGAND_LOGFC_UP   <- 0.05
DYSFUNCTIONAL_LIGAND_LOGFC_DOWN <- -0.05
DYSFUNCTIONAL_THRESH_PC     <- 0.1
DYSFUNCTIONAL_THRESH_FC     <- 0.05
DYSFUNCTIONAL_DO_FAST       <- TRUE

# --- 10_dysfunctional_viz.R ---
STEP13_DIR                   <- "step13"
DYSVIZ_BUBBLE_WIDTH          <- 14
DYSVIZ_BUBBLE_HEIGHT         <- 12
DYSVIZ_BUBBLE_DPI            <- 150
DYSVIZ_CHORD_SIZE            <- 10
DYSVIZ_WORDCLOUD_WIDTH       <- 10
DYSVIZ_WORDCLOUD_HEIGHT      <- 8

# --- 11_pathway_viz.R ---
STEP14_DIR               <- "step14"
PATHWAY_VIZ_RES          <- 150
PATHWAY_VIZ_WIDTH_1      <- 900     # single panel
PATHWAY_VIZ_HEIGHT_1     <- 900
PATHWAY_VIZ_WIDTH_2      <- 1800    # two panels
PATHWAY_VIZ_HEIGHT_2     <- 900
PATHWAY_VIZ_WIDTH_N      <- 900     # width per panel for 3+ panels
PATHWAY_VIZ_HEIGHT_N     <- 1800

# --- 12_gene_expression_viz.R ---
STEP15_DIR             <- "step15"
GENE_EXPR_VIZ_WIDTH    <- 14
GENE_EXPR_VIZ_HEIGHT   <- 10
GENE_EXPR_VIZ_DPI      <- 150
GENE_EXPR_VIZ_TYPE     <- "violin"   # "violin" or "dot"

# --- 13_export_objects.R ---
STEP16_DIR <- "step16"
