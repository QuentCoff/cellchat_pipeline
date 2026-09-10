#!/bin/bash
# run_all.sh
# Full ALL_clusters CellChat pipeline orchestrator.
#
# What it does:
#   Runs the complete analysis in order:
#     Procedure 1: per-condition CellChat objects (00_prepare_cellchat.R)
#     Procedure 2: multi-group merge (01) plus all visualization steps (02-13).
#   Pairwise-only scripts (03, 07, 09, 10) are called after the multi-group block.
#   Steps currently disabled in this runner are commented out; they can be enabled
#   by uncommenting the corresponding Rscript line.
#
# Inputs:
#   - Seurat RDS file: data_input/<PREP_INPUT_RDS>
#   - Configuration file: ALL_clusters/src/Rscript/config.R
#
# Outputs:
#   - Result/data/  : intermediate and final RData/RDS objects
#   - Result/plot/  : all generated figures
#
# Usage:
#   bash run_all.sh              # Procedure 1 + Procedure 2
#   bash run_all.sh skip_p1      # Procedure 2 only (re-use existing prep objects)

set -euo pipefail

SKIP_P1=false
if [[ "${1:-}" == "skip_p1" ]]; then
  SKIP_P1=true
fi

BASE_DIR="/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"
P1_SCRIPT="${BASE_DIR}/ALL_clusters/src/Rscript/00_prepare_cellchat.R"
P2_DIR="${BASE_DIR}/ALL_clusters/src/Rscript"

echo "============================================"
echo "  Pipeline: ALL_clusters (trimean)"
echo "  Conditions: Healthy, Crypto, Immuno"
echo "============================================"

# --- Procedure 1 ---
if [[ "$SKIP_P1" == true ]]; then
  echo ""
  echo "=== Skipping Procedure 1 ==="
else
  echo ""
  echo "=== Procedure 1: prepare all conditions ==="
  Rscript "${P1_SCRIPT}"
fi

# --- Procedure 2: multi-group ---
echo ""
echo "=== Procedure 2: merge all conditions + pairwise merges ==="
Rscript "${P2_DIR}/01_merge_cellchat.R"

echo ""
echo "=== Step 5: compare interactions (multi-group) ==="
Rscript "${P2_DIR}/02_compare_interactions.R"

echo ""
echo "=== Step 7A: circle per dataset (multi-group) ==="
Rscript "${P2_DIR}/04_circle_per_dataset.R"

#echo ""
#echo "=== Step 7B: circle coarse celltypes (multi-group) ==="
#Rscript "${P2_DIR}/05_circle_coarse_celltypes.R"

echo ""
echo "=== Step 9: signaling role scatter (multi-group) ==="
Rscript "${P2_DIR}/06_signaling_role_scatter.R"

echo ""
echo "=== Step 11: rankNet + heatmap (multi-group) ==="
Rscript "${P2_DIR}/08_ranknet_heatmap.R"

echo ""
echo "=== Step 14: pathway viz (multi-group) ==="
Rscript "${P2_DIR}/11_pathway_viz.R"

echo ""
echo "=== Step 15: gene expression viz (multi-group) ==="
Rscript "${P2_DIR}/12_gene_expression_viz.R"

echo ""
echo "=== Step 16: export objects (multi-group) ==="
Rscript "${P2_DIR}/13_export_objects.R"

# --- Pairwise comparisons ---
echo ""
echo "=== Step 6: diff interactions (pairwise) ==="
Rscript "${P2_DIR}/03_diff_interactions.R"

echo ""
echo "=== Step 10: net similarity (pairwise) ==="
Rscript "${P2_DIR}/07_net_similarity.R"

#echo ""
#echo "=== Step 10bis: robustness (100 iterations UMAP + Wilcoxon, pairwise) ==="
#Rscript "${P2_DIR}/07_robustness.R"

#echo ""
#echo "=== Step 12: bubble dysfunctional (pairwise) ==="
#Rscript "${P2_DIR}/09_bubble_dysfunctional.R"

#echo ""
#echo "=== Step 13: dysfunctional viz (pairwise) ==="
#Rscript "${P2_DIR}/10_dysfunctional_viz.R"

echo ""
echo "============================================"
echo "  Pipeline complete!"
echo "  Data:   ${BASE_DIR}/ALL_clusters/Result/data/"
echo "  Plots:  ${BASE_DIR}/ALL_clusters/Result/plot/"
echo "============================================"
