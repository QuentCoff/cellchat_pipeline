#!/bin/bash
# run_all.sh
# Full Vtruncated pipeline: truncatedMean (trim=0.04)
# Procedure 1 (Healthy + Crypto) → Procedure 2 (merge + steps 02-13)
#
# Usage:
#   bash run_all.sh              # P1 + P2
#   bash run_all.sh skip_p1      # P2 only

set -euo pipefail

SKIP_P1=false
if [[ "${1:-}" == "skip_p1" ]]; then
  SKIP_P1=true
fi

BASE_DIR="/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"
P1_SCRIPT="${BASE_DIR}/Vtruncated/src/Rscript/00_prepare_cellchat.R"
P2_DIR="${BASE_DIR}/Vtruncated/src/Rscript"

echo "============================================"
echo "  Pipeline: Vtruncated (truncatedMean, trim=0.04)"
echo "============================================"

# --- Procedure 1 ---
if [[ "$SKIP_P1" == true ]]; then
  echo ""
  echo "=== Skipping Procedure 1 ==="
else
  echo ""
  echo "=== Procedure 1: Healthy ==="
  Rscript "${P1_SCRIPT}" Healthy

  echo ""
  echo "=== Procedure 1: Crypto ==="
  Rscript "${P1_SCRIPT}" Crypto
fi

# --- Procedure 2 ---
echo ""
echo "=== Procedure 2: Merge ==="
Rscript "${P2_DIR}/01_merge_cellchat.R" healthy crypto

echo ""
echo "=== Step 5: compare interactions ==="
Rscript "${P2_DIR}/02_compare_interactions.R" healthy crypto

echo ""
echo "=== Step 6: diff interactions ==="
Rscript "${P2_DIR}/03_diff_interactions.R" healthy crypto

echo ""
echo "=== Step 7A: circle per dataset ==="
Rscript "${P2_DIR}/04_circle_per_dataset.R" healthy crypto

echo ""
echo "=== Step 7B: circle coarse celltypes ==="
Rscript "${P2_DIR}/05_circle_coarse_celltypes.R" healthy crypto

echo ""
echo "=== Step 9: signaling role scatter ==="
Rscript "${P2_DIR}/06_signaling_role_scatter.R" healthy crypto

echo ""
echo "=== Step 10: net similarity ==="
Rscript "${P2_DIR}/07_net_similarity.R" healthy crypto

echo ""
echo "=== Step 10bis: riverplot clusters ==="
Rscript "${P2_DIR}/07b_riverplot_clusters.R" healthy crypto

echo ""
echo "=== Step 10bis: robustness (100 iterations UMAP + Wilcoxon) ==="
Rscript "${P2_DIR}/07_robustness.R" functional 100

echo ""
echo "=== Step 11: rankNet + heatmap ==="
Rscript "${P2_DIR}/08_ranknet_heatmap.R" healthy crypto

echo ""
echo "=== Step 12: bubble dysfunctional ==="
Rscript "${P2_DIR}/09_bubble_dysfunctional.R" healthy crypto

echo ""
echo "=== Step 13: dysfunctional viz ==="
Rscript "${P2_DIR}/10_dysfunctional_viz.R" healthy crypto

echo ""
echo "=== Step 14: pathway viz ==="
Rscript "${P2_DIR}/11_pathway_viz.R" healthy crypto

echo ""
echo "=== Step 15: gene expression viz ==="
Rscript "${P2_DIR}/12_gene_expression_viz.R" healthy crypto

echo ""
echo "=== Step 16: export objects ==="
Rscript "${P2_DIR}/13_export_objects.R" healthy crypto

echo ""
echo "============================================"
echo "  Pipeline complete!"
echo "  Data:   ${BASE_DIR}/Vtruncated/Result/data/"
echo "  Plots:  ${BASE_DIR}/Vtruncated/Result/plot/"
echo "============================================"
