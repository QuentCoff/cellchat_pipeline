#!/usr/bin/env bash
# run_all.sh
# Pipeline complet SSC + Sertoli : Procedure 1 + Procedure 2
#
# Usage:
#   bash run_all.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat"

# ============================================
# Step 0: Create output directories
# ============================================
mkdir -p "$BASE_DIR/results/SSC_Sertoli/Procedure_1/healthy"
mkdir -p "$BASE_DIR/results/SSC_Sertoli/Procedure_1/crypto"
mkdir -p "$BASE_DIR/results/SSC_Sertoli/Procedure_2/healthy_vs_crypto/comparison"

# ============================================
# Step 1: Procedure 1 — Create CellChat objects per condition
# ============================================
echo "========================================"
echo "Step 1a: CellChat for Healthy (Descended)"
echo "========================================"
Rscript "$SCRIPT_DIR/01_prepare_cellchat.R" Healthy

echo ""
echo "========================================"
echo "Step 1b: CellChat for Crypto"
echo "========================================"
Rscript "$SCRIPT_DIR/01_prepare_cellchat.R" Crypto

# ============================================
# Step 2: Procedure 2 — Merge and compare
# ============================================
echo ""
echo "========================================"
echo "Step 2: Merge CellChat objects"
echo "========================================"
Rscript "$SCRIPT_DIR/02_merge_cellchat.R" healthy crypto

echo ""
echo "========================================"
echo "Step 3: Compare interactions (count + weight)"
echo "========================================"
Rscript "$SCRIPT_DIR/03_compare_interactions.R" healthy crypto

echo ""
echo "========================================"
echo "Step 4: Differential interactions (circle + heatmap)"
echo "========================================"
Rscript "$SCRIPT_DIR/04_diff_interactions.R" healthy crypto

echo ""
echo "========================================"
echo "Step 5: Circle plots per dataset (normalized scale)"
echo "========================================"
Rscript "$SCRIPT_DIR/05_circle_per_dataset.R" healthy crypto

echo ""
echo "========================================"
echo "=== PIPELINE SSC+SERTOLI COMPLETE ==="
echo "Results in: $BASE_DIR/results/SSC_Sertoli/"
echo "========================================"
