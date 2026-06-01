#!/usr/bin/env bash
# run_pipeline.sh
# Run CellChat Procedure 1 (prepare + visualize) in one shot
#
# Usage:
#   bash run_pipeline.sh <detailed_group> [receiver_indices]
#
# Examples:
#   bash run_pipeline.sh Healthy
#   bash run_pipeline.sh Healthy 1,2,3,4
#   bash run_pipeline.sh Crypto       # regroupe Uni-Crypto + Bi-Crypto
#   bash run_pipeline.sh Immuno

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash run_pipeline.sh <detailed_group> [receiver_indices]"
    echo "Examples:"
    echo "  bash run_pipeline.sh Healthy"
    echo "  bash run_pipeline.sh Healthy 1,2,3,4"
    exit 1
fi

DETAILED_GROUP="$1"
RECEIVER_IDX="${2:-}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "========================================"
echo "CellChat Pipeline"
echo "  Group  : ${DETAILED_GROUP}"
echo "  Receivers: ${RECEIVER_IDX:-none}"
echo "========================================"
echo ""

# Step 1: Prepare CellChat object
echo "[1/2] Running 01_prepare_cellchat.R..."
Rscript "${SCRIPT_DIR}/01_prepare_cellchat.R" "${DETAILED_GROUP}"
echo "[1/2] Done."
echo ""

# Step 2: Visualize pathways
echo "[2/2] Running 02_visualize_pathways.R..."
if [ -n "${RECEIVER_IDX}" ]; then
    Rscript "${SCRIPT_DIR}/02_visualize_pathways.R" "$(echo "${DETAILED_GROUP}" | tr '[:upper:]' '[:lower:]')" "${RECEIVER_IDX}"
else
    Rscript "${SCRIPT_DIR}/02_visualize_pathways.R" "$(echo "${DETAILED_GROUP}" | tr '[:upper:]' '[:lower:]')"
fi
echo "[2/2] Done."
echo ""
echo "========================================"
echo "Pipeline complete for: ${DETAILED_GROUP}"
echo "========================================"
