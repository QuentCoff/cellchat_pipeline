#!/bin/bash
# One-line setup for the CellChat pipeline on GLiCID.
# Usage:
#   bash setup.sh [ENV_NAME]
#   Default ENV_NAME is "env_cellchat".

# Exit immediately if any command fails.
set -e

# In a non-interactive shell the 'micromamba' shell function is not inherited.
# MAMBA_EXE is exported by .bashrc and points to the real executable.
MAMBA_CMD="${MAMBA_EXE:-micromamba}"
if ! command -v "$MAMBA_CMD" &> /dev/null; then
  echo "Error: micromamba not found. Please source ~/.bashrc or set MAMBA_EXE."
  exit 1
fi

# Directory where this script is located.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Environment name can be passed as the first argument; defaults to env_cellchat.
ENV_NAME="${1:-env_cellchat}"

echo "==> Creating/Updating micromamba environment: ${ENV_NAME}"
echo "    Usage: bash setup.sh [ENV_NAME]    (default: env_cellchat)"

# Create the environment if it does not exist; otherwise update it from the YAML spec.
if "$MAMBA_CMD" env list | grep -q "${ENV_NAME}"; then
  echo "Environment ${ENV_NAME} already exists. Updating..."
  "$MAMBA_CMD" update -n "${ENV_NAME}" -f "${SCRIPT_DIR}/environment.yml" -y
else
  echo "Environment ${ENV_NAME} does not exist. Creating..."
  "$MAMBA_CMD" create -n "${ENV_NAME}" -f "${SCRIPT_DIR}/environment.yml" -y
fi

# Install CellChat from GitHub inside the target environment.
# We use 'micromamba run' so that conda compilers are available in PATH.
echo "==> Installing CellChat from GitHub"
"$MAMBA_CMD" run -n "${ENV_NAME}" Rscript "${SCRIPT_DIR}/install_cellchat.R"

echo "==> Setup complete. Activate with: micromamba activate ${ENV_NAME}"
