# CellChat Pipeline

CellChat-based analysis pipeline for horse testicular single-cell data.

## Clone the repository

```bash
git clone https://github.com/QuentCoff/cellchat_pipeline.git
cd cellchat_pipeline
```

## Environment setup

Run the provided setup script from a GLiCID login node or compute node. It creates the micromamba environment and installs CellChat from GitHub.

```bash
bash setup.sh
```

To use a custom environment name:

```bash
bash setup.sh my_env_name
```

## Prepare a compute session

Start an interactive compute node on GLiCID:

```bash
srun --cluster=nautilus -N1 --qos=quick --cpus-per-task=2 --mem=100G --time=2:00:00 --pty bash
```

Activate the environment and load the required modules:

```bash
micromamba activate my_env_name
```

## Run the full pipeline

To run all numbered R scripts in order:

```bash
bash ALL_clusters/src/run_all.sh
```

## Run a specific script

Example:

```bash
Rscript ALL_clusters/src/Rscript/00_prepare_cellchat.R
```

Replace the script name with the step you want to run.

## Project structure

```
cellchat_pipeline/
├── README.md                  # This file
├── setup.sh                   # One-line environment + CellChat installer
├── environment.yml            # Micromamba environment specification
├── install_cellchat.R         # Installs the CellChat R package from GitHub
├── .gitignore                 # Excludes generated results (.rds, .png, ...)
├── data_input/                # Static input data
│   ├── genes_manquants.txt
│   └── ortho_complete_horse_human.csv
├── docs/                      # Additional documentation
│   └── clustering_pathways.md
└── ALL_clusters/              # Main analysis pipeline
    ├── src/Rscript/           # Numbered R scripts (run in order)
    │   ├── config.R           # Global parameters (colors, comparisons, paths)
    │   ├── 00_prepare_cellchat.R
    │   ├── 01_merge_cellchat.R
    │   ├── 02_compare_interactions.R
    │   ├── 03_diff_interactions.R
    │   ├── 04_circle_per_dataset.R
    │   ├── 05_circle_coarse_celltypes.R
    │   ├── 06_signaling_role_scatter.R
    │   ├── 07_net_similarity.R
    │   ├── 07_robustness.R
    │   ├── 08_ranknet_heatmap.R
    │   ├── 09_bubble_dysfunctional.R
    │   ├── 10_dysfunctional_viz.R
    │   ├── 11_pathway_viz.R
    │   ├── 12_gene_expression_viz.R
    │   ├── 13_export_objects.R
    │   └── ranknet_source.R
    └── Result/                # Generated outputs (not tracked by Git)
        ├── data/
        │   ├── cellchat_prep/
        │   └── merged/
        └── plot/
```

Step numbers in `Result/plot/` correspond to the numbered analysis scripts (e.g. `step9/` contains outputs from the signaling-changes part of `06_signaling_role_scatter.R`).


