# Rscript pipeline overview

This folder contains the R scripts that run the CellChat analysis for the `ALL_clusters` project.
Each script corresponds to a CellChat tutorial step and is designed to be run independently or via `../run_all.sh`.

All scripts source `config.R` for paths, conditions, colors and parameters.

## Script map

| Order | Script | CellChat step | What it does | Inputs | Outputs | Previous step |
|-------|--------|---------------|--------------|--------|---------|---------------|
| 1 | `00_prepare_cellchat.R` | Procedure 1 | Loads the Seurat object, splits it by condition, creates one CellChat object per condition and runs preprocessing + communication inference. | `data_input/<PREP_INPUT_RDS>` | `Result/data/cellchat_prep/<condition>/cellchat_<condition>.rds` | None |
| 2 | `01_merge_cellchat.R` | Procedure 2 Steps 1-3 | Merges per-condition CellChat objects into a multi-group object and into all pairwise objects defined in `config.R`. | Per-condition RDS from step 1 | `Result/data/merged/<MERGED_DIR_NAME>/<merge>/cellchat_merged_*.RData` and `cellchat_object.list_*.RData` | `00_prepare_cellchat.R` |
| 3 | `02_compare_interactions.R` | Procedure 2 Step 5 | Compares total number and strength of interactions across all conditions (bar plots). | Multi-group merged RData | `Result/plot/<MERGED_DIR_NAME>/step5/compare_interactions_*.png` | `01_merge_cellchat.R` |
| 4 | `03_diff_interactions.R` | Procedure 2 Step 6 | Differential interaction circle plots between each pairwise comparison. | Pairwise merged RData | `Result/plot/<MERGED_DIR_NAME>/<pair>/step6/diff_circle_*.png` | `01_merge_cellchat.R` |
| 5 | `04_circle_per_dataset.R` | Procedure 2 Step 7 | Circle plots per dataset on a common scale (count and strength). | Multi-group `object.list` RData | `Result/plot/<MERGED_DIR_NAME>/step7/circle_per_dataset_*.png` | `01_merge_cellchat.R` |
| 6 | `05_circle_coarse_celltypes.R` | Procedure 2 Step 8 | Same as step 7 but with coarse cell-type groups defined in `config.R`. | Multi-group `object.list` RData | `Result/plot/<MERGED_DIR_NAME>/step8/circle_coarse_*.png` | `01_merge_cellchat.R` |
| 7 | `06_signaling_role_scatter.R` | Procedure 2 Step 9 | 9A: signaling-role scatter plots per dataset. 9B: signaling-changes scatter plots per cell type for each pairwise comparison. | Multi-group + pairwise RData | `Result/plot/<MERGED_DIR_NAME>/step9/signaling_role_scatter.png` and `step9/<pair>/signaling_changes_*.png` | `01_merge_cellchat.R` |
| 8 | `07_net_similarity.R` | Procedure 2 Step 10 | Computes signaling-pathway similarity (functional or structural), runs UMAP + clustering and rank similarity for each pairwise comparison. | Pairwise merged RData | `Result/plot/<MERGED_DIR_NAME>/step10/<pair>/embedding_*.png`, `rank_similarity_*.png`, `cellchat_clustered_*.RData` | `01_merge_cellchat.R` |
| 9 | `08_ranknet_heatmap.R` | Procedure 2 Step 11 | RankNet information flow plots and outgoing/incoming/all ComplexHeatmaps per signaling pathway. | Multi-group RData | `Result/plot/<MERGED_DIR_NAME>/step11/ranknet_*.png`, `heatmap_*.pdf` | `01_merge_cellchat.R` |
| 10 | `09_bubble_dysfunctional.R` | Procedure 2 Step 12 | Identifies up- and down-regulated L-R pairs for each pairwise comparison and saves DEG-related CSVs. | Pairwise RData | `Result/plot/<MERGED_DIR_NAME>/<pair>/step12/bubble_*.png`, `net_*.csv`, `cellchat_deg.RData` | `01_merge_cellchat.R` |
| 11 | `10_dysfunctional_viz.R` | Procedure 2 Step 13 | Visualizes up/down-regulated signaling with bubble, chord and wordcloud plots. | Pairwise RData + step 12 CSVs | `Result/plot/<MERGED_DIR_NAME>/<pair>/step13/*.png`, `*.pdf` | `09_bubble_dysfunctional.R` |
| 12 | `11_pathway_viz.R` | Procedure 2 Step 14 | Circle plots for individual signaling pathways across all conditions. Accepts optional pathway list as argument. | Multi-group `object.list` RData | `Result/plot/<MERGED_DIR_NAME>/step14/circle_<pathway>.png` | `01_merge_cellchat.R` |
| 13 | `12_gene_expression_viz.R` | Procedure 2 Step 15 | Violin plots of gene-expression distribution for signaling genes per pathway. Accepts optional pathway list as argument. | Merged or DEG CellChat RData | `Result/plot/<MERGED_DIR_NAME>/step15/<type>_<pathway>.png` | `09_bubble_dysfunctional.R` (optional) |
| 14 | `13_export_objects.R` | Procedure 2 Step 16 | Copies final `object.list` and merged CellChat objects to the export folder. | Multi-group RData | `Result/plot/<MERGED_DIR_NAME>/step16/cellchat_*.RData` | `01_merge_cellchat.R` |

## Running a single script

```bash
micromamba activate env_cellchat
Rscript /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/ALL_clusters/src/Rscript/00_prepare_cellchat.R
```

## Running the full pipeline

```bash
micromamba activate env_cellchat
bash /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/ALL_clusters/src/run_all.sh
```

## Passing a pathway list

Some scripts accept a comma-separated pathway list:

```bash
Rscript 11_pathway_viz.R CXCL,BMP,WNT
Rscript 12_gene_expression_viz.R CXCL,BMP,WNT
```

If no pathways are provided, all significant pathways are processed.
