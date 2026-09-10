# Configuration guide for the ALL_clusters CellChat pipeline

This document describes every parameter defined in `src/Rscript/config.R`.
Each script reads the relevant variables from this central file, so you can
parameterize the whole pipeline without editing the R scripts themselves.

---

## 1. Project paths

| Parameter | Value / example | Description |
|-----------|-----------------|-------------|
| `PROJECT_NAME` | `"ALL_clusters"` | Name of the analysis folder under `BASE_DIR`. |
| `BASE_DIR` | `"/LAB-DATA/GLiCID/users/.../Cellchat"` | Root directory of the project. All input / output paths are derived from it. |

These two variables are used by every script to locate data and save results.

---

## 2. Biological conditions

| Parameter | Value / example | Description |
|-----------|-----------------|-------------|
| `CONDITIONS` | `c("Healthy", "Crypto", "Immuno")` | Biological groups to compare. They must match the values stored in `seurat_obj$detailed_group`. |
| `REFERENCE` | `"Healthy"` | Reference condition for pairwise comparisons. |
| `SEURAT_CONDITION_COL` | `"detailed_group"` | Name of the Seurat metadata column containing the condition labels (`CONDITIONS`). |
| `SEURAT_CELL_TYPE_COL` | `"cell_type"` | Name of the Seurat metadata column containing **all cluster / cell-type labels that you want to analyze**. This column must already contain the final annotations (e.g. `Sertoli`, `Leydig`, `Myoid`…), not raw cluster numbers. |

If you rename the conditions or cell-type column in the Seurat object, update
`CONDITIONS`, `REFERENCE`, `SEURAT_CONDITION_COL` and `SEURAT_CELL_TYPE_COL`
accordingly.

---

## 3. Cell-type subsetting

| Parameter | Value / example | Description |
|-----------|-----------------|-------------|
| `DO_SUBSET` | `FALSE` | If `TRUE`, each per-condition CellChat object is subset to `TARGET_CELLTYPES` before merging. If `FALSE`, all cell types are kept. |
| `TARGET_CELLTYPES` | `c("SSC", "Sertoli", "Leydig")` | Cell types to retain when `DO_SUBSET = TRUE`. Must be values present in `SEURAT_CELL_TYPE_COL`. |

`MERGED_DIR_NAME` is automatically derived:

- `DO_SUBSET = FALSE`  → `"ALL"`
- `DO_SUBSET = TRUE`   → paste of `TARGET_CELLTYPES` (e.g. `"SSC_Sertoli_Leydig"`)

---

## 4. Output naming

| Parameter | Derivation | Description |
|-----------|------------|-------------|
| `MERGED_DIR_NAME` | `ifelse(DO_SUBSET, paste(TARGET_CELLTYPES, collapse="_"), "ALL")` | Sub-folder name under `Result/data/merged` and `Result/plot`. |
| `MERGE_PREFIX` | `paste(tolower(CONDITIONS), collapse="_vs_")` | Multi-group output prefix, e.g. `healthy_vs_crypto_vs_immuno`. |

Do not edit these directly unless you want to override the automatic naming.

---

## 5. Plot colors

### Cell-type colors

| Parameter | Description |
|-----------|-------------|
| `CELLTYPE_COLORS` | Named vector `c("CellType" = "#hex")` used by all circle plots, heatmaps and chord diagrams. |

Every cell type present in the analysis must have an entry here. If a name is
missing, scripts that rely on colors may fail or fall back to default colors.

### Condition / dataset colors

| Parameter | Description |
|-----------|-------------|
| `CONDITION_COLORS` | Named vector `c("Condition" = "#hex")` used for condition-level plots. |

The names must exactly match the values in `CONDITIONS`.

---

## 6. Pairwise comparisons

| Parameter | Value / example | Description |
|-----------|-----------------|-------------|
| `PAIRWISE` | `list(c("Healthy", "Crypto"), c("Healthy", "Immuno"))` | List of two-condition vectors for scripts that only support pairwise comparisons. |

Each inner vector must contain two values that also appear in `CONDITIONS`.

---

## 7. Script-specific parameters

### `00_prepare_cellchat.R`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `PREP_INPUT_RDS` | `"Final_Annotated_Object_HUMAN.rds"` | Input Seurat object, read from `BASE_DIR/data_input/`. |
| `SEURAT_CONDITION_COL` | `"detailed_group"` | Seurat metadata column used to filter conditions. |
| `SEURAT_CELL_TYPE_COL` | `"cell_type"` | Seurat metadata column containing all cluster / cell-type labels to analyze, passed to `createCellChat(group.by = ...)`. |
| `PREP_CELLCHAT_DB` | `"human"` | Database to use: `"human"` or `"mouse"`. |
| `PREP_PROB_TYPE` | `"triMean"` | Method for `computeCommunProb`: `"triMean"` or `"truncatedMean"`. |
| `PREP_PROB_TRIM` | `0.30` | Trim value used only when `PREP_PROB_TYPE = "truncatedMean"`. |
| `PREP_RAW_USE` | `TRUE` | `raw.use` argument for `computeCommunProb`. |
| `PREP_POPULATION_SIZE` | `FALSE` | `population.size` argument for `computeCommunProb`. |
| `PREP_MIN_CELLS` | `10` | `min.cells` argument for `filterCommunication`. |

### `01_merge_cellchat.R`

No dedicated variables. This script is controlled by `DO_SUBSET` and
`TARGET_CELLTYPES` (section 3).

### `02_compare_interactions.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP5_DIR` | `"step5"` | folder | Output sub-folder under `Result/plot/MERGED_DIR_NAME/`. |
| `COMPARE_WIDTH` | `600` | pixels | PNG width. |
| `COMPARE_HEIGHT` | `500` | pixels | PNG height. |
| `COMPARE_RES` | `120` | dpi | PNG resolution. |

### `03_diff_interactions.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP6_DIR` | `"step6"` | folder | Output sub-folder. |
| `DIFF_CIRCLE_WIDTH` | `1800` | pixels | Width of the combined differential circle plot. |
| `DIFF_CIRCLE_HEIGHT` | `800` | pixels | Height of the combined differential circle plot. |
| `DIFF_CIRCLE_RES` | `120` | dpi | PNG resolution. |

### `04_circle_per_dataset.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP7_DIR` | `"step7"` | folder | Output sub-folder. |
| `CIRCLE_PER_DS_WIDTH` | `700` | pixels | Width per panel. Total width = `CIRCLE_PER_DS_WIDTH * n_datasets`. |
| `CIRCLE_PER_DS_HEIGHT` | `700` | pixels | PNG height. |
| `CIRCLE_PER_DS_RES` | `120` | dpi | PNG resolution. |
| `CIRCLE_EDGE_WIDTH_COUNT` | `8` | arbitrary | `edge.width.max` for interaction count plots. |
| `CIRCLE_EDGE_WIDTH_WEIGHT` | `6` | arbitrary | `edge.width.max` for interaction strength plots. |

### `05_circle_coarse_celltypes.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP8_DIR` | `"step8"` | folder | Output sub-folder. |
| `COARSE_WIDTH` | `600` | pixels | Width per panel. Total width = `COARSE_WIDTH * n_datasets`. |
| `COARSE_HEIGHT` | `600` | pixels | PNG height. |
| `COARSE_RES` | `120` | dpi | PNG resolution. |
| `COARSE_EDGE_WIDTH` | `12` | arbitrary | `edge.width.max` for coarse circle plots. |
| `COARSE_GROUP_DEFAULT` | `"Somatic"` | group | Coarse group assigned to any fine cell type not listed in `COARSE_GROUP_MAP`. |
| `COARSE_GROUP_MAP` | named vector | mapping | Fine cell type → coarse group mapping. |
| `COARSE_GROUP_LEVELS` | `c("SSC", "Germ", "Sertoli", "Leydig", "Somatic")` | vector | Ordered levels of the coarse groups. |

To adapt this script to another dataset, update the map with your own fine and
coarse group names.

### `06_signaling_role_scatter.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP9_DIR` | `"step9"` | folder | Output sub-folder. |
| `SIGNALING_ROLE_CELL_TYPES` | `c("Sertoli", "Leydig", "SSC")` | vector | Cell types analyzed in Step 9B. They must exist in the object. |
| `SIGNALING_ROLE_WIDTH_PER_PANEL` | `5` | inches | Width per dataset panel. |
| `SIGNALING_ROLE_HEIGHT` | `5` | inches | Plot height. |
| `SIGNALING_ROLE_DPI` | `150` | dpi | Resolution. |

### `07_net_similarity.R`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `STEP10_DIR` | `"step10"` | Output sub-folder. |
| `NET_SIM_TYPE` | `"functional"` | Similarity type: `"functional"` or `"structural"`. |
| `NET_SIM_UMAP_METHOD` | `"uwot"` | UMAP method passed to `netEmbedding`. |
| `NET_SIM_SEED` | `6` | Random seed for reproducibility. |
| `NET_SIM_WIDTH` | `10` | Figure width in inches. |
| `NET_SIM_HEIGHT` | `8` | Figure height in inches. |
| `NET_SIM_DPI` | `150` | Resolution. |

### `08_ranknet_heatmap.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP11_DIR` | `"step11"` | folder | Output sub-folder. |
| `RANKNET_WIDTH` | `10` | inches | Width of the rankNet bar plots. |
| `RANKNET_HEIGHT_PER_PW` | `0.25` | inches / pathway | Dynamic height per pathway for rankNet plots. |
| `HEATMAP_HEIGHT_PER_PW` | `0.30` | inches / pathway | Dynamic height per pathway for ComplexHeatmap. |
| `HEATMAP_HEIGHT_PER_LR` | `0.18` | inches / L-R pair | Dynamic height per L-R pair for L-R heatmaps. |
| `HEATMAP_PATTERNS` | `c("outgoing", "incoming", "all")` | vector | Signaling patterns to draw with `netAnalysis_signalingRole_heatmap`. |

### `09_bubble_dysfunctional.R`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `STEP12_DIR` | `"step12"` | Output sub-folder. |
| `DYSFUNCTIONAL_SOURCES` | `c("Sertoli", "SSC")` | Source cell types for bubble plots and DEG extraction. |
| `DYSFUNCTIONAL_TARGETS` | `c("Sertoli", "SSC")` | Target cell types for bubble plots and DEG extraction. |
| `DYSFUNCTIONAL_LIGAND_LOGFC_UP` | `0.05` | `ligand.logFC` threshold for up-regulated L-R pairs. |
| `DYSFUNCTIONAL_LIGAND_LOGFC_DOWN` | `-0.05` | `ligand.logFC` threshold for down-regulated L-R pairs. |
| `DYSFUNCTIONAL_THRESH_PC` | `0.1` | `thresh.pc` for `identifyOverExpressedGenes`. |
| `DYSFUNCTIONAL_THRESH_FC` | `0.05` | `thresh.fc` for `identifyOverExpressedGenes`. |
| `DYSFUNCTIONAL_DO_FAST` | `TRUE` | `do.fast` argument; set to `FALSE` if `presto` is not installed. |

### `10_dysfunctional_viz.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP13_DIR` | `"step13"` | folder | Output sub-folder. |
| `DYSVIZ_BUBBLE_WIDTH` | `14` | inches | Bubble plot width. |
| `DYSVIZ_BUBBLE_HEIGHT` | `12` | inches | Bubble plot height. |
| `DYSVIZ_BUBBLE_DPI` | `150` | dpi | Bubble plot resolution. |
| `DYSVIZ_CHORD_SIZE` | `10` | inches | Chord diagram width and height. |
| `DYSVIZ_WORDCLOUD_WIDTH` | `10` | inches | Wordcloud width. |
| `DYSVIZ_WORDCLOUD_HEIGHT` | `8` | inches | Wordcloud height. |

This script reuses `DYSFUNCTIONAL_*` parameters from section `09` for DEG
thresholds and source/target cell types.

### `11_pathway_viz.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP14_DIR` | `"step14"` | folder | Output sub-folder. |
| `PATHWAY_VIZ_RES` | `150` | dpi | PNG resolution. |
| `PATHWAY_VIZ_WIDTH_1` | `900` | pixels | Width for a single panel. |
| `PATHWAY_VIZ_HEIGHT_1` | `900` | pixels | Height for a single panel. |
| `PATHWAY_VIZ_WIDTH_2` | `1800` | pixels | Width for two panels. |
| `PATHWAY_VIZ_HEIGHT_2` | `900` | pixels | Height for two panels. |
| `PATHWAY_VIZ_WIDTH_N` | `900` | pixels | Width per panel when 3+ datasets. |
| `PATHWAY_VIZ_HEIGHT_N` | `1800` | pixels | Height when 3+ datasets. |

### `12_gene_expression_viz.R`

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `STEP15_DIR` | `"step15"` | folder | Output sub-folder. |
| `GENE_EXPR_VIZ_WIDTH` | `14` | inches | Plot width. |
| `GENE_EXPR_VIZ_HEIGHT` | `10` | inches | Plot height. |
| `GENE_EXPR_VIZ_DPI` | `150` | dpi | Resolution. |
| `GENE_EXPR_VIZ_TYPE` | `"violin"` | Plot type passed to `plotGeneExpression`: `"violin"` or `"dot"`. |

### `13_export_objects.R`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `STEP16_DIR` | `"step16"` | Output sub-folder. |

This script also reads `STEP12_DIR` to look for the DEG-updated CellChat object.

---

## Adapting the pipeline to a new dataset

1. Update `BASE_DIR` and `PROJECT_NAME` if the project location changes.
2. Set `CONDITIONS` and `REFERENCE` to match the new metadata.
3. Update `SEURAT_CONDITION_COL` and `SEURAT_CELL_TYPE_COL` if the Seurat metadata columns have different names.
4. Define all cell-type colors in `CELLTYPE_COLORS`.
5. If you subset, set `DO_SUBSET = TRUE` and fill `TARGET_CELLTYPES`.
6. Update `COARSE_GROUP_MAP` and `COARSE_GROUP_LEVELS` for coarse circle plots.
7. Update `SIGNALING_ROLE_CELL_TYPES`, `DYSFUNCTIONAL_SOURCES` and
   `DYSFUNCTIONAL_TARGETS` with the cell types you want to focus on.
8. Adjust figure sizes (`WIDTH`, `HEIGHT`, `RES`, `DPI`) if you have many cell
   types or pathways.

---

## Notes

- Do **not** edit derived variables (`MERGED_DIR_NAME`, `MERGE_PREFIX`) unless
  you want to override the automatic naming.
- Any parameter value change in `config.R` is picked up by all scripts on the next
  run, without touching the individual R files.
