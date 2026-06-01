# Procedure 2 Checklist — Comparative Analysis

**Objectif :** comparer la communication cellulaire entre deux conditions via CellChat (Procedure 2 du protocole Nature Protocols).  
**Statuts :** `[ ]` = non validé / à corriger `[x]` = validé `[~]` = présent avec divergence

---

## Prérequis

- [x] **RDS créés** — `cellchat_<prefix1>.rds` et `cellchat_<prefix2>.rds` générés par `Procedure_1/01_prepare_cellchat.R`
  - Ex: `results/Procedure_1/healthy/cellchat_healthy.rds`
  - Ex: `results/Procedure_1/crypto/cellchat_crypto.rds`

---

## Steps 1–3: Load and Merge

### Step 1 — Generate CellChat object for each dataset

- [x] **Step 1** Exécuté via `Procedure_1/01_prepare_cellchat.R` pour chaque condition.  
  *OK* — les RDS existent dans `results/Procedure_1/<prefix>/`.

### Step 2 — (Optional) updateCellChat

- [x] **Step 2** Optionnel — non implémenté.  
  *OK* — CellChat >= 1.6.0 supposé.

### Step 3 — mergeCellChat

- [x] **Step 3** `mergeCellChat(object.list, add.names = names(object.list))`  
  **Code actuel :** `01_merge_cellchat.R` lignes 62–75  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/Procedure_2/01_merge_cellchat.R:62-75`  
  *OK* — merge des deux objets + export RDS + RData.

---

## Step 5: Identify altered interactions and cell populations

### Description
**Objectif :** établir si la communication cell–cell est renforcée ou affaiblie entre deux conditions en comparant :
- **Option A** : le **nombre total d'interactions** (`compareInteractions`, measure = "count")
- **Option B** : la **force totale des interactions** (`compareInteractions`, measure = "weight")

Ces barplots globaux donnent la vue d'ensemble avant d'explorer les détails (pathways, LR pairs).

### Check

- [ ] **Step 5A** `compareInteractions(cellchat, show.legend = F, group = c(1,2))`  
  Comptage total d'interactions par condition.
- [ ] **Step 5B** `compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")`  
  Force totale des interactions par condition.

### Outputs (`step5/`)
- `compare_interactions_count.png` — barplot du nombre total d'interactions par condition
- `compare_interactions_weight.png` — barplot de la force totale des interactions par condition

---

## Step 6: Differential interactions between cell populations

### Description
**Objectif :** identifier les interactions **substantiellement altérées** entre paires de types cellulaires entre deux conditions.
- **Option A** : circle plot différentiel (positif = augmenté, négatif = diminué)
- **Option B** : heatmap différentielle

Les valeurs positives = interactions plus fortes dans la 2ème condition, négatives = plus fortes dans la 1ère.

### Check

- [ ] **Step 6A-i** `netVisual_diffInteraction(cellchat, weight.scale = T)` — circle plot, count
- [ ] **Step 6A-ii** `netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")` — circle plot, strength
- [ ] **Step 6B-i** `netVisual_heatmap(cellchat)` — heatmap, count
- [ ] **Step 6B-ii** `netVisual_heatmap(cellchat, measure = "weight")` — heatmap, strength

### Outputs (`step6/`)
- `diff_circle_count.png` — circle plot différentiel (rouge = + dans condition 2, bleu = - dans condition 2), count
- `diff_circle_weight.png` — idem, force
- `diff_heatmap_count.png` — heatmap différentielle, count
- `diff_heatmap_weight.png` — heatmap différentielle, force

---

## Step 7: Circle plots across multiple datasets (normalized)

### Description
**Objectif :** comparer les réseaux de communication entre conditions avec des **circle plots normalisés** sur la même échelle.
- `getMaxWeight` calcule le max de cellules par groupe et le max d'interactions sur tous les datasets → échelle commune
- Un circle plot par dataset, affichés côte à côte

### Check

- [ ] **Step 7-i** `getMaxWeight(object.list, attribute = c("idents", "count"))` — normalisation échelle commune
- [ ] **Step 7-ii** `netVisual_circle` en boucle sur `object.list` — un plot par condition

### Outputs (`step7/`)
- `circle_per_dataset_count.png` — panel : un circle plot par condition, même échelle, count
- `circle_per_dataset_weight.png` — idem, force

---

## Step 8: Circle plots at coarse cell type level

### Description
**Objectif :** simplifier le réseau en agrégeant les types cellulaires en grandes catégories, puis comparer les interactions entre conditions.

Groupes définis :
- **Germ** : SSC, Spermatocyte, Spermatid
- **Somatic** : Sertoli, Leydig, Myoid, Fibroblast, Endothelial

### Check

- [ ] **Step 8-i** Définir `group.cellType` — Germ / Somatic
- [ ] **Step 8-ii** `mergeInteractions` + `mergeCellChat` sur les groupes coarse
- [ ] **Step 8-iii** `netVisual_circle` sur `count.merged` par dataset — échelle commune
- [ ] **Step 8-iv** `netVisual_diffInteraction(measure = "count.merged", label.edge = T)` — diff count
- [ ] **Step 8-v** `netVisual_diffInteraction(measure = "weight.merged", label.edge = T)` — diff strength

### Outputs (`step8/`)
- `circle_coarse_per_dataset_count.png` — panel Germ/Somatic par condition, même échelle
- `diff_coarse_count.png` — circle plot différentiel coarse, count
- `diff_coarse_weight.png` — circle plot différentiel coarse, force

---

## Step 9: Compare major sources and targets in 2D space

### Description
**Objectif :** identifier les types cellulaires qui changent le plus leur rôle d'envoi/réception de signaux entre conditions.
- **Option A** : scatter plot de tous les types cellulaires (outgoing vs incoming), sur échelle commune
- **Option B** : scatter plot de changements pour **un type cellulaire spécifique** (optionnel, arg 3)

### Check

- [ ] **Step 9A-i** `num.link` + `weight.MinMax` — normalisation échelle commune
- [ ] **Step 9A-ii** `netAnalysis_signalingRole_scatter` en boucle — un panel par condition
- [ ] **Step 9B** `netAnalysis_signalingChanges_scatter(cellchat, idents.use = "...")` — optionnel

### Outputs (`step9/`)
- `signaling_role_scatter.png` — scatter outgoing vs incoming par condition, même échelle
- `signaling_changes_<cell_type>.png` — changements de signalisation pour un type cellulaire spécifique (optionnel)

---

## Step 10: Signaling pathway similarity (functional/structural)

### Description
**Objectif :** identifier les pathways avec les plus grandes différences entre conditions via similarité fonctionnelle ou structurelle + manifold learning.
- `computeNetSimilarityPairwise` → `netEmbedding` → `netClustering` → visualisation 2D
- `rankSimilarity` classe les pathways par distance entre conditions

### Check

- [ ] **Step 10-i** `computeNetSimilarityPairwise(cellchat, type = "functional")`
- [ ] **Step 10-ii** `netEmbedding(cellchat, type = "functional")`
- [ ] **Step 10-iii** `netClustering(cellchat, type = "functional")`
- [ ] **Step 10-iv** `netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)`
  ⚠️ `netVisual_embedding` est buggée sur merged objects (cherche la clé `"single"` alors que les données sont `"1-2"`) — utiliser la version `Pairwise`.
- [ ] **Step 10-v** `netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2)` — optionnel
  ⚠️ Même raison : `netVisual_embeddingZoomIn` hardcode `"single"` → inutilisable sur merged object.
- [ ] **Step 10-vi** `rankSimilarity(cellchat, slot.name = "netP", type = "functional", comparison2 = c(1,2))`

### Outputs (`step10/`)
- `embedding_functional.png` — visualisation 2D des pathways (grouped by similarity)
- `embedding_zoomin_functional.png` — zoom par groupe de pathways
- `rank_similarity_functional.png` — classement des pathways par distance entre conditions

---

## Step 11: Identify altered signaling with distinct interaction strength

### Description
**Objectif :** comparer le flux d'information (force des interactions) de chaque voie de signalisation ou paire L–R entre deux conditions. CellChat identifie les pathways qui s'éteignent, diminuent, s'activent ou augmentent.

### Option A — rankNet (overall information flow)

- [ ] **Step 11A-i** `rankNet(slot.name = "netP", mode = "comparison", measure = "weight", stacked = T, do.stat = F)`
  Barplot empilé des pathways, sans test statistique.
- [ ] **Step 11A-ii** `rankNet(slot.name = "netP", mode = "comparison", measure = "weight", stacked = T, do.stat = T)`
  Barplot empilé des pathways, **test de Wilcoxon apparié**.
- [ ] **Step 11A-iii** `rankNet(slot.name = "net", mode = "comparison", measure = "weight", stacked = T, do.stat = T)`
  Barplot empilé des **paires L–R**, test de Wilcoxon apparié.
- [ ] **Step 11A-iv** `rankNet(slot.name = "netP", mode = "comparison", measure = "weight", stacked = F, do.stat = F)`
  Barplot groupé des pathways, sans test statistique.

### Outputs (`step11/` Option A)
- `ranknet_pathways_stacked.png`
- `ranknet_pathways_stacked_stat.png`
- `ranknet_lr_stacked_stat.png`
- `ranknet_pathways_grouped.png`

### Option B — ComplexHeatmap (outgoing / incoming / all)

- [ ] **Step 11B-i/ii** `netAnalysis_signalingRole_heatmap(pattern = "outgoing")` — côte à côte pour les deux datasets
- [ ] **Step 11B-iii/iv** `netAnalysis_signalingRole_heatmap(pattern = "incoming")` — côte à côte pour les deux datasets
- [ ] **Step 11B-v/vi** `netAnalysis_signalingRole_heatmap(pattern = "all")` — côte à côte pour les deux datasets (optionnel)

### Outputs (`step11/` Option B)
- `heatmap_outgoing.png`
- `heatmap_incoming.png`
- `heatmap_all.png` (optionnel)

---

## Step 12: Identify dysfunctional signaling

### Description
**Objectif :** identifier les voies de signalisation dysfonctionnelles en comparant les probabilités de communication (Option A) ou par analyse d'expression différentielle (Option B).

### Option A — Bubble plots (communication probabilities)

- [ ] **Step 12A-i** `netVisual_bubble` — toutes les communications L–R entre les deux datasets
- [ ] **Step 12A-ii** `netVisual_bubble(max.dataset = 2, remove.isolate = T)` — L–R **up-régulés** dans le dataset 2
- [ ] **Step 12A-iii** `netVisual_bubble(max.dataset = 1, remove.isolate = T)` — L–R **down-régulés** dans le dataset 2

### Outputs (`step12/` Option A)
- `bubble_all.png`
- `bubble_up.png`
- `bubble_down.png`

### Option B — Differential expression analysis

- [ ] **Step 12B-i/ii** `identifyOverExpressedGenes` — analyse DE entre conditions (`do.fast = TRUE` avec presto)
- [ ] **Step 12B-iii (opt)** `identifyOverExpressedGenes(group.DE.combined = TRUE)` — DE en ignorant les groupes cellulaires
- [ ] **Step 12B-iv** `netMappingDEG` — mappe les résultats DE sur les communications
- [ ] **Step 12B-v** `subsetCommunication(datasets = pos.dataset, ligand.logFC = 0.05)` — L–R **up** dans le dataset 2
- [ ] **Step 12B-vi** `subsetCommunication(datasets = dataset1, ligand.logFC = -0.05)` — L–R **down** dans le dataset 2
- [ ] **Step 12B-vii (opt)** `extractGeneSubsetFromPair` — gènes individuels up/down
- [ ] **Step 12B-viii (opt)** `findEnrichedSignaling` — signaux enrichis selon gènes custom

### Outputs (`step12/` Option B)
- `net_up.csv`, `net_down.csv` — paires L–R up/down
- `gene_up.csv`, `gene_down.csv` — gènes individuels (optionnel)
- `enriched_signaling_up.csv` — signaux enrichis (optionnel)
- `cellchat_deg.RData` — objet CellChat avec résultats DE

## Step 13: Visualize up/down regulated signaling events

### Description
**Objectif :** visualiser les interactions L–R up-régulées (Step 12B-v) et down-régulées (Step 12B-vi) avec bubble plot, chord diagram ou wordcloud.

### Option A — Bubble plot

- [ ] **Step 13A-i/ii** `netVisual_bubble(pairLR.use = net.up)` — up-regulated signaling in dataset 2
- [ ] **Step 13A-iii/iv** `netVisual_bubble(pairLR.use = net.down)` — down-regulated signaling in dataset 2
- [ ] **Step 13A-v** `gg1 + gg2` — combine both bubble plots (optionnel)

### Option B — Chord diagram

- [ ] **Step 13B-i** `netVisual_chord_gene(object.list[[2]], net = net.up)` — up-regulated in dataset 2
- [ ] **Step 13B-ii** `netVisual_chord_gene(object.list[[1]], net = net.down)` — down-regulated in dataset 2 (visualisé dans dataset 1)

### Option C — Wordcloud plot

- [ ] **Step 13C-i** `computeEnrichmentScore(net.up, species = 'human')` — ligands enrichis dataset 2
- [ ] **Step 13C-ii** `computeEnrichmentScore(net.down, species = 'human')` — ligands enrichis dataset 1

### Outputs (`step13/`)
- `bubble_up_dysfunctional.png`, `bubble_down_dysfunctional.png`
- `chord_up.pdf`, `chord_down.pdf`
- `wordcloud_up.pdf`, `wordcloud_down.pdf`

## Step 14: Visually compare inferred cell-cell communication networks

### Description
**Objectif :** comparer visuellement les réseaux de communication pour des **pathways spécifiques** entre les deux conditions via circle plots ou heatmaps.

### Option A — Circle plots

- [ ] **Step 14A** `getMaxWeight` + `netVisual_aggregate(layout = "circle")` pour chaque pathway
  - Échelle commune (`edge.weight.max`) pour comparer les deux datasets côte à côte

### Option B — Heat map plots

- [ ] **Step 14B** `netVisual_heatmap` pour chaque pathway, côte à côte avec `ComplexHeatmap::draw`

### Pathways par défaut
`CXCL`, `BMP`, `WNT`, `VEGF`, `NOTCH` — modifiables via le 3ème argument du script.

### Outputs (`step14/`)
- `circle_<pathway>.png` — circle plots comparatifs (1 page = 2 panels)
- `heatmap_<pathway>.png` — heatmaps comparatifs (1 page = 2 panels)

## Step 15: Gene expression distribution of signaling genes

### Description
**Objectif :** visualiser la distribution d'expression des gènes impliqués dans chaque pathway de signaling à l'aide de `plotGeneExpression` (wrapper Seurat).

### Option A — Violin plot

- [ ] **Step 15A** `plotGeneExpression(type = "violin", split.by = "datasets")` — comparaison de l'expression des L-R genes entre conditions

### Option B — Dot plot (bonus)

- [ ] **Step 15B** `plotGeneExpression(type = "dot", split.by = "datasets")` — dot plot alternatif pour la même information

### Optional — Reorder datasets

- [ ] **Step 15-i (opt)** `cellchat@meta$datasets <- factor(..., levels = c(prefix1, prefix2))` — contrôle l'ordre d'apparition dans les plots

### Pathways par défaut
`CXCL`, `BMP`, `WNT`, `VEGF`, `NOTCH` — modifiables via le 3ème argument du script.

### Outputs (`step15/`)
- `violin_<pathway>.png` — violin plots comparatifs par dataset
- `dot_<pathway>.png` — dot plots comparatifs par dataset

## Step 16: Export CellChat objects

### Description
**Objectif :** exporter l'objet `cellchat` mergé et la liste `object.list` sous forme de fichiers `.RData` portables.

- [ ] **Step 16A** `save(object.list, ...)` — exporte la liste des deux objets séparés
- [ ] **Step 16B** `save(cellchat, ...)` — exporte l'objet mergé (avec DEG si disponible)

### Outputs (`step16/`)
- `cellchat_object.list_<prefix1>_<prefix2>.RData` — liste des deux objets CellChat
- `cellchat_merged_<prefix1>_<prefix2>.RData` — objet mergé (DEG inclus si dispo)

---

## Prochaines étapes (à implémenter)

- [ ] **Step 17+** River plots, sankey, autres visualisations comparatives avancées

---

**Résumé :** Steps 1–16 implémentés.
**Scripts :**
- `Rscript/Procedure_2/08_ranknet_heatmap.R` (Step 11)
- `Rscript/Procedure_2/09_bubble_dysfunctional.R` (Step 12)
- `Rscript/Procedure_2/10_dysfunctional_viz.R` (Step 13)
- `Rscript/Procedure_2/11_pathway_viz.R` (Step 14)
- `Rscript/Procedure_2/12_gene_expression_viz.R` (Step 15)
- `Rscript/Procedure_2/13_export_objects.R` (Step 16)
