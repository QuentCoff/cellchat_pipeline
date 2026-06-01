# Procedure 1 Checklist — Data input, Preprocessing & Inference

**Objectif :** vérifier que `Rscript/interaction/descended.R` suit exactement la Procedure 1 du protocole CellChat (Nature Protocols).  
**Statuts :** `[ ]` = non validé / à corriger `[x]` = validé `[~]` = présent avec divergence

---

## A. Data input and preprocessing

### Step 1 — Option B (Seurat object) → SKIPPED

- [x] **1B** Non applicable — passage direct à l'Option B du Step 2.

### Step 2 — Option B (createCellChat from Seurat)

- [x] **2B** `createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")`  
  **Code actuel :** `cellchat <- createCellChat(object = seurat_descended, group.by = "cell_type", assay = "RNA")`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:39`  
  *OK* — Option B directe depuis Seurat, `group.by = "cell_type"`.

### Step 3 — (Optional) addMeta / setIdent

- [x] **Step 3** Optionnel — non présent.  
  *OK* — le `meta` est déjà passé lors du `createCellChat`.

### Step 4 — Select L–R database

- [x] **Step 4** `CellChatDB <- CellChatDB.human` puis `showDatabaseCategory()` + `dplyr::glimpse()`  
  **Code actuel :**  
  `CellChatDB <- CellChatDB.human`  
  `showDatabaseCategory(CellChatDB)`  
  `dplyr::glimpse(CellChatDB$interaction)`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:59-63`  
  *OK* — commandes d'inspection ajoutées.

### Step 5 — Select subset of DB

- [x] **5C** `CellChatDB.use <- CellChatDB` (all DB)  
  **Code actuel :** `CellChatDB <- CellChatDB.human` puis assignation directe  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:67-68`  
  *OK* — correspond à l'Option C (all DB).

### Step 6 — Set DB and subsetData

- [x] **Step 6** `cellchat@DB <- CellChatDB.use` + `cellchat <- subsetData(cellchat)`  
  **Code actuel :** lignes 68 et 75  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:68`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:75`  
  *OK*.

### Step 7 — identifyOverExpressedGenes

- [x] **Step 7** `future::plan("multisession", workers = 4)`  
  **Code actuel :** `future::plan("multisession", workers = 4)`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:76`  
  *OK* — ajouté.

- [x] **Step 7** `cellchat <- identifyOverExpressedGenes(cellchat)`  
  **Code actuel :** `cellchat <- identifyOverExpressedGenes(cellchat)`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:78`  
  *OK* — `presto` installé, `do.fast = FALSE` retiré.

### Step 8 — identifyOverExpressedInteractions

- [x] **Step 8** `cellchat <- identifyOverExpressedInteractions(cellchat)`  
  **Code actuel :** ligne 77  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:77`  
  *OK*.

### Step 9 — (Optional) smoothData

- [x] **Step 9** Optionnel — non présent.  
  *OK*.

---

## B. Inference of cell–cell communication networks

### Step 10 — (Optional) computeAveExpr

- [x] **Step 10** Optionnel — non présent.  
  *OK*.

### Step 11 — computeCommunProb (L–R pair level)

- [~] **Step 11** `computeCommunProb(cellchat, type = "triMean", trim = NULL, raw.use = TRUE)`  
  **Code actuel :** `computeCommunProb(cellchat, type = "triMean", trim = NULL, raw.use = TRUE, population.size = TRUE)`  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:86`  
  *OK* — `type` et `trim` ajoutés explicitement.  
  *Note :* `population.size = TRUE` ajouté (paramètre valide pour compenser le déséquilibre des tailles de clusters).

### Step 12 — filterCommunication

- [x] **Step 12** `filterCommunication(cellchat, min.cells = 10)`  
  **Code actuel :** ligne 85  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:85`  
  *OK*.

### Step 13 — computeCommunProbPathway

- [x] **Step 13** `computeCommunProbPathway(cellchat)`  
  **Code actuel :** ligne 88  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:88`  
  *OK*.

### Step 14 — aggregateNet

- [x] **14A** `aggregateNet(cellchat)` (across all groups)  
  **Code actuel :** ligne 91  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:91`  
  *OK*.

### Step 15 — saveRDS

- [x] **Step 15** `saveRDS(cellchat, file = "...")`  
  **Code actuel :** lignes 409–410  
  `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/interaction/descended.R:409-410`  
  *OK*.

---

## Résumé des divergences à corriger

| # | Élément | Statut |
|---|---------|--------|
| — | `future::plan("multisession", workers = 4)` | ✅ Corrigé |
| — | `type = "triMean"`, `trim = NULL` explicités | ✅ Corrigé |
| — | `do.fast = FALSE` retiré (`presto` installé) | ✅ Corrigé |
| — | `population.size = TRUE` | ⏳ **À valider** — garder ou retirer ? |

---

**Prochaine étape :** valider `population.size`, puis isolation des blocs en fichiers séparés.
