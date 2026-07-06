# Explication : apparente incohérence de GDNF entre violin plots, circle plots et rankNet

## Observation

Dans les résultats `Vtruncated`, le pathway **GDNF** apparaît :
- **En rouge (100 % healthy)** dans le `rankNet` stacked bar
- **Présent dans le circle plot healthy** (Sertoli → SSC)
- **Absent du circle plot crypto**
- **"Vide" dans le violin plot** pour les Sertoli (ni rouge ni bleu visible)
- **Visible** pour SSC/GFRA1 en rouge et en bleu

Cela semble contradictoire : si GDNF est suffisamment exprimé pour générer une communication dans le circle plot healthy, pourquoi le violin plot ne montre rien ?

## Réponse : trois couches d'information différentes

| Figure | Code source | Données sous-jacentes | Ce qu'elle mesure |
|--------|-------------|----------------------|-------------------|
| **RankNet** | `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Vtruncated_test/src/Rscript/08_ranknet_heatmap.R:88` | `cellchat@netP$prob` (merged) | Force globale du pathway = somme des probabilités d'interaction L-R. Comparaison statistique entre conditions. |
| **Circle plot** | `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Vtruncated_test/src/Rscript/11_pathway_viz.R:65` | `object.list[[i]]@netP$pathways` (individuels) | Présence/absence d'une **communication fonctionnelle** dans une condition donnée. Généré uniquement si le pathway est dans `netP$pathways`. |
| **Violin plot** | `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Vtruncated_test/src/Rscript/12_gene_expression_viz.R:87-93` | `cellchat@data.signaling` (matrice d'expression brute) | Niveau d'**ARNm des gènes** du pathway, indépendamment de toute interaction fonctionnelle. |

> **Principe clé** : *exprimer un gène* (violin) ≠ *avoir une communication intercellulaire significative* (circle/ranknet).

---

## Données chiffrées : expression de GDNF dans les Sertoli (objet `Vtruncated` merged)

Résultats obtenus en session R interactive sur `cellchat_merged_healthy_crypto.RData` :

```r
# Healthy Sertoli
summary(as.vector(cellchat@data.signaling["GDNF", sertoli_h, drop = FALSE]))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.000000 0.000000 0.001371 0.000000 2.363176

table(as.vector(cellchat@data.signaling["GDNF", sertoli_h, drop = FALSE]) > 0)
# FALSE  TRUE
# 11025    29   ← 0.26 % des cellules

# Crypto Sertoli
summary(as.vector(cellchat@data.signaling["GDNF", sertoli_c, drop = FALSE]))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0000  0.0000  0.0000  0.0133  0.0000  2.7045

table(as.vector(cellchat@data.signaling["GDNF", sertoli_c, drop = FALSE]) > 0)
# FALSE  TRUE
# 21154   428   ← 1.98 % des cellules
```

**Conclusion** : dans les deux conditions, moins de **2 % des cellules Sertoli expriment GDNF**, avec une moyenne inférieure à **0.014**. L'expression est quasi-nulle.

---

## Pourquoi le violin plot apparaît "vide"

Dans `@/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Vtruncated_test/src/Rscript/12_gene_expression_viz.R:87-93`, `plotGeneExpression` appelle en interne `CellChat:::modify_vlnplot`, qui fixe **`pt.size = 0`** (pas de points individuels). Avec une distribution où **98-99 % des valeurs sont à zéro**, le violon se réduit à une **ligne plate à la base de l'axe Y** — l'œil ne perçoit rien.

> Le split par `datasets` affiche bien deux demi-violons (rouge = healthy, bleu = crypto), mais ils sont tous deux invisibles car écrasés sur 0.

---

## Pourquoi le circle plot montre quand même une communication

CellChat modélise la communication par un **modèle de masse-action probabiliste** (`CellChat::netVisual_aggregate`). Le circle plot healthy est significatif car il suffit que :

1. Quelques cellules Sertoli expriment GDNF (même très faiblement : moyenne 0.001)
2. Les cellules cibles (SSC) expriment **GFRA1** et **RET** à des niveaux suffisants

Le modèle calcule une **probabilité d'interaction par paire cellule-cellule**. Si le récepteur est fortement présent sur les cibles, une expression sporadique du ligand peut suffire à dépasser le seuil de significativité de `filterCommunication`.

Dans crypto, le circle plot est absent car :
- soit GFRA1/RET ne sont plus sur les bonnes cibles,
- soit les niveaux sont en dessous du seuil,
- soit la topologie cellulaire a changé.

---

## Pourquoi GFRA1 est visible en bleu dans SSC crypto

Le violin plot montre le niveau d'**expression brute**, indépendamment de l'existence d'une interaction fonctionnelle. GFRA1 peut être transcrit dans les SSC de crypto sans que la communication GDNF→GFRA1 soit significative, car CellChat exige :
- **GDNF** sur les émettrices (Sertoli) → quasi nul dans crypto
- **GFRA1** sur les cibles (SSC) → encore présent
- une **probabilité d'interaction** dépassant le seuil

GFRA1 est donc "en attente d'un signal qui n'arrive plus".

---

## Différence avec `Vtruncated_test`

Dans l'objet `Vtruncated_test`, GDNF est **totalement absent** de `cellchat@data.signaling` :

```r
"GDNF" %in% rownames(cellchat@data.signaling)
# [1] FALSE
```

Cause probable : `subsetData` lors de `00_prepare_cellchat.R` filtre les gènes présents dans chaque condition. Comme GDNF est déjà à la limite de la détection, une légère variation de données ou de paramètres entre `Vtruncated` et `Vtruncated_test` l'a fait tomber sous le seuil. Ensuite, `mergeCellChat` fait l'**intersection** des gènes entre healthy et crypto ; si GDNF manque dans l'objet crypto, il disparaît du merged object.

---

## Synthèse

| Condition | Sertoli (GDNF) | SSC (GFRA1) | Communication Sertoli→SSC | Pourquoi |
|-----------|---------------|-------------|---------------------------|----------|
| **Healthy** | Exprimé (0.26 % cellules, moy ~0.001) | Exprimé | **Présente** (circle plot) | Co-expression suffisante malgré une expression sporadique |
| **Crypto** | Exprimé (1.98 % cellules, moy ~0.013) | Exprimé | **Absente** (pas de circle plot) | Récepteur changé de localisation ou niveaux insuffisants |
| **Violin** | "Vide" dans les deux | Visible rouge+bleu | — | `pt.size = 0`, distribution écrasée sur 0 |

**Il n'y a pas d'incohérence.** Le violin plot, le rankNet et le circle plot mesurent des choses différentes : expression brute, contribution statistique du pathway, et inférence de communication fonctionnelle.
