# Clustering des Pathways dans CellChat — `computeNetSimilarityPairwise`

## Vue d'ensemble

CellChat regroupe (clustérise) les voies de signalisation (pathways) selon la **similarité de leurs réseaux de communication cellule-cellule**. Deux modes existent : `functional` et `structural`.

---

## 1. Functional (défaut)

**Principe** : deux pathways sont similaires s'ils communiquent entre les **mêmes paires de types cellulaires**.

### Comment ça marche
1. Pour chaque pathway, on extrait une matrice binaire `source × target` :
   - `1` = communication présente (probabilité > 0)
   - `0` = communication absente
2. On compare deux pathways avec l'**index de Jaccard** :
   ```
   J(A, B) = |A ∩ B| / |A ∪ B|
   ```
   - `A` et `B` = ensembles de paires source→target actives
   - `J = 1.0` → identiques
   - `J = 0.0` → aucune paire en commun
3. Filtrage SNN : on ne garde que les `k` plus proches voisins de chaque pathway.

### Exemple
- NOTCH active : Sertoli→SSC, Leydig→Myoid
- CDP active : Sertoli→SSC, Leydig→Myoid
→ **Même cluster** (Jaccard élevé)

---

## 2. Structural

**Principe** : deux pathways sont similaires si leurs réseaux de communication ont la **même topologie / architecture globale**, indépendamment des identités cellulaires.

### Comment ça marche
1. Même matrice binaire `source × target` que functional
2. On calcule une **distance structurelle** entre graphes (`computeNetD_structure`)
3. Similarité = `1 − distance`
4. Même filtrage SNN que functional

### Exemple
- NOTCH : réseau très centralisé (un hub dominant)
- CDP : réseau très centralisé (un hub dominant)
→ **Même cluster** (même architecture), même si les cellules impliquées diffèrent

---

## 3. Filtrage SNN (Shared Nearest Neighbor)

Appliqué aux deux modes :
- Pour chaque pathway, on garde seulement ses **k plus proches voisins**
- `k` par défaut : `ceiling(sqrt(nombre_total_pathways))`
- Le pruning (`prune.SNN = 1/15`) élimine les liens faibles unidirectionnels

**But** : réduire le bruit, ne conserver que les vraies proximités avant le clustering final.

---

## Récap comparatif

| | Functional | Structural |
|---|---|---|
| **Compare** | Paires source→target en commun | Architecture/topologie du réseau |
| **Métrique** | Jaccard index | Distance de graphe |
| **Interprétation** | "Mêmes routes cellulaires" | "Même forme de réseau" |
| **Usage** | Défaut, plus biologique | Plus abstrait, forme globale |

---

## Référence

Fichier source du pipeline :
- `Rscript/Procedure_2/07_net_similarity.R`
- Fonction CellChat : `computeNetSimilarityPairwise()`

# STEP 13 wordcloud
En résumé
Mot gros = ligand impliqué dans beaucoup d'interactions différentielles et/ou avec des changements d'expression importants
Mot petit = ligand peu ou modérément dysrégulé
Dans le contexte du script :

wordcloud_up = ligands up-régulés dans crypto (vs healthy)
wordcloud_down = ligands down-régulés dans crypto (vs healthy)



# STEP 15
e split.by = "datasets" du script est passé à StackedVlnPlot via les ....

Conclusion : le rouge et le bleu correspondent à tes deux conditions (datasets) :

Rouge = healthy
Bleu = crypto
Chaque violin montre la distribution d'expression d'un gène de signaling (ligand ou récepteur du pathway, ex: NOTCH) dans un type cellulaire donné, séparée par condition. Tu peux comparer visuellement si l'expression est plus élevée dans healthy ou crypto pour chaque gène dans chaque type cellulaire. 

Si tu veux vérifier l'ordre exact des couleurs :

r
levels(cellchat@meta$datasets)
Feedback submitted



la barre gradué entre 0 et 3 indique quoi ?

La barre graduée de 0 à 3 représente le niveau d'expression génique normalisé (expression level).

Ce que ça signifie
Valeur	Signification
0	Pas d'expression (ou très faible) du gène dans ces cellules
1-2	Expression modérée
3	Expression élevée
C'est quoi exactement ?
Dans plotGeneExpression, CellChat extrait les données d'expression des gènes de signaling (ligands/récepteurs) depuis l'objet Seurat. L'expression est log-normalisée (typiquement log1p(counts)).

Le violin plot montre :

Forme du violin = distribution de l'expression dans chaque groupe de cellules
Couleur/fill = niveau moyen d'expression (la barre graduée 0→3)
En pratique
Quand tu compares le rouge (healthy) et le bleu (crypto) côte à côte :

Si le côté rouge est plus foncé (vers 3) → le gène est plus exprimé dans healthy
Si le côté bleu est plus foncé (vers 3) → le gène est plus exprimé dans crypto
Note : c'est une échelle relative (log-normalisée), pas des comptes bruts. Une valeur de 3 ne veut pas dire "3 transcrits", mais "expression élevée sur cette échelle normalisée".