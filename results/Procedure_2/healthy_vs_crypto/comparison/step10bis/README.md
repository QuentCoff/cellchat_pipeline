# `07_robustness.R` — analyse de robustesse du classement des pathways

## Problème adressé

`rankSimilarity` (step 10-vi) ordonne les pathways par distance entre Healthy
et Crypto dans l'embedding UMAP. Mais **UMAP est stochastique** : chaque
exécution donne un layout 2D légèrement différent → les distances et donc le
classement varient d'un run à l'autre.

La question : *quels pathways sont réellement différents entre conditions, et
pas juste par hasard de l'UMAP ?*

## Pipeline

### 1. Étape déterministe (une seule fois)

```r
computeNetSimilarityPairwise(cellchat, type = "functional")
```

Calcule la matrice de similarité Jaccard entre pathways. **Ne dépend pas du
hasard.**

### 2. Boucle stochastique (N=100 par défaut)

À chaque itération :

```r
cc_i <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
# Sans set.seed -> UMAP utilise une graine différente à chaque appel
```

Puis extrait pour chaque pathway commun aux 2 conditions :

- **distance** = distance euclidienne entre coord Healthy et Crypto dans
  l'embedding 2D
- **rank** = rang du pathway (1 = plus grande distance)

→ Table longue `(iter, pathway, distance, rank)` sauvée dans
`iterations_long_functional.csv`.

### 3. Résumé par pathway

Pour chaque pathway, calcule sur les 100 runs :

- **Distance** : mean, sd, cv, median, IQR (q25/q75)
- **Rank** : median, IQR, min/max, fréquence de top-1/3/5/10

→ Mesure la stabilité : un pathway avec faible CV et IQR de rang étroit est
**robuste**.

### 4. Tests statistiques pairés (Wilcoxon signed-rank)

Pour chaque paire (A, B) de pathways :

$$H_0 : \text{median}(d_A^{(i)} - d_B^{(i)}) = 0$$

Test apparié sur les 100 itérations (même UMAP pour A et B → contrôle la
variabilité de l'UMAP). Correction **BH** (Benjamini-Hochberg) sur les
K(K-1)/2 tests.

→ 3 matrices :

- `wilcoxon_pvalues_*.csv` : p-values brutes
- `wilcoxon_padj_BH_*.csv` : p-values ajustées
- `median_diff_*.csv` : différence médiane (signe = direction)

### 5. Score "wins / losses / ties"

Pour chaque pathway A :

- **wins** = nombre d'autres pathways que A bat significativement
  (padj < 0.05 ET median_diff > 0)
- **losses** = nombre d'autres pathways qui battent A
- **ties** = pas de différence significative

→ Classement final dans `per_pathway_summary_functional.csv` (trié par wins
décroissant).

### 6. Visualisations

- `boxplot_distance_functional.png` : distribution des distances sur les 100 runs
- `boxplot_rank_functional.png` : distribution des rangs sur les 100 runs

## À quoi ça sert

Un pathway "top du classement" en step 10 ne mérite confiance que si :

1. Son **rang médian** sur 100 runs est stable (IQR étroit)
2. Il **bat significativement** un grand nombre d'autres pathways
   (`wins` élevé)

Sinon, son classement initial peut être un artefact d'une seule réalisation
UMAP.

## Usage

```bash
# Default: type=functional, n_iter=100
Rscript 07_robustness.R

# Custom
Rscript 07_robustness.R functional 200
Rscript 07_robustness.R structural 100
```

## Fichiers de sortie

| Fichier | Contenu |
|---------|---------|
| `iterations_long_<type>.csv` | Table longue (iter, pathway, distance, rank) |
| `rank_frequency_<type>.csv` | P(rang = k) par pathway |
| `per_pathway_summary_<type>.csv` | Résumé + wins/losses/ties (trié) |
| `wilcoxon_pvalues_<type>.csv` | Matrice p-values brutes (paires de pathways) |
| `wilcoxon_padj_BH_<type>.csv` | Matrice p-values ajustées BH |
| `median_diff_<type>.csv` | Différence médiane par paire |
| `boxplot_distance_<type>.png` | Distribution des distances |
| `boxplot_rank_<type>.png` | Distribution des rangs |

Dans ce contexte, **"battre"** signifie : *avoir une distance significativement plus grande* entre Healthy et Crypto.

### Exemple concret

Prenons **CADM vs APP** sur les 100 runs UMAP :

| Métrique | CADM | APP | Différence (CADM - APP) |
|----------|------|-----|------------------------|
| Distance moyenne | 4.85 | 4.65 | +0.20 |
| Médiane des différences | | | **+0.18** |
| Wilcoxon paired test | | | p = 0.003 |

**Résultat** : `median(diff) > 0` et `padj < 0.05` → CADM **bat** APP.

### Ce que ça veut dire biologiquement

| Pathway A bat Pathway B | Signification |
|------------------------|---------------|
| **distance(A) > distance(B)** | Le réseau de communication du pathway A est **plus différent** entre Healthy et Crypto que celui du pathway B |
| **→ A est plus "dysrégulé"** | Le profil de qui parle à qui (senders/receivers) change plus entre les deux conditions pour A que pour B |

### Wins/losses en termes biologiques

| Score | Interprétation |
|-------|---------------|
| **22 wins** | Ce pathway est significativement plus dysrégulé que 22 autres pathways |
| **1 loss** | 1 seul pathway est significativement plus dysrégulé que lui |
| **0 ties** | Aucune ambiguïté — il domine clairement |

### Exemple inverse

Un pathway avec **0 wins, 23 losses, 0 ties** :
- Distance faible et stable sur les 100 runs
- → Son réseau de communication est **quasi identique** entre Healthy et Crypto
- → Il est **conservé** (pas dysrégulé)

### En résumé

**Battre = être plus dysrégulé entre conditions** (distance plus grande de façon significative sur les 100 UMAP).