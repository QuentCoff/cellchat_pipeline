# Résultats CellChat — Communications Intercellulaires (Testicule de cheval)

## Contexte

- **Espèce** : Cheval (Equus caballus), gènes convertis en orthologues humains via Ensembl
- **Orthologie** : 16 577 gènes mappés sur 33 687 (normal — gènes non codants, pseudogènes et gènes espèce-spécifiques n'ont pas d'orthologue)
- **Base CellChat** : CellChatDB.human (ligands, récepteurs, co-facteurs humains)
- **8 types cellulaires** annotés par module scoring :

| Type cellulaire | Marqueurs | Couleur |
|-----------------|-----------|---------|
| SSC | DAZL, DMRT1, GFRA1 | Vert foncé `#417505` |
| Spermatocyte | SETX, HORMAD2, PIWIL1 | Bleu `#4A90E2` |
| Spermatid | CCDC168, EFCAB3, CREM | Rouge `#D0021B` |
| Sertoli | ERBB4, SYNE2, FSHR, SOX9 | Orange `#F5A623` |
| Leydig | LHCGR, STAR | Violet `#880FF3` |
| Myoid | ACTA2, MYH11 | Vert clair `#7ED321` |
| Fibroblast | DCN, TCF21, EBF1 | Gris `#5A5A5A` |
| Endothelial | VWF, ETS1, FLI1 | Corail `#F16B54` |

- **3 conditions** :
  - **Descended** (42 008 cellules) : testicules descendus
  - **Crypto** (26 796 cellules) : testicules cryptorchides
  - **Immuno** (18 038 cellules) : testicules immunocastrés

---

## Analyse globale (`cellchat_test.R`)

### 1. `comparison_nb_interactions.png` — Descended vs Crypto

**Ce qu'on voit** : Deux barplots côte à côte.
- **Barplot gauche** : nombre total d'interactions inférées
  - Descended (rose/corail) : **821 interactions**
  - Crypto (turquoise) : **536 interactions**
- **Barplot droit** : force totale des interactions (somme des probabilités)
  - Descended : **22.478**
  - Crypto : **17.333**

**Interprétation** : Les testicules descendus (sains) ont ~53% plus d'interactions et ~30% plus de force de signalisation que les cryptorchides. La cryptorchidie perturbe globalement les communications intercellulaires.

---

### 2. `comparison_desc_vs_immuno.png` — Descended vs Immuno

**Ce qu'on voit** : Même format que ci-dessus mais comparant Descended à Immuno.

**Interprétation** : Permet de voir si l'immunocastration impacte les communications différemment de la cryptorchidie.

---

### 3. `diff_interactions_network.png` — Réseau différentiel (Descended vs Crypto)

**Ce qu'on voit** : Deux réseaux en cercle.
- Chaque **nœud** = un type cellulaire (SSC, Sertoli, Leydig, etc.)
- Chaque **lien** = la différence d'interactions entre les deux conditions
- **Lien rouge** = plus d'interactions dans Descended que dans Crypto
- **Lien bleu** = plus d'interactions dans Crypto que dans Descended
- **Épaisseur du lien** = amplitude de la différence
- **Réseau gauche** : différentiel en nombre d'interactions
- **Réseau droit** : différentiel en force (poids) des interactions

**Interprétation** : Permet d'identifier quelles paires de types cellulaires sont les plus affectées par la cryptorchidie. Par exemple, si le lien Sertoli→SSC est rouge et épais, cela signifie que la communication Sertoli→SSC est fortement réduite chez les cryptorchides.

---

### 4. `diff_interactions_heatmap.png` — Heatmap différentielle

**Ce qu'on voit** : Deux heatmaps (matrice carrée).
- **Lignes** = types cellulaires émetteurs (source du signal)
- **Colonnes** = types cellulaires récepteurs (cible du signal)
- **Couleur rouge** = plus d'interactions dans Descended
- **Couleur bleu** = plus d'interactions dans Crypto
- **Intensité** = amplitude de la différence
- **Heatmap gauche** : différentiel en nombre
- **Heatmap droite** : différentiel en force

**Interprétation** : Version matricielle du réseau différentiel, plus facile à lire pour quantifier précisément les différences entre émetteurs et récepteurs.

---

### 5. `network_count_descended.png` / `crypto` / `immuno` — Réseaux individuels (nombre)

**Ce qu'on voit** : Un réseau en cercle par condition.
- **Nœuds** = types cellulaires, taille proportionnelle au nombre de cellules
- **Liens** = interactions entre types cellulaires
- **Épaisseur** = nombre d'interactions entre la paire
- **Couleur des liens** = correspond au type cellulaire source

**Interprétation** : Visualise l'architecture globale des communications pour chaque condition séparément. On peut comparer visuellement la "densité" de communication entre conditions.

---

### 6. `network_weight_descended.png` / `crypto` / `immuno` — Réseaux individuels (force)

**Ce qu'on voit** : Même format que les réseaux précédents mais l'épaisseur des liens représente la **force** (probabilité) des interactions au lieu du nombre brut.

**Interprétation** : Distingue les interactions rares mais fortes (biologiquement importantes) des interactions fréquentes mais faibles.

---

## Fichiers RDS sauvegardés

| Fichier | Contenu |
|---------|---------|
| `cellchat_descended.rds` | Objet CellChat complet pour Descended — contient toutes les interactions, voies de signalisation, centralité |
| `cellchat_crypto.rds` | Idem pour Crypto |
| `cellchat_immuno.rds` | Idem pour Immuno |
| `cellchat_merged_desc_vs_crypto.rds` | Objet mergé pour analyses comparatives Descended vs Crypto |
| `cellchat_merged_desc_vs_immuno.rds` | Objet mergé pour analyses comparatives Descended vs Immuno |

Ces fichiers `.rds` permettent de reprendre l'analyse dans R sans tout recalculer (ex: explorer des voies spécifiques, faire des bubble plots, etc.).

---

## Analyse Descended détaillée (`D_analysis.R`)

### Échantillons

| Échantillon | Groupe |
|-------------|--------|
| D-CV19 | Sain |
| D-CV25 | Sain |
| D-CV30 | Sain |
| D-CV42 | Sain |
| D-CV22 | Uni-Descended |
| D-CV44 | Uni-Descended |

### Résultats dans `results/descended_analysis/`

#### Visualisations individuelles
- `network_count_D-CVXX.png` : réseau des interactions (nombre) pour chaque échantillon
- `network_weight_D-CVXX.png` : réseau des interactions (force) pour chaque échantillon

#### Comparaisons sains vs sains (`sains_vs_sains/`)
- `compare_D-CVXX_vs_D-CVYY.png` : barplots nb + force
- `diff_net_D-CVXX_vs_D-CVYY.png` : réseau différentiel
- `diff_heatmap_D-CVXX_vs_D-CVYY.png` : heatmap différentielle
- Permet de voir la **variabilité inter-individuelle** entre échantillons sains

#### Comparaisons Uni vs sains (`uni_vs_sains/`)
- Même format mais compare chaque Uni-Descended (D-CV22, D-CV44) à chaque sain
- Labels : `D-CV22_Uni` vs `D-CV19_Sain`
- Permet de voir si les testicules descendants d'un animal uni-cryptorchide se comportent comme un testicule sain ou non

#### Résumé global
- `summary_interactions.csv` : tableau avec nb cellules, nb interactions et force par échantillon
- `summary_barplot_interactions.png` : barplot comparatif coloré :
  - **Bleu** (`#4E79A7`) = échantillons sains
  - **Orange** (`#F28E2B`) = échantillons Uni-Descended
- `summary_barplot_strength.png` : idem mais pour la force des interactions
