# Rapport d'Orthologie : Cheval -> Humain

## 1. Statistiques Globales (Mapping Local)
- **Total des gènes de l'objet Seurat (Cheval)** : 33687
- **Gènes avec une correspondance Humaine** : 15555 (46.18%)
- **Gènes perdus (sans correspondance)** : 18132 (53.82%)

### Détails des Doublons
- **Doublons côté Cheval** : 305 (Un même gène cheval pointe vers plusieurs gènes humains. Seul le premier est conservé lors du mapping.)
- **Doublons côté Humain** : 149 (Plusieurs gènes de cheval différents convergent vers le même gène humain. Cela peut entraîner une fusion des comptes de ces gènes lors de la création de l'objet Seurat.)

## 2. Vérification en direct (Ensembl via biomaRt)
*Un test a été effectué en interrogeant directement les serveurs d'Ensembl aujourd'hui pour vérifier si la base locale est toujours à jour.*

- **Total gènes cherchés dans Seurat** : 33687
- **Total gènes trouvés sur Ensembl aujourd'hui** : 15010
- **Correspondances parfaites avec le fichier local** : 15010 (100.00% des gènes trouvés)

=================================================================
RÉSULTATS DE LA VÉRIFICATION DES MARQUEURS (23 testés)
=================================================================

Marqueurs PRÉSENTS dans l'objet (20) :
 [1] "DMRT1"   "GFRA1"   "SETX"    "HORMAD2" "PIWIL1"  "EFCAB3"  "CREM"   
 [8] "ERBB4"   "SYNE2"   "FSHR"    "SOX9"    "LHCGR"   "ACTA2"   "MYH11"  
[15] "DCN"     "TCF21"   "EBF1"    "VWF"     "ETS1"    "FLI1"   

Marqueurs ABSENTS de l'objet (3) :
[1] "DAZL" "CCDC168" "STAR"   

DAZL => SSC
CCDC168 => Spermatid
STAR => Leydig


Tout les gènes marqueurs manquants ne sont pas les plus exprimé (avrage expression/ percent expression) => positif mais à voir 
=================================================================

**Conclusion** : Le mapping local est **100% identique** aux données d'Ensembl actuelles. Mais pertes d'un peu moins e 54% des gènes
