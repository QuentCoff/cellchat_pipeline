# Installation de l'environnement CellChat
## Une seule commande (depuis une node de calcul ou le login)
bash /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/setup.sh

Pour choisir un autre nom d'environnement :
bash /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/setup.sh mon_env_cellchat

# Command preparation:
## Lancement de la node de calucle:
srun --cluster=nautilus -N1 --qos=quick --cpus-per-task=2 --mem=100G --time=2:00:00 --pty bash

## Activation de l'environnement:
micromamba activate env_cellchat

## Chargement des modules:
module load gcc/13.1.0 cmake/3.26.4

## Lancement d'un script R (exemple):
Rscript /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/ALL_clusters/src/Rscript/00_prepare_cellchat.R


# Commandes utiles:
## Transfert de fichier local => Glicid:
scp /home/quentin/Documents/ETUDE/stage/Stage_2026/dataset/AnnotationV1.R glicid:/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/



