# Command preparation: 
## Lancement de la node de calucle: 
srun --cluster=nautilus -N1 --qos=quick --cpus-per-task=2 --mem=32G --time=2:00:00 --pty bash

## Activation de l'environnement:
micromamba activate env_cellchat

## Chargement des modules:
module load gcc/13.1.0 cmake/3.26.4

## Lancement de script R:
Rscript /LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/


# Commandes utiles:
## Transfert de fichier local => Glicid:
scp /home/quentin/Documents/ETUDE/stage/Stage_2026/dataset/AnnotationV1.R glicid:/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/Rscript/



