=== Chargement des données ===
Dimensions totales: 86842 cellules
Conditions disponibles:

Descended    Crypto    Immuno 
    42008     26796     18038 

=== Filtrage sur Descended ===
Cellules Descended: 42008
Cell types dans Descended:

         SSC Spermatocyte    Spermatid      Sertoli       Leydig        Myoid 
        2702         8187         4989        16616         5004         2314 
  Fibroblast  Endothelial 
        1575          621 

=== Création CellChat ===
[1] "Create a CellChat object from a Seurat object"
The `data` slot in the default assay is used. The default assay is RNA 
The `meta.data` slot in the Seurat object is used as cell meta information 
Set cell identities for the new CellChat object 
The cell groups used for CellChat analysis are  SSC, Spermatocyte, Spermatid, Sertoli, Leydig, Myoid, Fibroblast, Endothelial 
Message d'avis :
Dans createCellChat(object = seurat_descended, group.by = "cell_type") :
  The 'meta$samples' is not a factor. We now force it as a factor! 

Cell groups: 8
Groupes: SSC, Spermatocyte, Spermatid, Sertoli, Leydig, Myoid, Fibroblast, Endothelial

=== Configuration CellChatDB ===

=== Pré-traitement ===
The number of highly variable ligand-receptor pairs used for signaling inference is 2265 

=== Inférence du réseau ===
triMean is used for calculating the average gene expression per cell group. 
[1] ">>> Run CellChat on sc/snRNA-seq data <<< [2026-04-13 13:29:09.547567]"
  |======================================================================| 100%
[1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2026-04-13 13:32:36.977589]"

=== Génération des visualisations ===
null device 
          1 
null device 
          1 
(env_cellchat) bash-4.4$