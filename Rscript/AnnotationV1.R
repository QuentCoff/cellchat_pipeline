 
# ============================================================
# 1️⃣ Setup & Load Data
# ============================================================
print("[1/5] Setup & Load Data - Chargement des librairies...")
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(openxlsx)
print("[1/5] Librairies chargées avec succès")

print("[1/5] Chargement de l'objet Seurat humanisé...")
obj <- readRDS("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/data_input/Horse_Humanized_Seurat.rds")
print(paste("[1/5] Objet chargé:", ncol(obj), "cellules," , nrow(obj), "gènes"))
obj$Horse <- obj$samples 

# Define Constants
cell_order   <- c("SSC", "Spermatocyte", "Spermatid", "Sertoli", "Leydig", "Myoid", "Fibroblast", "Endothelial")
cond_order   <- c("Descended", "Crypto", "Immuno")
group_order  <- c("Healthy", "Uni-Descended", "Uni-Crypto", "Bi-Crypto", "Immuno")

celltype_colors <- c("SSC"="#417505", "Spermatocyte"="#4A90E2", "Spermatid"="#D0021B", 
                     "Sertoli"="#F5A623", "Leydig"="#880FF3", "Myoid"="#7ED321", 
                     "Fibroblast"="#5A5A5A", "Endothelial"="#F16B54")
cond_colors  <- c("Descended" = "#D477A0", "Crypto" = "#158364", "Immuno" = "#3E8AB1")
group_colors <- c("Healthy" = "#4E79A7", "Uni-Descended" = "#F28E2B", "Uni-Crypto" = "#76B7B2", 
                  "Bi-Crypto" = "#1B7F5A", "Immuno" = "#7B4FA3")

# Metadata Setup
obj@meta.data <- obj@meta.data %>%
  mutate(
    detailed_group = factor(case_when(
      Horse %in% c("D-CV19", "D-CV25", "D-CV30", "D-CV42") ~ "Healthy",
      Horse %in% c("D-CV22", "D-CV44")                     ~ "Uni-Descended",
      Horse %in% c("C-CV22", "C-CV44")                     ~ "Uni-Crypto",
      Horse %in% c("C-CV23")                               ~ "Bi-Crypto",
      Horse %in% c("I-CV24", "I-CV39", "I-CV40")            ~ "Immuno"
    ), levels = group_order),
    condition = factor(case_when(
      grepl("^D-", Horse) ~ "Descended",
      grepl("^C-", Horse) ~ "Crypto",
      grepl("^I-", Horse) ~ "Immuno"
    ), levels = cond_order)
  )

# ============================================================
# 2️⃣ Processing & Clustering
# ============================================================
print("[2/5] Processing & Clustering - Normalisation...")
obj <- NormalizeData(obj) %>% FindVariableFeatures(nfeatures=2000)
print("[2/5] Scaling...")
obj <- ScaleData(features=VariableFeatures(.))
print("[2/5] PCA...")
obj <- RunPCA(features=VariableFeatures(.))
print("[2/5] Harmony integration (peut prendre 5-10 min)...")
obj <- RunHarmony(obj, "Horse", "pca")
print("[2/5] Clustering et UMAP...")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.3) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30)
print("[2/5] Processing terminé!")

print("[2/5] Sauvegarde des UMAP bruts...")
# Save Raw Cluster UMAPs
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/Raw_Combined.png", DimPlot(obj, group.by="seurat_clusters"), width=7, height=6)
print("[2/5] Raw_Combined.png sauvegardé")
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/Raw_Condition_Split.png", DimPlot(obj, group.by="seurat_clusters", split.by="condition"), width=15, height=5)
print("[2/5] Raw_Condition_Split.png sauvegardé")
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/Raw_Group_Split.png", DimPlot(obj, group.by="seurat_clusters", split.by="detailed_group"), width=20, height=5)
print("[2/5] Raw_Group_Split.png sauvegardé")

# ============================================================
# 3️⃣Marker-Based Annotation (Module Scoring)
# ============================================================
print("[3/5] Marker-Based Annotation (Module Scoring)...")
markers_panel <- list(SSC=c("DAZL", "DMRT1", "GFRA1"), Spermatocyte=c("SETX", "HORMAD2", "PIWIL1"),
                      Spermatid=c("CCDC168", "EFCAB3", "CREM"), Sertoli=c("ERBB4", "SYNE2", "FSHR", "SOX9"),
                      Leydig=c("LHCGR", "STAR"), Myoid=c("ACTA2", "MYH11"),
                      Fibroblast=c("DCN", "TCF21", "EBF1"), Endothelial=c("VWF", "ETS1", "FLI1"))

for (ct in names(markers_panel)) {
  print(paste("[3/5] Calcul du score pour", ct, "..."))
  obj <- AddModuleScore(obj, features = list(markers_panel[[ct]]), name = paste0(ct, "_Score"))
}
obj$cell_type <- apply(obj@meta.data[, paste0(names(markers_panel), "_Score1")], 1, 
                       function(x) names(markers_panel)[which.max(x)])
obj$cell_type <- factor(obj$cell_type, levels = cell_order)
print("[3/5] Annotation terminée!")
print(table(obj$cell_type))

# ============================================================
# 4️⃣  Reporting & Metrics
# ============================================================
print("[4/5] Reporting & Metrics - Sauvegarde des plots annotés...")
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/plot_annotated/Annotated_Combined.png", DimPlot(obj, group.by="cell_type", cols=celltype_colors), width=8, height=6)
print("[4/5] Annotated_Combined.png sauvegardé")
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/plot_annotated/Annotated_Condition_Split.png", DimPlot(obj, group.by="cell_type", split.by="condition", cols=celltype_colors), width=15, height=5)
print("[4/5] Annotated_Condition_Split.png sauvegardé")
ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/plot_annotated/Annotated_Group_Split.png", DimPlot(obj, group.by="cell_type", split.by="detailed_group", cols=celltype_colors), width=20, height=5)
print("[4/5] Annotated_Group_Split.png sauvegardé")

export_stats <- function(group_var, folder_name) {
  outdir <- paste0("C:/snRNA/Report_", folder_name)
  dir.create(outdir, showWarnings = FALSE)
  prop_table <- prop.table(table(obj$cell_type, obj@meta.data[[group_var]]), margin = 2)
  qc_stats <- obj@meta.data %>% group_by(.data[[group_var]], cell_type) %>% 
    summarise(Total_Nuclei = n(), Mean_Genes_Per_Nucleus = mean(nFeature_RNA), .groups = 'drop')
  
  wb <- createWorkbook()
  addWorksheet(wb, "Proportions"); writeData(wb, "Proportions", as.data.frame.matrix(prop_table), rowNames = TRUE)
  addWorksheet(wb, "QC_Metrics"); writeData(wb, "QC_Metrics", qc_stats)
  saveWorkbook(wb, file.path(outdir, "Metrics.xlsx"), overwrite = TRUE)
  
  p <- ggplot(as.data.frame(prop_table), aes(x=Var2, y=Freq, fill=Var1)) +
    geom_bar(stat="identity", position="stack") + scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis()
  ggsave(file.path(outdir, "Composition_Barplot.png"), p, width=7, height=6)
}

export_stats("condition", "Condition")
export_stats("detailed_group", "Group")
saveRDS(obj, "/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/final_output/FFinal_Annotated_Object.rds")
print("[4/5] Objet RDS final sauvegardé: FFinal_Annotated_Object.rds")
print("[4/5] Génération des rapports Excel et plots de composition...")




###For the metrics 
#
#
export_stats <- function(group_var, folder_name) {
  outdir <- paste0("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/Report", folder_name)
 # dir.create(outdir, showWarnings = FALSE)
  
  # 1. Calculate the Summary Table
  
  summary_stats <- obj@meta.data %>% 
    group_by(.data[[group_var]], cell_type) %>% 
    summarise(
      Total_Nuclei = n(),
      Mean_Genes_Per_Nucleus = mean(nFeature_RNA),
      Total_Genes_In_Group = sum(nFeature_RNA),
      .groups = 'drop'
    )
  
  # 2. Calculate Proportions
  prop_table <- prop.table(table(obj$cell_type, obj@meta.data[[group_var]]), margin = 2)
  prop_df <- as.data.frame(prop_table)
  colnames(prop_df) <- c("CellType", "Group", "Proportion")
  
  # 3. Save Excel Workbook (2 Sheets)
  wb <- createWorkbook()
  addWorksheet(wb, "Proportions"); writeData(wb, "Proportions", as.data.frame.matrix(prop_table), rowNames = TRUE)
  addWorksheet(wb, "Detailed_QC"); writeData(wb, "Detailed_QC", summary_stats)
  saveWorkbook(wb, file.path(outdir, paste0(folder_name, "_Metrics_Report.xlsx")), overwrite = TRUE)
  
  # --- PLOT A: Cell Type Proportions (Stacked Bar) ---
  p1 <- ggplot(prop_df, aes(x=Group, y=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position="stack") + 
    scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + ggtitle(paste(folder_name, "Proportions"))
  ggsave(file.path(outdir, "Plot_Proportions.png"), p1, width=8, height=6)
  
  # --- PLOT B: Total Nuclei Count (Grouped Bar) ---
  p2 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Total_Nuclei, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + ggtitle(paste(folder_name, "Total Nuclei"))
  ggsave(file.path(outdir, "Plot_Total_Nuclei.png"), p2, width=8, height=6)
  
  # --- PLOT C: Mean Genes Per Nucleus (Grouped Bar) ---
  # This shows the "complexity" or "quality" of the nuclei
  p3 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Mean_Genes_Per_Nucleus, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + ggtitle(paste(folder_name, "Mean Gene Detection"))
  ggsave(file.path(outdir, "Plot_Mean_Genes.png"), p3, width=8, height=6)
}

# Run the updated function
export_stats("condition", "Condition")
export_stats("detailed_group", "Group")



##
##
#

export_detailed_stats <- function(group_var, folder_name) {
  outdir <- paste0("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/Report", folder_name)
  #dir.create(outdir, showWarnings = FALSE)
  
  # 1. Calculate the core metrics per Horse first (to get the mean per sample)
  sample_stats <- obj@meta.data %>%
    group_by(Horse, .data[[group_var]], cell_type) %>%
    summarise(
      Nuclei_Count = n(),
      Mean_Genes = mean(nFeature_RNA),
      .groups = 'drop'
    )
  
  # 2. Calculate the Mean across Horses for each group
  summary_stats <- sample_stats %>%
    group_by(.data[[group_var]], cell_type) %>%
    summarise(
      Mean_Nuclei_Per_Horse = mean(Nuclei_Count),
      Total_Nuclei_All_Samples = sum(Nuclei_Count),
      Mean_Genes_Per_Cell = mean(Mean_Genes),
      Total_Genomic_Activity = sum(Mean_Genes * Nuclei_Count),
      .groups = 'drop'
    )
  
  # 3. Save Excel Report
  write.xlsx(summary_stats, file.path(outdir, paste0(folder_name, "_Statistical_Summary.xlsx")))
  
  # --- FIGURE A: Total Nuclei (The absolute count) ---
  p1 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Total_Nuclei_All_Samples, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + labs(title=paste(folder_name, "Total Nuclei (Sum)"))
  ggsave(file.path(outdir, "1_Total_Nuclei.png"), p1, width=9, height=7)
  
  # --- FIGURE B: Mean Nuclei Per Horse (The fair comparison) ---
  p2 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Mean_Nuclei_Per_Horse, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + labs(title=paste(folder_name, "Mean Nuclei per Horse"))
  ggsave(file.path(outdir, "2_Mean_Nuclei_Per_Horse.png"), p2, width=9, height=7)
  
  # --- FIGURE C: Mean Gene Complexity (Quality check) ---
  p3 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Mean_Genes_Per_Cell, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + labs(title=paste(folder_name, "Mean Genes per Nucleus"))
  ggsave(file.path(outdir, "3_Mean_Genes_Per_Cell.png"), p3, width=9, height=7)
  
  # --- FIGURE D: Total Gene Activity (Transcriptional Output) ---
  p4 <- ggplot(summary_stats, aes(x=.data[[group_var]], y=Total_Genomic_Activity, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=celltype_colors) +
    theme_minimal() + RotatedAxis() + labs(title=paste(folder_name, "Total Transcriptional Output"))
  ggsave(file.path(outdir, "4_Total_Gene_Activity.png"), p4, width=9, height=7)
}

# Execute
export_detailed_stats("condition", "Condition")
export_detailed_stats("detailed_group", "Group")
print("[4/5] Rapports générés!")

print("[5/5] Validation - Génération des plots de validation...")
print("[5/5] Sauvegarde du DotPlot de validation...")




#
#
# ============================================================
#







#
# ============================================================
# 5️ Validation
# ============================================================
dir.create("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/validation/Individual_Horses", showWarnings = FALSE)
dir.create("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/validation/Validation_Plots", showWarnings = FALSE)

# --- FIGURE A: The Biological Validation (DotPlot) ---
#
p_val <- DotPlot(obj, features = markers_panel, group.by = "cell_type") + 
  RotatedAxis() + 
  scale_color_gradient(low = "lightgrey", high = "#D0021B") +
  theme_minimal() +
  labs(title = "Marker Gene Validation by Cell Type", x = "Marker Genes", y = "Annotated Cell Type")

ggsave("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/validation/Validation_Plots/Validation_DotPlot.png", p_val, width = 12, height = 6)
print("[5/5] Validation_DotPlot.png sauvegardé")

print(paste("[5/5] Génération des UMAP individuels pour", length(unique(obj$Horse)), "échantillons..."))

# --- FIGURE B: Loop to Save One Separate UMAP for Each Horse ---
all_horses <- unique(obj$Horse)

for (h in all_horses) {
  print(paste("[5/5] Traitement de", h, "..."))
  obj_sub <- subset(obj, subset = Horse == h)
  p <- DimPlot(obj_sub, 
               group.by = "cell_type", 
               cols = celltype_colors, 
               reduction = "umap") + 
    theme_minimal() +
    labs(title = paste("Cell Composition: Horse", h))
  file_name <- paste0("/LAB-DATA/GLiCID/users/e244507t@univ-nantes.fr/Cellchat/results/annotation_humanised/validation/Individual_Horses/UMAP_", h, ".png")
  ggsave(file_name, p, width = 8, height = 7)
  rm(obj_sub)
}
print("[5/5] Tous les UMAP individuels sauvegardés!")
print("=== SCRIPT TERMINÉ AVEC SUCCÈS ===")
### Thankssss







