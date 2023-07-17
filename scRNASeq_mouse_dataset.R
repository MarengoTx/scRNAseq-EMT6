options(stringsAsFactors = F)

setwd("/Path/")


#Libraries:
library(SingleCellExperiment)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(ggExtra)
library(Matrix)
library(patchwork)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
library(scales)
library(dplyr)
library(cowplot)
library(circlize)
library(clusterProfiler)
library(enrichplot)
library(ggridges)
library(Rmagic)
library(ggpmisc)
library(glmnet)
library(foreach)
library(ggrepel)
library(EnvStats)
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)
library(ggplot.multistats)
library(ggpubr)
library(survminer)
library(harmony)
library(knitr)
library(kableExtra)
library(TILPRED)
library(ProjecTILs)

#Color schemes
c25 <- c("dodgerblue2",
         "#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2",
         "#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")


# load data ####
load("scRNASeq_mouse_2_samples_w_correct_annotation_TRBV_07192022.RData")  # < TRBV data + degs here

# read in raw files ####
sc_ctr <- Read10X(data.dir = "vehicle/filtered_feature_bc_matrix/")
sc_ctr <- CreateSeuratObject(counts = sc_ctr, project = "CTR", min.cells = 3, min.features = 200)
sc_T <- Read10X(data.dir = "Treated/filtered_feature_bc_matrix/")
sc_T <- CreateSeuratObject(counts = sc_T, project = "Treated", min.cells = 3, min.features = 200)


# prep data for integration
list_all <- merge(sc_ctr, sc_T)

list_all[["percent.mt"]] <- PercentageFeatureSet(list_all, pattern = "^mt-")
head(list_all)
VlnPlot(list_all, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
        cols = c("blue", "red"),
        pt.size = 0)

list_all <- subset(list_all, cells = which(list_all@meta.data$nFeature_RNA > 400 & list_all@meta.data$nFeature_RNA < 3500 & list_all@meta.data$percent.mt < 5)) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 12000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 40, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) 

head(list_all)
FeatureScatter(list_all, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA", 
               col=c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                     '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                     '12' = c25[13],'13' = c25[14]), 
               pt.size = 0.25) +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap( ~ list_all$orig.ident)

table(list_all$orig.ident, as.character(sapply(rownames(list_all@meta.data), function(xx) strsplit(xx, "-")[[1]][2]))) %>% 
  kable() %>% kable_styling()
DimPlot(list_all, 
        reduction = "umap", 
        order = T,
        split.by = "orig.ident",
        cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                 '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                 '12' = c25[13],'13' = c25[14]))

FeaturePlot(sc_harmony,
            features = c("Cd8b1", "Cd4", "Klrg1"), 
            pt.size = 0.25, 
            order = T, 
            split.by = "orig.ident",
            min.cutoff = "q75")

# run harmony to merge data
head(list_all@meta.data)
table(list_all$orig.ident)
sc_harmony <- list_all %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
DimPlot(sc_harmony, reduction = "harmony")

# UMAP and Clustering ####
sc_harmony <- sc_harmony %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) 

head(sc_harmony@meta.data)
table(sc_harmony@meta.data$seurat_clusters)

all_markers <- FindAllMarkers(sc_harmony,
                              logfc.threshold = 0.5,
                              only.pos = T, 
                              min.pct = 0.5)
head(all_markers)
all_markers %>% 
  group_by(cluster) %>% 
  filter(p_val_adj <= 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(20)) %>%
  View
top_genes <- all_markers %>% 
  group_by(cluster) %>% 
  filter(p_val_adj <= 0.05) %>% 
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(2)) %>%
  View

heat <- t(Seurat::AverageExpression(sc_harmony, 
                                    assays = "MAGIC_RNA",
                                    group.by = "seurat_clusters",
                                    verbose = T,
                                    features = top_genes$gene)$MAGIC_RNA)
breaks <- seq(-1,1,length.out=50)
clrs <- colorRampPalette(colors = c("#87CEEB", "white", "#F08080"))(length(breaks))
rws <- rownames(heat)
cls <- colnames(heat)
heat <- apply(heat,2,scale)
colnames(heat) <- cls
rownames(heat) <- rws
Heatmap(t(heat), 
        cluster_rows = F, 
        cluster_columns = F, 
        row_dend_side = "right",row_names_side = "left",
        name = paste("Scaled", "Mean", "Expression", sep="\n"),
        col = clrs)

# run magic ####
sc_harmony <- magic(sc_harmony)
DefaultAssay(sc_harmony) <- "MAGIC_RNA"


# start checking markers 
head(sc_harmony@meta.data)
rownames(sc_harmony)[grep("ccr7", rownames(sc_harmony), ignore.case = T)]

eatureScatter(sc_harmony, 
              feature1 = "Ptprc",
              feature2 = "Ccr7", 
              plot.cor = F,
              cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                       '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                       '12' = c25[13])) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2))
DimPlot(sc_harmony, reduction = "umap", 
        split.by = "orig.ident",
        cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                 '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                 '12' = c25[13]))

VlnPlot(sc_harmony, features = c("Cd8a", "Cd4", "Cd14",
                                 "Ly6c1", "Ccr7","Sell", 
                                 "Gzma", "Ifng", "Lamp1",
                                 "Foxp3", "Ctla4", "Pdcd1"), 
        split.plot = T,
        group.by = "seurat_clusters",
        # log = T,
        split.by = "orig.ident",
        pt.size = 0, 
        ncol = 3, 
        cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                 '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                 '12' = c25[13])
)


# >> use Project TILs to classify cells ####
# BiocManager::install(c("AUCell"))
# remotes::install_github("carmonalab/TILPRED")
# remotes::install_github("carmonalab/UCell")
# remotes::install_github("carmonalab/scGate")
# install.packages("ggparty")
# remotes::install_github("carmonalab/ProjecTILs")
ref <- load.reference.map(ref = "~/Google Drive/Protocols_bfx/ref_LCMV_CD4_mouse_release_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref,label = T, cols = refCols)
markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
DefaultAssay(sc_harmony) <- "RNA"
query.projected <- make.projection(sc_harmony, ref=ref, filter.cells = T)

plot.projection(ref, query.projected, linesize = 0.5, pointsize = 0.5) 
query.projected <- cellstate.predict(ref=ref, query=query.projected)
table(query.projected$functional.cluster)
boxplot(query.projected$functional.cluster.conf ~ query.projected$functional.cluster)
View(cbind(head(query.projected),head(sc_harmony)))


# merge the three TILPRED references ####
sc_harmony$TILPRED_TIL <- query.projected$functional.cluster
sc_harmony$TILPRED_TIL_score <- query.projected$functional.cluster.conf
DimPlot(sc_harmony, group.by = "TILPRED_TIL")
FeaturePlot(sc_harmony, features = "TILPRED_TIL_score", min.cutoff = 0.75)

sc_harmony$TILPRED_CD8_viral <- query.projected$functional.cluster
sc_harmony$TILPRED_CD8_viral_score <- query.projected$functional.cluster.conf
DimPlot(sc_harmony, group.by = "TILPRED_CD8_viral")
FeaturePlot(sc_harmony, features = "TILPRED_CD8_viral_score", min.cutoff = 0.75)

sc_harmony$TILPRED_CD4_viral <- query.projected$functional.cluster
sc_harmony$TILPRED_CD4_viral_score <- query.projected$functional.cluster.conf
DimPlot(sc_harmony, group.by = "TILPRED_CD4_viral")
FeaturePlot(sc_harmony, features = "TILPRED_CD4_viral_score", min.cutoff = 0.75)

# cleanup columns I don't need from annotation to prep for cell type final annotation ####
head(sc_harmony)
sc_harmony$X <- NULL
sc_harmony$functional.cluster <- NULL
sc_harmony$TILPRED_TILs_ref <- NULL

# cleanup and harmonize labels for TILPRED ####
table(sc_harmony$TILPRED_CD4_viral)
table(sc_harmony$TILPRED_TIL)
table(sc_harmony$TILPRED_CD8_viral)

sc_harmony$TILPRED_CD8_viral.2 <- ifelse(sc_harmony$TILPRED_CD8_viral == "Eff Cycling", "CD8_Cycling_Effector",
                                         ifelse(sc_harmony$TILPRED_CD8_viral == "Eff Early", "CD8_Early_Activated",
                                                ifelse(sc_harmony$TILPRED_CD8_viral == "Eff Interm", "CD8_Intermediate_Memory",
                                                       ifelse(sc_harmony$TILPRED_CD8_viral == "Memory Prec", "CD8_Precursor_Memory",
                                                              ifelse(sc_harmony$TILPRED_CD8_viral == "SLEC", "CD8_ShortLived_Effector",
                                                                     ifelse(sc_harmony$TILPRED_CD8_viral == "Tex", "CD8_Exhausted",
                                                                            ifelse(sc_harmony$TILPRED_CD8_viral == "Tpex", "CD8_Precursor_Exhausted","Other")))))))
table(sc_harmony$TILPRED_CD8_viral.2)
table(sc_harmony$TILPRED_CD8_viral.2, sc_harmony$TILPRED_CD8_viral)

sc_harmony$TILPRED_TIL.2 <- ifelse(sc_harmony$TILPRED_TIL == "CD4_NaiveLike", "CD4_Naive_CM",
                                   ifelse(sc_harmony$TILPRED_TIL == "CD8_EarlyActiv", "CD8_Early_Activated",
                                          ifelse(sc_harmony$TILPRED_TIL == "CD8_EffectorMemory", "CD8_Effector_Memory",
                                                 ifelse(sc_harmony$TILPRED_TIL == "CD8_NaiveLike", "CD8_Naive_CM",
                                                        ifelse(sc_harmony$TILPRED_TIL == "CD8_Tex", "CD8_Exhausted",
                                                               ifelse(sc_harmony$TILPRED_TIL == "CD8_Tpex", "CD8_Precursor_Exhausted",
                                                                      ifelse(sc_harmony$TILPRED_TIL == "Tfh", "CD4_Follicular_Helper",
                                                                             ifelse(sc_harmony$TILPRED_TIL == "Th1", "CD4_Th1_Helper",
                                                                                    ifelse(sc_harmony$TILPRED_TIL == "Treg", "CD4_Treg","Other")))))))))

table(sc_harmony$TILPRED_TIL.2)
table(sc_harmony$TILPRED_TIL.2, sc_harmony$TILPRED_TIL)

sc_harmony$TILPRED_CD4_viral.2 <- ifelse(sc_harmony$TILPRED_CD4_viral == "Eomes_HI", "CD4_Eomese_High",
                                         ifelse(sc_harmony$TILPRED_CD4_viral == "INFI_stimulated", "CD4_Follicular_Helper_Effector",
                                                ifelse(sc_harmony$TILPRED_CD4_viral == "Tcm", "CD4_CM",
                                                       ifelse(sc_harmony$TILPRED_CD4_viral == "Tcmp", "CD4_Precursor_Memory",
                                                              ifelse(sc_harmony$TILPRED_CD4_viral == "Tfh_Effector", "CD4_Follicular_Helper_Effector",
                                                                     ifelse(sc_harmony$TILPRED_CD4_viral == "Tfh_Memory", "CD4_Follicular_Helper_Memory",
                                                                            ifelse(sc_harmony$TILPRED_CD4_viral == "Th1_Effector", "CD4_Th1_Effector",
                                                                                   ifelse(sc_harmony$TILPRED_CD4_viral == "Th1_Memory", "CD4_Th1_Memory",
                                                                                          ifelse(sc_harmony$TILPRED_CD4_viral == "Treg", "CD4_Treg","Other")))))))))
table(sc_harmony$TILPRED_CD4_viral.2)
table(sc_harmony$TILPRED_CD4_viral.2, sc_harmony$TILPRED_CD4_viral)

# calculate highest score per each cell type and keep that for TILPRED ####
colnames(sc_harmony@meta.data)
sc_harmony$TILPRED_harmonized <- unlist(apply(sc_harmony@meta.data[,18:26], 1, function(xx){
  temp_scores <- xx[c(2,4,6)]
  temp_max <- which(temp_scores == max(temp_scores, na.rm = T))[1]
  temp_cell_ind <- ifelse(temp_max == 1, 9, 
                          ifelse(temp_max == 2, 8, 7))
  temp_cell <- xx[temp_cell_ind]
  temp_assigned <- ifelse(temp_scores[temp_max] >= quantile(sc_harmony@meta.data[,(temp_max + 6)], 0.5)[[1]], temp_cell, paste0("other_", temp_cell))
  return(temp_cell)
}))
table(sc_harmony$TILPRED_harmonized)
DimPlot(sc_harmony, 
        group.by = "TILPRED_harmonized", 
        cols = c25,
        # label = T, 
        # label.box = T, 
        # repel = T, 
        pt.size = 1, 
        split.by = "TILPRED_harmonized", ncol = 10) 

table(sc_harmony$TILPRED_harmonized, sc_harmony$TILPRED_CD8_viral.2)

DefaultAssay(sc_harmony) <- "MAGIC_RNA"
DimPlot(sc_harmony, group.by = "CellType_Seurat", cols = c25)

table(sc_harmony$TILPRED_harmonized.2, sc_harmony$seurat_clusters)
rownames(sc_harmony)[grep("cd103", rownames(sc_harmony), ignore.case = T)]
FeaturePlot(sc_harmony, 
            cols = c("lightgray","red"),
            features = c("Cd4", "Cd8a", "Klrb1a", "Cd14"))
FeatureScatter(sc_harmony,
               feature1 = "Cd8a",
               feature2 = "Cd4",
               group.by = "seurat_clusters")
VlnPlot(sc_harmony,
        features = c("Cd8a", "Cd4"),
        group.by = "seurat_clusters",
        stack = T,
        pt.size = 0)
VlnPlot(sc_harmony,
        features = c("Cd8a", "Cd4","Gzmk", "Cd69", "Ccr7", "Sell", "Itgae","Cx3cr1", "Lag3", "Pdcd1", "Cxcr6", "Ifng"),
        group.by = "TILPRED_harmonized.2",
        stack = T,
        pt.size = 0)
sc_harmony@assays$MAGIC_RNA@data["Foxp3",1:4]
ggplot(gg_t[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory"),], 
       aes(x=Ifng, y=Gzmb)) +
  geom_point(size=0.3) +
  guides(color = guide_legend(override.aes = list(size=5)))  +
  facet_wrap(~ factor(sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory")]))
# adjust final labels ####
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "9")] <- "CD8_Effector_Exhauted"
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "0")] <- "CD8_Effector_Memory"
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "1")] <- "CD8_Memory_Exhausted"
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "5")] <- "CD8_Effector"
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "7" &
                                        sc_harmony@assays$MAGIC_RNA@data["Cd8a",] > 1 )] <- "CD8_CM"
sc_harmony$TILPRED_harmonized.2[which(sc_harmony$seurat_clusters == "7" &
                                        sc_harmony@assays$MAGIC_RNA@data["Cd4",] > 0.5 )] <- "CD4_CM"
# +>+> save temp tilpred.2 before step below 12152021 ####
temp_tilpred_harmonized.2 <- sc_harmony$TILPRED_harmonized.2
sc_harmony$TILPRED_harmonized.2 <- ifelse(sc_harmony$TILPRED_harmonized == "CD4_Th1_Effector", "CD8_Effector_Memory",
                                          ifelse(sc_harmony$TILPRED_harmonized == "CD4_Eomese_High", "CD8_Precursor_Exhausted", 
                                                 ifelse(sc_harmony$seurat_clusters == "11", "Monocytes", 
                                                        ifelse(sc_harmony$seurat_clusters == "0" & sc_harmony$TILPRED_harmonized == "CD8_ShortLived_Effector", "CD8_Effector", 
                                                               ifelse(sc_harmony$TILPRED_harmonized == "CD4_Treg" & 
                                                                        sc_harmony@assays$MAGIC_RNA@data["Foxp3",] > 1 & 
                                                                        sc_harmony@assays$MAGIC_RNA@data["Gzmb",] > 1.5, "CD4_Treg_Effector", 
                                                                      ifelse(sc_harmony$TILPRED_harmonized == "CD8_Precursor_Exhausted" & 
                                                                               sc_harmony$seurat_clusters == "3", "CD4_Th1_Helper",
                                                                             ifelse(sc_harmony$TILPRED_harmonized == "CD8_ShortLived_Effector" & 
                                                                                      sc_harmony$seurat_clusters == "3", "CD4_Th1_Helper",
                                                                                    ifelse(sc_harmony$TILPRED_harmonized == "CD8_Precursor_Exhausted" & 
                                                                                             sc_harmony$seurat_clusters == "2", "CD4_Treg_Effector", 
                                                                                           ifelse(sc_harmony$TILPRED_harmonized == "CD4_Th1_Memory" & 
                                                                                                    sc_harmony$seurat_clusters == "0", "CD8_Effector_Memory", 
                                                                                                  ifelse(sc_harmony$TILPRED_harmonized == "CD8_Precursor_Exhausted" & 
                                                                                                           sc_harmony$seurat_clusters == "7", "CD8_Effector_Memory",
                                                                                                         ifelse(sc_harmony$TILPRED_harmonized == "CD8_Intermediate_Memory" & 
                                                                                                                  sc_harmony$seurat_clusters == "9", "CD8_Effector_Memory",
                                                                                                                ifelse(sc_harmony@assays$MAGIC_RNA@data["Ccr7",] > 0.75 & 
                                                                                                                         sc_harmony@assays$MAGIC_RNA@data["Sell",] > 0.75, "CD4_CM",
                                                                                                                       ifelse(sc_harmony@assays$MAGIC_RNA@data["Cd8a",] > 1 & 
                                                                                                                                sc_harmony@assays$MAGIC_RNA@data["Gzmb",] > 0.75 &
                                                                                                                                sc_harmony$seurat_clusters == "5", "CD8_Effector",
                                                                                                                              ifelse(sc_harmony$TILPRED_harmonized == "CD4_Th1_Helper" & 
                                                                                                                                       sc_harmony$seurat_clusters %in% c("0", "5"), "CD8_Effector_Memory",
                                                                                                                                     ifelse(sc_harmony$TILPRED_harmonized == "CD4_Th1_Helper" & 
                                                                                                                                              sc_harmony$seurat_clusters == "2", "CD8_CM", 
                                                                                                                                            ifelse(sc_harmony$TILPRED_harmonized == "CD8_Naive_CM" & 
                                                                                                                                                     sc_harmony$seurat_clusters == "0", "CD8_Effector_Memory", 
                                                                                                                                                   ifelse(sc_harmony$TILPRED_harmonized == "CD4_Th1_Helper" & 
                                                                                                                                                            sc_harmony$seurat_clusters == "1", "CD8_CM", 
                                                                                                                                                          ifelse(sc_harmony$TILPRED_harmonized == "CD8_Effector_Memory" & 
                                                                                                                                                                   sc_harmony$seurat_clusters == "2", "CD8_Exhausted", 
                                                                                                                                                                 ifelse(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory" & 
                                                                                                                                                                          sc_harmony$seurat_clusters == "1", "CD8_Exhausted", sc_harmony$TILPRED_harmonized)))))))))))))))))))
# +>+> save temp tilpred.2 before step below 12152021 ####
temp_tilpred_harmonized.2_after_step_above <- sc_harmony$TILPRED_harmonized.2
table(sc_harmony$TILPRED_harmonized.2)

colnames(gg_t)[grep("Cd25",colnames(gg_t), ignore.case = T)]
ggplot(gg_t[which(sc_harmony$TILPRED_harmonized.2 == "CD4_CM"),], 
       aes(x=Cd4, y=Cd44, col=Il7r)) +
  geom_point(size=1.3) +
  scale_color_gradient(low = "cyan",  high = "red") +
  theme_cowplot() +
  facet_wrap(~ factor(sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD4_CM")]))

DimPlot(sc_harmony, 
        group.by = "TILPRED_harmonized.2", 
        cols = c25,
        pt.size = 1, 
        split.by = "TILPRED_harmonized.2", ncol = 10) 
DimPlot(sc_harmony, 
        group.by = "TILPRED_harmonized.2", 
        cols = c25,
        pt.size = 1)



table(sc_harmony$TILPRED_CD4_viral.2)
table(sc_harmony$TILPRED_CD4_viral.2, sc_harmony$TILPRED_CD4_viral)
DimPlot(sc_harmony, group.by = "seurat_clusters", label = T, label.box = T, repel = T)

FeatureScatter(sc_harmony,
               feature1 = "Sell",
               feature2 = "Cd8a", 
               cols = c25,
               group.by = "TILPRED_harmonized.2")
FeatureScatter(sc_harmony,
               feature1 = "Cd4",
               feature2 = "TILPRED_TIL_score", 
               group.by = "TILPRED_harmonized")
Heatmap(table(sc_harmony$CellType_Seurat, sc_harmony$TILPRED_harmonized.2), border = T,
        col = colorRamp2(breaks = c(0,5,20,200), colors = c("white", "#F8F8FF", "#AFEEEE", "#4169E1")))

Heatmap(table(sc_harmony$TILPRED_harmonized.2, sc_harmony$seurat_clusters), 
        border = T, 
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,10,20,100), colors = c("white", "#F8F8FF", "#AFEEEE", "#4169E1")))

tt <- table(sc_harmony$TILPRED_harmonized.2, sc_harmony$orig.ident)
tt <- melt(tt)
ggplot(tt, aes(x=Var1, y=log10(value+1), fill=Var2)) +
  geom_bar(stat="identity", position="dodge", col="darkgray") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size=15)) +
  # theme(axis.text.x = element_blank()) +
  scale_fill_manual(name="Group",
                    values = c("CTR" = "lightgray", "Treated" = "tomato")) 

# make heatmap cell types and markers ####
genes2check <- c("Cd8a", "Cd8b1", "Gzmk", "Gzmb", "Cd69", "Ccr7", "Sell",  "Ifng",
                 "Cd4", "Itgae","Cx3cr1", "Ctla4", "Lag3", "Pdcd1", "Cxcr6", 
                 "Cd14")
heat <- gg_t[,which(colnames(gg_t) %in%  genes2check)]
heat <- aggregate(heat, by=list(sc_harmony$TILPRED_harmonized.2), FUN=mean)
heat <- t(heat)
colnames(heat) <- heat[1,]
heat <- heat[-1,]
cls <- colnames(heat)
heat <- t(apply(heat,1,as.numeric))
heat <- t(apply(heat, 1, scale))
colnames(heat) <- cls
heat <- heat[genes2check,c(grep("CD8", colnames(heat)),
                           grep("CD4", colnames(heat)),
                           grep("Mono", colnames(heat)))]

heat[1:5,1:5]
Heatmap(heat, 
        cluster_rows = F,
        cluster_columns = F,
        col = colorRamp2(breaks = c(-2,0,2), 
                         colors = c("white", "#F5F5F5", "black")),
        border = T)


ggplot(gg_t, 
       aes(x=Ifng, y=Gzmb)) +
  theme_bw() +
  geom_point(size=0.3) +
  guides(color = guide_legend(override.aes = list(size=5)))  +
  facet_wrap(sc_harmony$TILPRED_harmonized.2 ~ factor(sc_harmony$orig.ident))

# run SCTransform ####
sc_harmony <- SCTransform(sc_harmony)

# run DEG ####
time_points <- as.character(unique(sc_harmony$orig.ident))
cell_pop <- unique(sc_harmony$TILPRED_harmonized.2)
for(ii in 1:length(cell_pop)){
  if(ii == 1)
    res_deg_all <- list()
  
  temp.cells <- subset(sc_harmony, cells=which(sc_harmony$TILPRED_harmonized.2 == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "RNA"
  Idents(temp.cells) <- as.character(temp.cells$orig.ident)
  
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs 
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.2, 
                                   min.pct = 0.2),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  res_deg_all[[ii]] <- data.frame(rbind(temp.res))
  names(res_deg_all)[ii] <- cell_pop[ii]
  
  message(paste(" >> Done for",cell_pop[ii], "<<"))
  message(paste("   - >> ",sum(temp.res$p_val_adj < 0.05), " DEGs "))
}

# plot DEGs ####
res_all <- data.table::rbindlist(res_deg_all, idcol = T)
head(res_all)
colnames(res_all)[1] <- "CellType"
res_all$Significant <- ifelse(res_all$p_val_adj < 0.05, "YES", "NO")
table(res_all$CellType[which(res_all$Significant == "YES")])
res_all$Direction <- ifelse(res_all$Significant == "YES" & res_all$avg_log2FC > 0, "UP",
                            ifelse(res_all$Significant == "YES" & res_all$avg_log2FC < 0, "DOWN", "NS"))
res_all$size <- ifelse(res_all$Direction == "UP", 0.5,
                       ifelse(res_all$Direction == "DOWN", 0.5, 0.2))
View(res_all)
as.data.frame.matrix(table(res_all$CellType, res_all$Direction)) %>%
  mutate(Ratio = apply(., 1, function(xx) round(sum(xx[1] + xx[3])/sum(xx), 3)))


tt_deg <- table(res_all$Direction, res_all$CellType)
Heatmap(tt_deg[-2,],
        col=colorRamp2(breaks = c(0,5,50,200), 
                       colors = c("#FFFFF0", "#FAF0E6","#DDA0DD", "#0000CD")),
        name="Number DEG")

ggplot(res_all,
       aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
  geom_vline(xintercept = 0, size=0.75, col="gray", linetype="dotted") +
  geom_hline(yintercept = -log10(0.05), size=0.75, col="gray", linetype="dotted") +
  geom_point(alpha = 0.7) +
  xlim(-3,3) +
  theme_minimal() +
  facet_wrap( ~ CellType, scales = "free") +
  scale_color_manual(values = c("DOWN" = "cyan",
                                "UP" = "tomato",
                                "NS" = "lightgray"))

res_all_sig <- res_all %>%
  group_by(CellType) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

res_all_sig[grep("Pdcd1", res_all_sig$Gene, ignore.case = T),]
set_ortho <- read.table("~/Google Drive/Protocols_bfx/Ortholog_Human_Mouse_genes.txt", header = T, sep="\t")
head(set_ortho)
set_ortho[grep("^Il15", set_ortho$human, ignore.case = T),]
res_all_sig[which( res_all_sig$Gene %in% c("Pdcd1","Lag3", "Tck", "Zap70", "Ifnar1", "Ifng","Gzmb", "Gzma")),]
label <- c("Pdcd1","Lag3", "Tck", "Zap70", "Ifnar1", "Ifng","Gzmb", "Gzma", "Ccr7", "Cd69", "Rela", "Nfkb1", "Nfkb2", "Ccr5",
           "Lat", "Cd28", "Lcp2", "Ccl22", "Il7", "Il15", "Cd27", "Tnf", "Ctla4", "Ccl5", "Cxcr3", "Cx3cr1")
label[which(label %in% res_all_sig$Gene)]
cells <- as.character(unlist(unique(res_all_sig[which( res_all_sig$Gene %in% label),1])))

# plot labelled DEGs ####
for(ii in 1:length(cells)){
  subs <- res_all_sig[which(res_all_sig$CellType == cells[ii]),]
  subs$labels <- ifelse(subs$Gene %in% label, subs$Gene, "")
  gg <- ggplot(subs,
               aes(x=avg_log2FC, y=-log10(p_val_adj), col=Direction)) +
    geom_vline(xintercept = 0, size=0.75, col="gray", linetype="dotted") +
    geom_hline(yintercept = -log10(0.05), size=0.75, col="gray", linetype="dotted") +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    scale_color_manual(values = c("DOWN" = "cyan",
                                  "UP" = "tomato",
                                  "NS" = "lightgray")) +
    geom_label_repel(aes(label = labels), color="black",
                     box.padding   = 0.35, 
                     point.padding = 0.5, 
                     segment.color = 'grey50', 
                     max.overlaps = 100) +
    ggtitle(cells[ii])
  print(gg)
}

# run GSEA using ortholog data ####
orth_db <- read.table("~/Google Drive/Protocols_bfx/human_mouse_hcop_fifteen_column_Jan2022.txt", sep = "\t", header = T, quote = "")
orth_db <- orth_db[,c("human_entrez_gene", "human_name", "mouse_entrez_gene", "mouse_symbol")]
colnames(orth_db)
orth_db[grep("1302", orth_db$mouse_entrez_gene),]
head(orth_db)
sig_all <- readRDS("~/Google Drive/Protocols_bfx/GSEA_signatures/Mm.h.all.v7.1.entrez.rds")
sum(lapply(sig_all, function(xx) xx[which(xx == "13024")]) > 0, na.rm = T)
lapply(sig_all, function(xx) xx <- orth_db$mouse_symbol[which(orth_db$mouse_entrez_gene %in% xx)])
load("~/Google Drive/Protocols_bfx/GSEA_signatures/Mm.c2.cp.reactome.v6.1.entrez.RData")
Mm.c2.cp.reactome.v6.1.entrez
load("~/Google Drive/Protocols_bfx/GSEA_signatures/Mm.c3.tft.v6.1.entrez.RData")
Mm.c3.tft.v6.1.entrez
library(dplyr)
run_gsea_sc <- function(Comparison, ref_signature){

  temp_comparison <- res_all %>% filter(CellType == Comparison) %>% as.data.frame()
  table(temp_comparison$Direction)
  
  message("Setting up signatures")
  set_ortho <- orth_db
  ref_signature <- lapply(ref_signature, function(xx) xx <- orth_db$mouse_symbol[match(xx,orth_db$mouse_entrez_gene)])
  message("  >> Done")

  ranked_ii <- temp_comparison$avg_log2FC
  names(ranked_ii) <- temp_comparison$Gene
  ranked_ii <- sort(ranked_ii, decreasing = T)
  
  message(paste("Running fgsea for",temp_comparison[1,1]))
  temp_fgsea <- tryCatch(fgsea::fgsea(pathways = ref_signature,
                                      stats = ranked_ii),
                         error = function(xx){
                           message(xx)
                           message("\nadding NAs\n")
                           dummy_list <- NA
                           return(dummy_list)
                         })

  temp_fgsea$CellType <- Comparison
  temp_fgsea$leadingEdge <- as.character(temp_fgsea$leadingEdge)
  
  message(paste("Done fgsea for >",Comparison))
  message("")
  return(temp_fgsea)
}

for(ii in 1:length(unique(res_all$CellType))){
  if(ii == 1){
    res_all_fgsea2 <- list()
    comp_ref <- unique(res_all$CellType)
  }
  temp <- run_gsea_sc(comp_ref[ii], Mm.c3.tft.v6.1.entrez)
  
  res_all_fgsea2[[ii]] <- temp
  names(res_all_fgsea2)[ii] <- comp_ref[ii]
}

head(res_all_fgsea2$CD8_Memory_Exhausted)

# extrat GSEA results ####
gsea_res <- rbindlist(res_all_fgsea2, idcol = T)
colnames(gsea_res)[1] <- "CellType_ref"
gsea_res$Direction <- ifelse(gsea_res$padj < 0.05 & gsea_res$NES > 0, "UP", 
                             ifelse(gsea_res$padj < 0.05 & gsea_res$NES < 0, "DOWN", "NS"))

gsea_res %>%
  filter(grepl("CD8_Effector", gsea_res$CellType, ignore.case = T)) %>%
  filter(padj < 0.05) %>%
  mutate(pval = NULL) %>%
  mutate(log2err = NULL) %>%
  mutate(ES = NULL) %>%
  mutate(size = NULL) %>%
  mutate(leadingEdge = NULL) %>%
  View

gsea_res_reactome <- gsea_res
gsea_res_hallmark <- gsea_res
gsea_res_TF <- gsea_res
gsea_res_all <- data.frame(rbind(gsea_res_reactome,
                                 gsea_res_hallmark,
                                 gsea_res_TF))
gsea_res_all %>%
  filter(grepl("CD8_Memory_Exhausted", gsea_res_all$CellType, ignore.case = T)) %>%
  filter(padj < 0.05) %>%
  mutate(pval = NULL) %>%
  mutate(log2err = NULL) %>%
  mutate(ES = NULL) %>%
  mutate(size = NULL) %>%
  mutate(leadingEdge = NULL) %>%
  View

write.table(gsea_res_all, "all_GSEA_results_01172022.txt", sep="\t", quote = F, row.names = F)

table(gsea_res$Direction[grep("CD8_Memory_Exhausted", gsea_res$CellType)])
table(old_gsea_res$Direction[grep("CD8_Effector_gmt", old_gsea_res$CellType)])

View(gsea_res[grep("CD8_Effector_gmt", gsea_res$CellType),])
View(old_gsea_res[grep("CD8_Effector_gmt", old_gsea_res$CellType),])

vv <- gplots::venn(list(NEW = gsea_res$pathway[grep("CD8_Effector_gmt", gsea_res$CellType)],
                        OLD= old_gsea_res$pathway[grep("CD8_Effector_gmt", old_gsea_res$CellType)]))
View(attributes(vv)[3])
table(gsea_res$CellType)

gsea_sign <- gsea_res %>%
  filter(padj < 0.05) %>%
  group_by(CellType) %>%
  arrange(desc(abs(NES)))
head(gsea_sign)

gsea_sign$leadingEdge <- as.character(gsea_sign$leadingEdge)
View(apply(gsea_sign$leadingEdge, 2, as.character))
# write.table(gsea_sign[,1:7], "GSEA_significant_pathways_last.txt", sep="\t", quote=F, row.names = F)
# gsea_sign <- read.table("GSEA_significant_pathways.txt", header = T, sep="\t")
table(gsea_sign$CellType)
View(gsea_sign[grep("CD8_Effector", gsea_sign$CellType, ignore.case = T),])


dd_gsea <- data.frame(rbind(gsea_sign[grep("cytoto", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("chemok", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("cytok", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("t_cell_re", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("alpha_beta_T", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("T_cell_proli", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("antige", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("leuko", gsea_sign$pathway, ignore.case = T),],
                            gsea_sign[grep("lymph", gsea_sign$pathway, ignore.case = T),]))[,c(1,2,3,4,7)]
unique(dd_gsea$CellType)
dd_gsea <- dd_gsea %>% 
  filter(CellType != "Monocytes_gmt_gobp")
# filter(CellType %in% unique(dd_gsea$CellType)[c(2:6,8:10,12,14,17,18)])
dd_gsea <- dd_gsea[order(dd_gsea$NES, decreasing = T),]
dd_gsea$pathway <- ifelse(duplicated(dd_gsea$pathway), paste0(dd_gsea$pathway,".1"), dd_gsea$pathway)
# ggplot(dd_gsea, aes(x=factor(pathway, levels = dd_gsea$pathway), y=rev(NES), fill=CellType)) +
# [grep("^GO_", dd_gsea$pathway, ignore.case = T)]
ggplot(dd_gsea, 
       aes(x=factor(pathway, levels = dd_gsea$pathway), 
           y=rev(NES), 
           fill=CellType)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = c25) +
  xlab("NES") + ylab("")

library(kableExtra)
dd_gsea %>%
  kbl() %>%
  kable_styling()

gsea_res %>%
  filter(NES > 0) %>%
  filter(grepl("CD8_Effector_Memory", CellType)) %>%
  filter(padj < 0.25)


features <- c("IFNAR1", "IFNAR2", "IRF3", "NFKB1" ,"NFKB2", "RELA", "RELB", "CXCL12", 
              "CCL22", "CCL5","IFNG", "GZMB", "TNF", "CXCL10", "STAT1", "STAT2", "IRF9",
              "Zap70", "Lck")
features[-which(features %in% colnames(gg_t))]
cell_type <- unique(sc_rna$SingleR_Covid)
time_point <- unique(sc_rna$orig.ident)
table(sc_rna$orig.ident, sc_rna$SingleR_Covid)
for(ii in 1:length(cell_type)){
  message(paste("Working on", cell_type[ii]))
  for(kk in 1:length(time_point)){
    temp_cell_type_ind <- which(rownames(gg_t) %in% rownames(sc_rna@meta.data)[which(sc_rna$SingleR_Covid == cell_type[ii])])
    temp_time_point_ind <- which(rownames(gg_t) %in% rownames(sc_rna@meta.data)[which(sc_rna$orig.ident == time_point[kk])])
    temp_cells_ind <- intersect(temp_cell_type_ind,temp_time_point_ind)
    
    if(length(temp_cells_ind) == 0){
      message(paste("skipping", cell_type[ii]))
      next
    }
    temp_features_ind <- which(colnames(gg_t) %in% unique(c(features, other_genes, "TBK1", "IFITM3", "IFIT1", "MX1", "EIF2AK2")))
    gg_temp <- gg_t[temp_cells_ind, temp_features_ind]
    cor_mat <- tryCatch(cor(gg_temp, method = "pearson", use = "complete.obs"),
                        error = function(xx){
                          message(xx)
                          dummy_df <- cor(gg_temp, method = "pearson", use = "complete.obs")
                          return(NA)
                        })
    #cor_mat[1:5,1:5]
    mm_cor_mat <- reshape2::melt(cor_mat)
    mm_cor_mat$CellType <- cell_type[ii]
    mm_cor_mat$TimePoint <- time_point[kk]
    head(mm_cor_mat)
    if(kk == 1){
      temp_res_cor <- mm_cor_mat
    } else {
      temp_res_cor <- rbind(temp_res_cor, mm_cor_mat)
    }
  }
  if(ii == 1){
    res_cor <- list(temp_res_cor)
    names(res_cor)[ii] <- cell_type[ii]
  } else {
    res_cor[[ii]] <- temp_res_cor
    names(res_cor)[ii] <- cell_type[ii]
  }
  message(paste("  > Done for", cell_type[ii]))
}

names(res_cor)
lapply(res_cor, dim)
head(res_cor)

# make network for correlation plots ####
res_cor_network <- rbindlist(res_cor)
colnames(res_cor_network)[1:3] <- c("Node1", "Node2", "CorrCoeff")
res_cor_network$plot <- ifelse(res_cor_network$CorrCoeff > 0.7 & res_cor_network$CorrCoeff < 1, "YES", "NO")
res_cor_network$clr <- ifelse(res_cor_network$plot == "YES", "tomato", "#FFF8DC")
table(res_cor_network$TimePoint, res_cor_network$CellType)
head(res_cor_network)
table(res_cor_network$plot)


time_point <- c("C1D1_pre", "C1D1_post", "C1D4_pre", "C1D4_post", "C1D7")
for(ii in 1:length(cell_type)){
  if(ii == 1){
    pdf("DyNA_cytotoxic_per_cell_v2.pdf", width = 14, height = 4)
    par(mfrow=c(1,5))
    # par(mar=c(2,2,1,2))
  }
  message(paste("Working on", cell_type[ii]))
  for(kk in 1:length(time_point)){
    net <- res_cor_network %>%
      subset(CellType == cell_type[ii]) %>%
      subset(TimePoint == time_point[kk])
    
    g <- graph.data.frame(net, directed = FALSE)
    V(g)$color <- ifelse(V(g)$name %in% net$Node1[net$plot == "YES"], "tomato", "#FFF8DC")
    E(g)$color <- ifelse(is.na(edge_attr(g)$plot), scales::alpha("white", 0.1), 
                         ifelse(edge_attr(g)$plot == "YES","black", scales::alpha("white", 0.1)))
  
    plot(g, 
         edge.width = ifelse(is.na(edge_attr(g)$plot), 0.001,
                             ifelse(edge_attr(g)$plot == "YES", 2, 0.001)),
         edge.curved = 0,
         layout = layout_in_circle, 
         main = paste(cell_type[ii], time_point[kk], sep="\n"),
         vertex.label.color="black",
         vertex.label.font=2)
    # print(nplot)
    message(paste("Done for", cell_type[ii], time_point[kk]))
    if(ii == 1 & kk == 1){
      res_cor_plot <- data.frame(CellType = cell_type[ii],
                                 TimePoint = time_point[kk],
                                 CorVertices = sum(net$plot == "YES", na.rm = T))
    } else {
      res_cor_plot <- data.frame(rbind(res_cor_plot, c(cell_type[ii],
                                                       time_point[kk],
                                                       sum(net$plot == "YES"))))
    }
  }
  if(ii == length(cell_type)){
    dev.off()
    message("  >> Done")
  }
}
head(res_cor_plot)
res_cor_plot$TP_number <- ifelse(res_cor_plot$TimePoint == "C1D1_pre", 1, 
                                 ifelse(res_cor_plot$TimePoint == "C1D1_post", 2, 
                                        ifelse(res_cor_plot$TimePoint == "C1D4_pre", 3, 
                                               ifelse(res_cor_plot$TimePoint == "C1D4_post", 4, 
                                                      ifelse(res_cor_plot$TimePoint == "C1D7", 5, -9)))))

table(res_cor_plot$CellType)
ggplot(res_cor_plot, 

       aes(x = as.numeric(TP_number), 
           y=as.numeric(CorVertices))) +
  geom_point(size=1.5) +
  geom_smooth(method = "loess", lwd=1.5) +
  theme_classic() +
  scale_x_continuous(n.breaks = 5,
                     labels = c("C1D1_pre", "C1D1post", "C1D4pre", "C1D4post", "C1D7")) +
  facet_wrap(~ CellType, scale="free", ncol=7) +
  scale_color_manual(values = c25) +
  theme(axis.text.x = element_text(angle=90))



# do ECDF plots ####
gg_t <- data.frame(t(sc_harmony@assays$MAGIC_RNA@data))
gg_t[1:5,1:5]

ggplot(gg_t, aes(x=Gzma, y=factor(as.character(sc_harmony@meta.data$seurat_clusters), levels = 0:20), 
                 fill=sc_harmony@meta.data$orig.ident)) +
  geom_density_ridges(alpha=0.3) +
  theme_bw() +
  ylab("")


rownames(sc_harmony)[grep("^rela", rownames(sc_harmony), ignore.case = T)]
cell_line <- c("TCR-signalling CD8+ T cells")
ggplot(gg_t, 
       aes(x=(Cd8a * Lck * Zap70), 
           col=factor(sc_harmony$orig.ident))) +
  stat_ecdf(geom = "step", size=1.5) +
  theme_bw() +
  labs(col = "") +
  ylab("cumulative frequency") +
  ggtitle(cell_line) +
  # ggtitle("CD4/Treg") +
  scale_color_manual(values = c("CTR" = "blue",
                                "Treated" = "tomato")) +
  facet_wrap(~ sc_harmony$seurat_clusters)

v_c <- ecdf(rowMeans(gg_t[ind_cells,
                          which(colnames(gg_t) %in% c("GZMK","IFNG"))]))
v_t <- ecdf(rowMeans(gg_t[ind_cells,
                          which(colnames(gg_t) %in% c("GZMK","IFNG"))]))
ks.test(knots(v_t),knots(v_c))

FeatureScatter(sc_harmony, 
               feature1 = "Cd8b1",
               feature2 = "Zap70", 
               pt.size = 0.75,
               plot.cor = F,
               cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                        '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                        '12' = c25[13])
) +
  facet_wrap(~ sc_harmony$orig.ident)

rownames(sc_harmony)[grep("^cd27", rownames(sc_harmony), ignore.case = T)]
VlnPlot(sc_harmony, features = c("Pdcd1", "Zap70", "Lck",
                                 "Il2", "Cd69", "Irf4",
                                 "Nfatc1", "Nfatc2", "Nfatc3"), 
        same.y.lims = T,
        split.plot = T,
        group.by = "seurat_clusters",
        # log = T,
        split.by = "orig.ident",
        pt.size = 0, 
        ncol = 3, 
        cols = c('0' = c25[1], '1' = c25[2], '2' = c25[3], '3' = c25[4], '4' = c25[5], '5' = c25[6],
                 '6' = c25[7], '7' = c25[8], '8' = c25[9], '9' = c25[10], '10' = c25[11],'11' = c25[12],
                 '12' = c25[13])
)

FeaturePlot(sc_harmony, 
            features = c(
              # "Pdcd1", "Zap70",
              "Ncr1", "Itga2", 
              "Cd27", "Klrg1"
            ),
            reduction = "umap", 
            pt.size = 0.25, 
            order = T, 
            split.by = "orig.ident", 
            min.cutoff = "q80")


# get median expression of the markers in each cluster
meta_median <- sc_all@meta.data
order_samples <- order(meta_median$seurat_clusters, decreasing = T)
meta_median <- meta_median[order_samples,]
head(meta_median)
median_markers <- sc_all[which(rownames(sc_all) %in% markers),]
median_markers <- data.matrix(median_markers@assays$RNA@counts[,order_samples])
# median_markers <- log1p(median_markers)
head(colnames(median_markers))
median_markers <- aggregate(t(median_markers), by=list(sc_all@meta.data$seurat_clusters[order_samples]), FUN=mean)
rownames(median_markers) <- median_markers$Group.1
head(median_markers)
median_markers$Group.1 <- NULL
heat <- apply(data.matrix(median_markers[,-1]), 2, scale)
# rownames(heat) <- order_unique_samples
rownames(heat) <- 1:nrow(heat)
Heatmap(t(heat), 
        cluster_rows = F, 
        cluster_columns = F, 
        row_dend_side = "right",row_names_side = "left",
        col = colorRamp2(c(-2, 0, 2), c("#437EA4","white","#E41A1C")))






# >> cell classification << ####
orTab <- singleCellNet::utils_loadObject("~/Google Drive/Protocols_bfx/human_mouse_genes_Jul_24_2018.rda")
# write.table(orTab, "~/Google Drive/Protocols_bfx/Ortholog_Human_Mouse_genes.txt", sep = "\t", row.names = F, quote = F)
# orTab <- read.table("~/Google Drive/Protocols_bfx/Ortholog_Human_Mouse_genes.txt", sep = "\t", header = T)
dim(orTab);head(orTab)
orTab$human <- as.character(orTab$human)
orTab$mouse <- as.character(orTab$mouse)
orTab[which(orTab$human == "CD45RA"),]
orTab[grep("ptprc", orTab$human, ignore.case = T),]

# run singleR for immune populations ####
library(SingleR)
# others to try are Single Cell Net, ClustifyR, CellAssign
ref_singler <- MouseRNAseqData(ensembl = FALSE, cell.ont = "none")
table(ref_singler$label.main)
ref_immgen <- ImmGenData(ensembl = FALSE, cell.ont = "nonna")
table(ref_immgen$label.fine,ref_immgen$label.main)
table(ref_immgen$label.main)
ref_immgen$label.main <- ifelse(ref_immgen$label.main == "NKT", "NK cells", ref_immgen$label.main)

counts_singlr <- LogNormalize(GetAssayData(sc_all, slot = "counts", assay = "RNA"))
counts_singlr[1:10,1:5]
tt_cells <- SingleR(test = counts_singlr,
                    ref = list(MD = ref_singler),
                    labels = list(ref_singler$label.main), 
                    de.method = "wilcox")
tt_cells <- SingleR(test = counts_singlr,
                    ref = list(MD = ref_immgen),
                    labels = list(ref_immgen$label.main))
head(tt_cells)

res_singler <- data.frame(IDs=rownames(tt_cells),
                          SingleR = tt_cells$pruned.labels)
res_singler$Cluster <- Idents(sc_all_eq)[match(res_singler$IDs, names(Idents(sc_all_eq)))]
res_singler$IDs_Seurat <- names(Idents(sc_all_eq))[match(res_singler$IDs, names(Idents(sc_all_eq)))]
res_singler$SingleR <- ifelse(res_singler$SingleR == "B cells, pro", "B cells",res_singler$SingleR)
tt_singler_res <- as.data.frame.matrix(table(res_singler$SingleR, res_singler$Cluster))
tt_singler_res <- tt_singler_res[-c(which(rowSums(tt_singler_res) <= 5),
                                    which(rownames(tt_singler_res) == "Microglia")),]
Heatmap(tt_singler_res,cluster_rows = T, cluster_columns = F,
        col = colorRamp2(c(0, 7, 15), c("#437EA4","white","#E41A1C")), name = "Abundance")

res_singler$GROUP <- sc_all@meta.data$GROUP[match(res_singler$IDs, names(Idents(sc_all)))]
table(res_singler$GROUP, res_singler$SingleR)
head(res_singler)


plotDeltaDistribution(tt_cells, 
                      labels.use = rownames(tt_singler_res),
                      size = 1, 
                      ncol = 4)
head(tt_cells)

# calculate percentage cell types from singleR ####
perc_table <- data.frame(t(apply(tt_singler_res, 2, function(xx) xx/sum(xx[1], xx[2]))))
rws <- rownames(perc_table)
perc_table <- melt(perc_table)
colnames(perc_table) <- c("GROUP", "CellsFraction")
perc_table$Gene <- rep(rws, 2)
ggplot(perc_table, aes(x=Gene, y=CellsFraction, fill=GROUP)) +
  geom_bar(stat="identity") +
  theme_bw() + xlab("") +
  scale_fill_manual(values = c("C" = alpha("blue", 0.75), "T" = alpha("red", 0.75))) +
  theme(axis.text.x = element_text(size=10, angle=70, vjust = 0.5),
        axis.title.y = element_text(size=10))

cell_for_seurat <- res_singler$SingleR[match(rownames(sc_all@meta.data), res_singler$IDs_Seurat)]
cell_for_seurat <- ifelse(cell_for_seurat == "B cells, pro", "B cells", cell_for_seurat)
cell_for_seurat <- ifelse(cell_for_seurat == "NKT", "NK cells", cell_for_seurat)
cell_for_seurat <- ifelse(cell_for_seurat == "Adipocytes", "Stromal cells", cell_for_seurat)
cell_for_seurat <- ifelse(cell_for_seurat == "Dendritic cells", "DC", cell_for_seurat)
cell_for_seurat <- ifelse(cell_for_seurat %in% rownames(tt_singler_res), cell_for_seurat, NA)
sc_all@meta.data$singleR <- cell_for_seurat
table(cell_for_seurat)
DimPlot(sc_all, reduction = "umap", group.by = "singleR",
        cols = c('B cells' = c25[1], 'DC' = c25[2], 'Fibroblasts' = c25[3], 'ILC' = c25[4], 'Macrophages' = c25[5], 'Monocytes' = c25[6],
                 'Neutrophils' = c25[7], 'NK cells' = c25[8], 'Stromal cells' = c25[9], 'T cells' = c25[10], "Endothelial cells"=c25[11],
                 "Erythrocytes" = c25[12], "Granulocytes" = c25[13]),
        split.by = "GROUP")

# calculated differential gene expresison from singleR ####

DefaultAssay(sc_all_eq) <- "RNA"
cell_type <- rownames(heat)
for(ii in 1:length(cell_type)){
  if(cell_type[ii] %in% c("Stromal cells", "Fibroblasts"))
    next
  
  temp.cells <- sc_all_eq
  Idents(temp.cells) <- temp.cells@meta.data$singleR
  temp.cells <- subset(temp.cells, idents = cell_type[ii])
  
  Idents(temp.cells) <- as.character(temp.cells@meta.data$GROUP)
  table(Idents(temp.cells))
  
  # find DEGs
  temp.res <- FindMarkers(temp.cells, 
                          ident.1 = "T", ident.2 = "C", 
                          logfc.threshold = 0.1, min.pct = 0.1,
                          test.use = "MAST",  
                          assay = "RNA")
  head(temp.res, n=20)
  temp.res$Cell_type <- cell_type[ii]
  temp.res$Gene <- rownames(temp.res)
  
  if(ii == 1){
    res_cluster_deg_cells <- temp.res
    
  } else {
    res_cluster_deg_cells <- data.frame(rbind(res_cluster_deg_cells, temp.res))
    
  }
  message(paste("Done for",cell_type[ii]))
  
}
head(res_cluster_deg_cells); tail(res_cluster_deg_cells)
res_cluster_deg_cells_sig <- res_cluster_deg_cells[which(res_cluster_deg_cells$p_val < 0.1),]

# run GSEA as in Gobin's paper ####
table(res_new_deg_all$Cell_type)
res_new_deg_all$Cell_type <- ifelse(res_new_deg_all$Cell_type == "M2.1", "M2",
                                    ifelse(res_new_deg_all$Cell_type == "M2.2", "M1.2", res_new_deg_all$Cell_type))
# take top 100 up/down per each cluster
# temp_mono <- res_cluster_deg_all[-which(res_cluster_deg_all$CellType %in% c("CD8 T", "CD4 T", "Mono", "NK","DC", "other","B")),]
head(res_new_deg_all)
temp_mono <- res_new_deg_all
for(ii in 1:length(unique(temp_mono$Cell_type))){
  # for(ii in 1:seq(1,24,4)){
  if(ii == 1){
    top100_genes <- data.frame(matrix(ncol = length(unique(temp_mono$Cell_type))*4,
                                      nrow=500))
    colnames(top100_genes) <- c(sapply(unique(temp_mono$Cell_Type), function(xx) paste0(xx,c("_UP_Gene","_UP_FC","_DOWN_Gene","_DOWN_FC"))))
  }
  temp_deg <- temp_mono[which(temp_mono$Cell_type == unique(temp_mono$Cell_type)[ii]),]
  # temp_deg <- temp_deg[which(temp_deg$p_val_adj <= 0.1),]
  temp_deg <- temp_deg[order(temp_deg$avg_log2FC, decreasing = T),]
  sup <- sum(temp_deg$avg_log2FC > 0)
  temp_names_up <- temp_deg$Gene[1:500]
  temp_up <- temp_deg$avg_log2FC[1:500]
  temp_deg <- temp_deg[order(temp_deg$avg_log2FC, decreasing = F),]
  sdown <- sum(temp_deg$avg_log2FC < 0)
  temp_names_down <- temp_deg$Gene[1:500]
  temp_down <- temp_deg$avg_log2FC[1:500]
  
  if(ii == 1){
    top100_genes <- data.frame(cbind(temp_names_up,
                                     temp_up,
                                     temp_names_down,
                                     temp_down))
    
  } else {
    top100_genes <- data.frame(cbind(top100_genes,
                                     temp_names_up,
                                     temp_up,
                                     temp_names_down,
                                     temp_down))
  }
  if(ii == length(unique(temp_mono$Cell_type))){
    colnames(top100_genes) <- c(sapply(unique(temp_mono$Cell_type), function(xx) paste0(xx,c("_UP_Gene","_UP_FC","_DOWN_Gene","_DOWN_FC"))))
  }
  # top100_genes[, ii] <- temp_names_up
  # top100_genes[, ii+1] <- temp_up
  # top100_genes[, ii+2] <- temp_names_down
  # top100_genes[, ii+3] <- temp_down
  
  message(paste("Done for", unique(temp_mono$Cell_type)[ii]))
}

# get gene signatures
gmt_custom <- c(qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/MousePath_GO_gmt.gmt"),
                qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/MousePath_Pathway_gmt.gmt"))
gmt_custom <- gmt_custom[grep("^GO_BP_", names(gmt_custom))]
gmt_custom <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/MousePath_Pathway_gmt.gmt")
gmt_custom <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/h.all.v6.2.symbols.gmt")
gmt_custom <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.reactome.v6.1.symbols.gmt")
gmt_custom <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c2.cp.kegg.v6.1.symbols.gmt")
gmt_custom <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/c5.bp.v6.1.symbols.gmt")
gmt_custom <- data.frame(ont = rep(names(gmt_custom), as.numeric(lapply(gmt_custom, length))),
                         gene = unlist(gmt_custom))
head(gmt_custom)

for(ii in 1:length(unique(res_new_deg_all$Cell_type))){
  cls_1 <- paste0(unique(res_new_deg_all$Cell_type)[ii], "_UP_FC")
  fc_1 <- top100_genes[,which(colnames(top100_genes) == cls_1)]
  cls_2 <- paste0(unique(res_new_deg_all$Cell_type)[ii], "_DOWN_FC")
  fc_2 <- top100_genes[,which(colnames(top100_genes) == cls_2)]
  
  cls_1 <- paste0(unique(res_new_deg_all$Cell_type)[ii], "_UP_Gene")
  gene_1 <- top100_genes[,which(colnames(top100_genes) == cls_1)]
  cls_2 <- paste0(unique(res_new_deg_all$Cell_type)[ii], "_DOWN_Gene")
  gene_2 <- top100_genes[,which(colnames(top100_genes) == cls_2)]
  
  ranked_genes <- as.numeric(c(fc_1, rev(fc_2)))
  names(ranked_genes) <- c(gene_1, rev(gene_2))
  fgsea_up <- fgsea::fgsea(pathways = gmt_custom,
                           stats = ranked_genes)
  if(ii == 1){
    list_res <- list()
    res_gsea_df <- data.frame()
  } 
  list_res[[ii]] <- fgsea_up[which(fgsea_up$padj < 0.1),c(1,3,6)]
  
  fgsea_up$CellType <- unique(res_new_deg_all$Cell_type)[ii]
  res_gsea_df <- data.frame(rbind(res_gsea_df, fgsea_up[which(fgsea_up$padj < 0.1),c(1,3,6,9)]))
  
  message(nrow(list_res[[ii]]))
  names(list_res)[ii] <- unique(res_new_deg_all$Cell_type)[ii]
  
  message(paste("Done for",unique(res_new_deg_all$Cell_type)[ii]))
}

View(res_gsea_df %>% filter(padj <= 0.1) %>% filter(grepl("immune", pathway, ignore.case = T)))
# write.table(res_gsea_df, "res_GSEA_all_01.txt", sep='\t', row.names = F, quote=F)
res_gsea_df <- read.table("GSEA_significant_pathways.txt", header = T, sep="\t")
res_gsea_df_sig <- res_gsea_df %>% filter(padj <=0.05)
length(unique(res_gsea_df_sig$pathway))

to_plot <- res_gsea_df_sig[unique(c(grep("immun", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("antigen", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("chemotaxis", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("chemok", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("cytok", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("migr", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("death", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("prolif", res_gsea_df_sig$pathway,ignore.case = T),
                                    grep("interleu", res_gsea_df_sig$pathway,ignore.case = T))),]



names(list_res)
list_res[[1]]
lapply(list_res, print)


fgsea_all <- data.frame(matrix(nrow=length(unique(to_plot$pathway)),
                               ncol=length(unique(to_plot$CellType))))
rownames(fgsea_all) <- unique(to_plot$pathway)
colnames(fgsea_all) <- unique(to_plot$CellType)


for(ii in 1:nrow(fgsea_all)){
  fgsea_all[ii,] <- rep(0,ncol(fgsea_all))
  
  temp_path <- to_plot[which(to_plot$pathway == rownames(fgsea_all)[ii]),]
  
  ind_cols <- which(colnames(fgsea_all) %in% unique(temp_path$CellType))
  
  if(length(ind_cols) > 1){
    temp_scores <- temp_path$NES[match(colnames(fgsea_all), temp_path$CellType)]
    temp_scores <- temp_scores[-which(is.na(temp_scores))]
    fgsea_all[ii,ind_cols] <- temp_scores
  } else {
    fgsea_all[ii,ind_cols] <- temp_path$NES
  }
  # colnames(fgsea_all)
  
}
fgsea_all$Pathway <- rownames(fgsea_all)
fgsea_all$Pathway[3] <- "CLASS_II_ANTIGEN_PROC_PRESENTATION"
fgsea_all$Pathway[12] <- "POS_REG_NK_CELL_CHEMOTAXIS"
fgsea_all$Pathway[22] <- "POS_REG_CELL_MIGRATION_ANGIOGENESIS"
fgsea_all$Pathway[16] <- "NEG_REG_CHEMOKINE_MEDIATED_SIGNALING"
fgsea_all$Pathway[17] <- "POS_REG_CHEMOKINE_CXC_LIGAND"
fgsea_all$Pathway[23] <- "NEG_REG_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION"
fgsea_all$Pathway[24] <- "POS_REG_SMOOTH_MUSCLE_CELL_MIGRATION"




mm <- melt(fgsea_all)
colnames(mm)[2:3] <- c("CellType", "NES")
head(mm)
mm$Pathway <- gsub("GO_BP_MM_", "", 
                   gsub("GO_CC_MM_", "", 
                        gsub("GO_MF_MM_", "", 
                             gsub("PANTHER_MM_", "", 
                                  gsub("WIKIPATHWAYS_MM_", "",
                                       gsub("INOH_MM_", "", 
                                            gsub("NETPATH_MM_","",
                                                 sub("BIOCARTA_MM_", "", mm$Pathway))))))))


mm$Direction <- ifelse(mm$NES > 0, "UP",
                       ifelse(mm$NES < 0, "DOWN","NoChange"))
head(mm)
ggplot(mm, aes(x=CellType, y=Pathway, col=Direction, size=abs(NES))) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = c("UP" = "#DC143C", "DOWN" = "#4682B4", "NoChange" = "white"),
                     breaks = c("UP", "DOWN")) +
  scale_size(breaks = c(0,1,2)) +
  guides(size = guide_legend(title = "NES"),
         color = guide_legend(override.aes = list(size=5))) +
  theme(axis.text.x = element_text(size=14,angle = 70,hjust = 1),
        axis.text.y = element_text(size=11)) +
  ylab("") + xlab("")

# > address Jonathan's questions ####
# Can you help make a new figure to just like the above, with the same UMAP. But say e.g. gray all the cell events. 
# And show these 2 genes TRBV13-2 and TRBV13-3 (or just TRBV13, if those are not available), over the gray background.
# This way I can see where these 2 genes show up in which cell type population.
# Maybe have 2 graphs, 1 for vehicle and 1 for treatment, to see the changes as well.
sc_harmony
table(sc_harmony$orig.ident)
rownames(sc_harmony)[grep("TRBV13", rownames(sc_harmony), ignore.case = T)]
FeaturePlot(sc_harmony,
            col=c("lightgray", "red"),
            order = T,
            min.cutoff = "q75") +
  facet_wrap(~ sc_harmony$orig.ident)

FeaturePlot(sc_harmony,
            features = "Trbv13-3",
            col=c("lightgray", "red"),
            order = T,
            min.cutoff = "q75") +
  facet_wrap(~ sc_harmony$orig.ident)

sc_harmony$DP <- sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",] + sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",]
# sc_harmony$DP <- (sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",]+0.1) * (sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",]+0.1)
FeaturePlot(sc_harmony,
            features = "DP",
            col=c("lightgray", "red"),
            # order = T, 
            min.cutoff = "q15") +
  facet_wrap(~ sc_harmony$orig.ident) +
  ggtitle("Trbv13-2 + Trbv13-3")


FeaturePlot(sc_harmony,
            cells = rownames(sc_harmony@meta.data)[which(sc_harmony$orig.ident == "Treated")],
            features = c("Trbv13-2", "Trbv13-3"),
            blend = T,
            col=c("red", "green"),
            order = T,
            min.cutoff = "q75"
)
# +
#   facet_wrap(~ sc_harmony$orig.ident)
# get percentage populations:
sc_harmony@meta.data <- read.table("sc_metadata_annotated.txt", sep='\t', header = T, row.names = 1)
tt <- as.data.frame.matrix(table(sc_harmony$TILPRED_harmonized.2, sc_harmony$orig.ident))
tt <- data.frame(cbind(tt,
                       apply(tt, 2, function(xx) round(xx/sum(xx), 4)*100)))
tt$CellType <- rownames(tt)

colnames(tt)[3:4] <- gsub("\\.1","",colnames(tt)[3:4])
order_rows <- tt$CellType[order(tt$Treated.1, decreasing = T)]
mm <- reshape2::melt(tt[,3:5])
colnames(mm)[2:3] <- c("Group", "Percentage")
mm$Group <- gsub("\\.1","",mm$Group)
head(mm)
ggplot(mm, aes(x=factor(CellType, levels = rev(order_rows)), y=Percentage, fill=factor(Group, levels=c("CTR", "Treated")))) +
  geom_bar(stat="identity", position="dodge", col="black") +
  theme_classic() +
  coord_flip() +
  labs(fill="Group") +
  ylab("Percentage Cells") + xlab("") +
  scale_fill_manual(values = c("CTR" = scales::alpha("gray", 0.5),
                               "Treated"= scales::alpha("orange", 0.5)))

tt[,3:4] %>%
  arrange(desc(Treated)) %>%
  kbl() %>%
  kable_styling(font_size = 20)

Idents(sc_harmony) <- sc_harmony$TILPRED_harmonized.2
DoHeatmap(sc_harmony, 
          slot = "data", 
          features = c("Trbv13-2", "Trbv13-3"), 
          cells = rownames(sc_harmony@meta.data)[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory")])
colnames(gg_t)[grep("trbv13",colnames(gg_t), ignore.case = T)]
cell_line <- c("CD8_Effector_Memory")
ggplot(gg_t[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory"),], 
       aes(x=(Trbv13.2 + Trbv13.3), 
           col=factor(sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory")]))) +
  stat_ecdf(geom = "step", size=1.5) +
  theme_bw() +
  labs(col = "") +
  ylab("cumulative frequency") +
  ggtitle(cell_line) +
  # ggtitle("CD4/Treg") +
  scale_color_manual(values = c("CTR" = "blue",
                                "Treated" = "tomato"))

ggplot(gg_t[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory"),],
       aes(x=Trbv13.2+Trbv13.3, y=Cd8b1, col=Trbv13.2)) +
  geom_point() +
  theme_classic() +
  scale_color_gradient2(low=muted("blue"), 
                        mid = "lightgray", 
                        high = muted("tomato"), 
                        midpoint = (median(gg_t$Trbv13.2)+median(gg_t$Trbv13.3))) +
  labs(color="Expression") +
  facet_wrap(~ sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector_Memory")])




VlnPlot(sc_harmony,
        idents = "CD8_Effector_Memory",
        features = c("Trbv13-2", "Trbv13-3", "Ccl5"), same.y.lims = T,
        split.by = "orig.ident",
        pt.size = 0)
table(sc_harmony$TILPRED_harmonized.2, sc_harmony$orig.ident)

# address other questions for Jonathan ####
FeatureScatter(sc_harmony, 
               cells = colnames(sc_harmony)[-which(sc_harmony$TILPRED_harmonized.2 == "Monocytes")],
               group.by = "orig.ident",
               feature1 = "Entpd1",
               feature2 = "Itgae",
               cols=c25) +
  facet_wrap( ~ sc_harmony$TILPRED_harmonized.2[-which(sc_harmony$TILPRED_harmonized.2 == "Monocytes")])

ggplot(gg_t[-which(sc_harmony$TILPRED_harmonized.2 == "Monocytes"),], 
       aes(x=(Entpd1 * Itgae * Tnfrsf9), 
           col=factor(sc_harmony$orig.ident[-which(sc_harmony$TILPRED_harmonized.2 == "Monocytes")]))) +
  stat_ecdf(geom = "step", size=1.5) +
  theme_bw() +
  labs(col = "") +
  ylab("cumulative frequency") +
  ggtitle(cell_line) +
  # ggtitle("CD4/Treg") +
  scale_color_manual(values = c("CTR" = "blue",
                                "Treated" = "tomato")) +
  facet_wrap(~ sc_harmony$orig.ident[-which(sc_harmony$TILPRED_harmonized.2 == "Monocytes")])

ggplot(NULL, 
       aes(x=sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector")], 
           y=sc_harmony@assays$SCT@data["Entpd1",which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector")])) +
  geom_boxplot() +
  theme_bw() +
  # labs(col = "") +
  ylab("Expression") + xlab("Groups") 

t.test(gg_t[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector"),"Itgae"] ~ sc_harmony$orig.ident[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector")])


Idents(sc_harmony) <- sc_harmony$TILPRED_harmonized.2
VlnPlot(sc_harmony, 
        group.by = "TILPRED_harmonized.2",
        feature = "Entpd1",
        idents = "CD8_Effector",
        split.by = "orig.ident",
        # cols = c25,
        cols = c("blue", "red"),
        pt.size = 0.5)
temp_meta <- read.table("test_meta_for_celltype_test_2_sample_rdata.txt", header = T, sep = "\t", row.names = 1)
temp_meta2 <- read.table("sc_metadata_annotated.txt", header = T, sep = "\t", row.names = 1)
Heatmap(as.data.frame.matrix(table(sc_harmony$TILPRED_harmonized, sc_harmony$TILPRED_harmonized.2)),
        col = colorRamp2(breaks = c(0,10,100,300),
                         colors = c("white", "cyan", "pink", "purple")),
        border = T)


# try visualization tools ####

# Create Vitessce view config
vc <- VitessceConfig$new("My config")
dataset <- vc$add_dataset("My dataset")$add_object(SeuratWrapper$new(sc_harmony, cell_set_meta_names = list("TILPRED_harmonized.2"), out_dir = "out"))
scatterplot <- vc$add_view(dataset, Component$SCATTERPLOT, mapping = "umap")
status <- vc$add_view(dataset, Component$STATUS)
desc <- vc$add_view(dataset, Component$DESCRIPTION)
desc <- desc$set_props(description = "Mouse scRNASeq data")
cell_sets <- vc$add_view(dataset, Component$CELL_SETS)
Genes <- vc$add_view(dataset, Component$GENES)
Spatial <- vc$add_view(dataset, Component$SPATIAL)
heatmap <- vc$add_view(dataset, Component$HEATMAP)
vc$layout(hconcat(
  vconcat(scatterplot, heatmap),
  vconcat(cell_sets, vconcat(desc, status)),
  vconcat(Genes, Spatial)
))

# Render the Vitessce widget
vc$widget(theme = "light")

# >> ? address Jonathan's questions - 01.10.2022 ####
# check GSSE results
# check CD8_effector_memory_exhausted results
head(res_all_fgsea$CD8_memo)
# >> this was done above, problem fixed with the analysis
## get avg expression by cell type
genes2check <- c("Il2", "Il2ra", "Ifng", "Gzma", "Gzmb", "Foxp3", "Plac8")
df_out <- t(sc_harmony@assays$SCT@data[which(rownames(sc_harmony) %in% genes2check),])
head(sc_harmony@meta.data)
colnames(sc_harmony@meta.data)
df_out <- merge(df_out, sc_harmony@meta.data[,c(1,28)], by="row.names")
colnames(df_out)[1] <- "Barcode"
head(df_out)
# write.table(df_out, "Expression_selected_genes_by_cell_type_all.txt", sep='\t', quote = F, row.names = F)
df_out <- aggregate(df_out[,2:8], by=list(ID=paste0(df_out$TILPRED_harmonized.2, "_", df_out$orig.ident)), FUN=mean)
# write.table(df_out, "Expression_selected_genes_by_cell_type_treatment_mean.txt", sep='\t', quote = F, row.names = F)

# violin plots of the above genes in selected populations
pop2use <- c("CD8_Effector", "CD8_Memory_Exhausted", "CD8_Effector_Memory", "CD8_Naive_CM", "CD4_Th1_Helper", "CD4_Treg_Effector", "CD8_CM", "CD4_CM")
pop2use <- c("CD8_CM")
VlnPlot(sc_harmony %>%
          subset(cells = colnames(sc_harmony)[which(sc_harmony$TILPRED_harmonized.2 %in% pop2use)]), 
        assay = "SCT", 
        group.by = "TILPRED_harmonized.2",
        features = "Foxp3", #genes2check,
        split.by = "orig.ident",
        # stack = T,
        # col=c25,
        pt.size = 0) +
  ylab("")
table(sc_harmony$TILPRED_harmonized.2, sc_harmony$orig.ident)

DefaultAssay(sc_harmony) <- "MAGIC_RNA"

FeatureScatter(sc_harmony %>%
                 subset(cells = colnames(sc_harmony)[which(sc_harmony$TILPRED_harmonized.2 == "CD8_Effector")]), 
               feature1 = "Cd8a",
               feature2 = "Cd4",
               group.by = "orig.ident",
               plot.cor = F)

DimPlot(sc_harmony, 
        group.by = "TILPRED_harmonized.2", 
        split.by = "orig.ident", 
        # label = T, label.box = T,repel = T,
        cells.highlight = colnames(sc_harmony)[which(sc_harmony$TILPRED_harmonized.2 == "CD8_CM")]) 

DimPlot(sc_harmony, 
        group.by = "TILPRED_harmonized.2", 
        split.by = "orig.ident", 
        # label = T, label.box = T,repel = T,
        cells.highlight = colnames(sc_harmony)[which(sc_harmony$TILPRED_harmonized.2 == "CD4_Treg_Effector")]) 

FeaturePlot(sc_harmony,
            # features = c("Il2ra", "Gzma", "Prf1", "Foxp3", "Ifng"),
            features = c("Cd8a", "Cd4"),
            # features = c("Gzmb", "Zap70", "Lck", "Lamp1", "Klf2"),
            split.by = "orig.ident",
            min.cutoff = "q50", 
            ncol=4)

res_all %>%
  filter(CellType == "CD8_Effector") %>%
  filter(Significant == "YES") %>%
  arrange(avg_log2FC) %>%
  slice(seq_len(20)) %>%
  ggplot(., aes(x=factor(Gene, levels = as.character(Gene)), y=avg_log2FC)) +
  geom_bar(stat="identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=12)) +
  xlab("")
View

# extract pathway data for Jonathan - email 02.01.2022####
path2get <- c("WIKIPATHWAYS_MM_ELECTRON_TRANSPORT_CHAIN-WP111", "BIOCARTA_MM_TH1_TH2_DIFFERENTIATION", "BIOCARTA_MM_IL_17_SIGNALING_PATHWAY", "NETPATH_MM_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "BIOCARTA_MM_SELECTIVE_EXPRESSION_OF_CHEMOKINE_RECEPTORS_DURING_T-CELL_POLARIZATION", "BIOCARTA_MM_LCK_AND_FYN_TYROSINE_KINASES_IN_INITIATION_OF_TCR_ACTIVATION", "BIOCARTA_MM_CTL_MEDIATED_IMMUNE_RESPONSE_AGAINST_TARGET_CELLS_", "INOH_MM_JAK_STAT_MOLECULARVARIATION_2", "INOH_MM_CD4_T_CELL_RECEPTOR_SIGNALING-ERK_CASCADE", "BIOCARTA_MM_ROLE_OF_TOB_IN_T-CELL_ACTIVATION", "INOH_MM_JAK_STAT_MOLECULARVARIATION_1", "BIOCARTA_MM_T_HELPER_CELL_SURFACE_MOLECULES", "BIOCARTA_MM_IL12_AND_STAT4_DEPENDENT_SIGNALING_PATHWAY_IN_TH1_DEVELOPMENT", "WIKIPATHWAYS_MM_TCR_SIGNALING_PATHWAY-WP69", "BIOCARTA_MM_T_CYTOTOXIC_CELL_SURFACE_MOLECULES", "BIOCARTA_MM_THE_CO-STIMULATORY_SIGNAL_DURING_T-CELL_ACTIVATION")
load("~/Google Drive/Protocols_bfx/GSEA_signatures/mouse_c2_v5p2.rdata")
load("~/Google Drive/Protocols_bfx/GSEA_signatures/Mm.c3.tft.v6.1.entrez.RData")
load("~/Google Drive/Protocols_bfx/GSEA_signatures/Mm.c2.cp.reactome.v6.1.entrez.RData")
paths <- qusage::read.gmt("~/Google Drive/Protocols_bfx/GSEA_signatures/MousePath_Pathway_gmt.gmt")
head(names(Mm.c2))
path2get[-which(path2get %in% names(paths))]
names(paths)[grep("biocarta", names(paths), ignore.case = T)]
lapply(paths[grep("biocarta", names(paths), ignore.case = T)], 
       function(xx) cat, 
       "\n", 
       file="~/Google Drive/BridgeBfx_SB/MarengoTx/pathways_for_jonathan.txt", 
       append=TRUE)
sink(file="~/Google Drive/BridgeBfx_SB/MarengoTx/pathways_for_jonathan2.txt")
paths[which(names(paths) %in% path2get)]
sink()

# address TRMs from email 02.18.2022 ####
rownames(sc_harmony)[grep("cd28", rownames(sc_harmony), ignore.case = T)]
VlnPlot(sc_harmony,
        assay = "MAGIC_RNA",
        features = c("Cd8a", "Cd4", "Cd69", "Ccr7", "Sell", "Itgae","Itga1", "Ifng", "Cd44", "Cxcr6", "Icos", "Cd28", "Pdcd1", "Il2ra", "Cx3cr1", "Klf2"),
        # group.by = "TILPRED_harmonized.2",
        group.by = "seurat_clusters",
        # split.by = "orig.ident",
        stack = T, flip = T,
        col=c25,
        pt.size = 0)

DefaultAssay(sc_harmony) <- "MAGIC_RNA"

# make biplots to label TRMs independently from current annotation
FeatureScatter(sc_harmony,
               feature1 = "Cd69",
               feature2 = "Itgae",
               pt.size = 0.5,
               cols = c25,
               # group.by = "TILPRED_harmonized.2",
               group.by = "seurat_clusters",
               plot.cor = F) +
  # facet_wrap(~ sc_harmony$orig.ident, ncol=2) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2))

sc_harmony$Cd69p_Cd103p  <- ifelse(sc_harmony@assays$MAGIC_RNA@data["Cd69",] > 1 &
                                     sc_harmony@assays$MAGIC_RNA@data["Itgae",] > 0.5, "YES", "NO")
VlnPlot(sc_harmony,
        features = c("Itga1", "Cd44"),
        group.by = "Cd69p_Cd103p",
        pt.size = 0)

FeatureScatter(sc_harmony,
               feature1 = "Cd4",
               feature2 = "Cd44",
               pt.size = 0.5,
               cols = c25,
               # group.by = "TILPRED_harmonized.2",
               group.by = "orig.ident",
               plot.cor = F) +
  facet_wrap(~ sc_harmony$Cd69p_Cd103p, ncol=2) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2))

sc_harmony$Cd8TRM  <- ifelse(sc_harmony$Cd69p_Cd103p == "YES" &
                               sc_harmony@assays$MAGIC_RNA@data["Cd44",] > 0.7 &
                               sc_harmony@assays$MAGIC_RNA@data["Cd8a",] > 0.7, "YES", "NO")

sc_harmony$Cd4TRM  <- ifelse(sc_harmony$Cd69p_Cd103p == "YES" &
                               sc_harmony@assays$MAGIC_RNA@data["Cd44",] > 0.7 &
                               sc_harmony@assays$MAGIC_RNA@data["Cd4",] > 0.7, "YES", "NO")

table(sc_harmony$Cd8TRM,
      sc_harmony$Cd4TRM)

sc_harmony$CellType_manual <- ifelse(sc_harmony$Cd8TRM == "YES", "CD8_TRM",
                                     ifelse(sc_harmony$Cd4TRM == "YES", "CD4_TRM", sc_harmony$TILPRED_harmonized.2))
DimPlot(sc_harmony,
        group.by = "CellType_manual",
        cols = c(c25, "cyan", muted("red")),
        split.by = "orig.ident")

FeaturePlot(sc_harmony,
            features = c("Cd8a", "Cd4"),
            pt.size = 0.5,
            min.cutoff = "q25")
FeatureScatter(sc_harmony,
               feature1 = "Cd4",
               feature2 = "Cd8a",
               pt.size = 0.5,
               cols = c25,
               # group.by = "TILPRED_harmonized.2",
               group.by = "seurat_clusters",
               plot.cor = F) +
  # facet_wrap(~ sc_harmony$Cd8TRM, ncol=2) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2))

rownames(sc_harmony)[grep("TRBV", rownames(sc_harmony), ignore.case = T)]
FeaturePlot(sc_harmony,
            features = c( "Trbv13-3"),
            # min.cutoff = "q15",
            pt.size = 0.5) +
  facet_wrap(~ sc_harmony$orig.ident, ncol=2)

# address Jonathan email 02242022 ####
# plot all TRBVs
trbv2check <- c("Trbv1", "Trbv12-1", "Trbv12-2", "Trbv13-1", "Trbv13-2", "Trbv13-3", "Trbv14", "Trbv15", "Trbv16", "Trbv17", "Trbv19", "Trbv2", "Trbv20", "Trbv21", "Trbv22", "Trbv23", "Trbv24", "Trbv26", "Trbv29", "Trbv3", "Trbv30", "Trbv31", "Trbv4", "Trbv5")
for(ii in 1:length(trbv2check)){
  if(ii == 1)
    pdf("UMAP_TRBVs_02242022.pdf", width = 6.5, height = 3.25)
  
  if(length(which(rownames(sc_harmony) == trbv2check[ii])) == 1){
    gg <- FeaturePlot(sc_harmony,
                      features = trbv2check[ii],
                      min.cutoff = "q15",
                      cols = c("lightgray", "navy"),
                      # order = T,
                      pt.size = 0.25) +
      facet_wrap(~ sc_harmony$orig.ident, ncol=2)
    print(gg)
    message(paste("Done for", trbv2check[ii]))
  } else {
    message(paste("No match for", trbv2check[ii]))
  }
}
dev.off()

# combine TRBV13-2 + TRBV13-3.
ggplot(sc_harmony@meta.data %>%
         mutate(Trbv13sum = sc_harmony@assays$SCT@data["Trbv13-2",] +
                  sc_harmony@assays$SCT@data["Trbv13-3",]) %>%
         arrange(Trbv13sum),
       aes(x=sc_harmony@reductions$umap@cell.embeddings[,1],
           y=sc_harmony@reductions$umap@cell.embeddings[,2],
           col=Trbv13sum)) +
  # arrange() +
  geom_point(size = 0.25) +
  theme_classic() +
  scale_color_gradient2(low = "lightgray", mid = "lightgray", high = "tomato") +
  facet_wrap(~ orig.ident) +
  xlab("UMAP_1") + ylab("UMAP_2") +
  labs(col=element_blank())

# other list of genes in email
rownames(sc_harmony)[grep("^il27",rownames(sc_harmony),ignore.case = T)]
genes2check <- c("Il7r", "Sell", "Ccr7", "Ifng", "Gzma", "Gzmb", "Foxp3", "Ctla4", "Cd69", "Fas", "Entpd1", "Itgae", "Tox", "Tcf7", "Eomes", "Nkg7", "Tbx21", "Prdm1", "Pdcd1", "Lag3", "Havcr2", "Klrg1", "Tigit", "Cd244", "Ctla4", "Cd28", "Tnfrsf9", "Cd40lg", "Cxcr4", "Ccr5", "Cxcl13", "Il2ra")
genes2check[-which(genes2check %in% rownames(sc_harmony))]

for(ii in 1:length(genes2check)){
  if(ii == 1)
    pdf("UMAP_genes_v2_02242022.pdf", width = 6.5, height = 3.25)
  
  if(length(which(rownames(sc_harmony) == genes2check[ii])) == 1){
    gg <- FeaturePlot(sc_harmony,
                      features = genes2check[ii],
                      min.cutoff = "q15",
                      cols = c("lightgray", "navy"),
                      # order = T,
                      pt.size = 0.25) +
      facet_wrap(~ sc_harmony$orig.ident, ncol=2)
    print(gg)
    message(paste("Done for", genes2check[ii]))
  } else {
    message(paste("No match for", genes2check[ii]))
  }
}
dev.off()

# re-label cells ####
table(sc_harmony$TILPRED_harmonized.2)
sc_harmony$new_labels <- ifelse(sc_harmony$TILPRED_harmonized.2 %in% c("CD8_Cycling_Effector", "CD8_ShortLived_Effector"), "CD8_Effector",
                                ifelse(sc_harmony$TILPRED_harmonized.2 %in% c("CD8_CM","CD8_Early_Activated", "CD8_Precursor_Memory", "CD8_Intermediate_Memory"), "CD8_Naive_CM",
                                       ifelse(sc_harmony$TILPRED_harmonized.2 %in% c("CD4_Treg_Effector"), "CD4_Treg",
                                              ifelse(sc_harmony$TILPRED_harmonized.2 %in% c("CD4_Naive_CM"), "CD4_CM",
                                                     ifelse(sc_harmony$TILPRED_harmonized.2 %in% c("CD4_Th1_Helper", "CD4_Th1_Memory", "CD4_Follicular_Helper_Effector", "CD4_Follicular_Helper", "CD4_Follicular_Helper_IFN_Stim", "CD4_Follicular_Helper_Memory"), "CD4_T_helper",sc_harmony$TILPRED_harmonized.2)))))

table(sc_harmony$new_labels,
      sc_harmony$orig.ident)

DimPlot(sc_harmony,
        group.by = "new_labels",
        cols = c(c25, "cyan", muted("red")),
        split.by = "orig.ident")

VlnPlot(sc_harmony %>% subset(cells = which(sc_harmony$new_labels == "CD4_T_helper")), 
        # assay = "SCT",
        ncol=5, same.y.lims = T,
        group.by = "new_labels",
        features = c("Bcl6", "Rorc", "Tbx21", "Gata3", "Il18r1"),
        pt.size = 0.123)

FeaturePlot(sc_harmony,
            c("Bcl6", "Rorc", "Tbx21", "Gata3", "Il2rb", "Stat1"),
            min.cutoff = "q50",
            order=T,
            ncol=3)

Idents(sc_harmony) <- sc_harmony$new_labels
cell_markers_top <- FindAllMarkers(sc_harmony,
                                   logfc.threshold = 0.5,
                                   only.pos = T, 
                                   min.pct = 0.5)

head(cell_markers_top)
cell_markers_sel <- cell_markers_top %>%
  filter(!grepl("^Gm",gene)) %>%
  filter(!grepl("^Hist",gene)) %>%
  filter(p_val_adj <= 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(30))
write.table(cell_markers_top, "CellMarkers_newAnnotaiton_all_02252022.txt", sep='\t', quote = F, row.names = F)  


DoHeatmap(sc_harmony, 
          assay = "MAGIC_RNA", 
          slot = "data",
          features = cell_markers_sel$gene, draw.lines = T)

# make heatmap for columns
set.seed(100000)
sample_cols <- sample(1:ncol(sc_harmony), 2000, replace = F)
heat <- sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% cell_markers_sel$gene),
                                         sample_cols]
heat[1:10,1:10]
rws <- rownames(heat)
cls <- colnames(heat)
heat <- t(apply(heat, 1, scale))
colnames(heat) <- cls

df_anno <- sc_harmony@meta.data[sample_cols, c("new_labels", "orig.ident")]
head(df_anno)
df_anno <- df_anno[order(df_anno$new_labels),]

Heatmap(heat[,match(rownames(df_anno), colnames(heat))], 
        show_column_names = F,
        cluster_columns = F,
        top_annotation = df_anno)



# add the number/% of TRBV31-2+TRBV13-3 positive per each of the T-cell subsets ####
plot(sc_harmony@assays$RNA@counts["Trbv13-2",],
     sc_harmony@assays$RNA@counts["Trbv13-3",])

hist(sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",])
abline(v = quantile(sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",], 0.2), col="red")
# --- label TRBV13_2_3 positive ####
sc_harmony$TRBV13_2_pos <- ifelse(sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",] >= 0.03, "YES", "NO")
sc_harmony$TRBV13_3_pos <- ifelse(sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",] >= 0.04, "YES", "NO")
sc_harmony$TRBV13_2and3_pos <- ifelse(sc_harmony$TRBV13_3_pos == "YES" | sc_harmony$TRBV13_2_pos == "YES", "YES", "NO")


ggplot(sc_harmony@meta.data,
       aes(x=sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",],
           y=sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",],
           col=orig.ident)) +
  theme_bw() +
  geom_point(size=0.25) +
  geom_vline(xintercept = 0.04, linetype = "dashed", size = 0.75) +
  geom_hline(yintercept = 0.03, linetype = "dashed", size = 0.75) +
  scale_color_manual(values = c("CTR" = "gray", 
                                "Treated" = "orange")) +
  xlab("Trbv13-3") + ylab("Trbv13-2") +
  labs(col="Group") +
  facet_wrap(~ orig.ident) +
  stat_quadrant_counts(xintercept = 0.04,
                       yintercept = 0.03) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 1))


ggplot(sc_harmony@meta.data,
       aes(x=sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",],
           y=orig.ident,
           fill=orig.ident)) +
  geom_density_ridges(alpha=0.5, scale = 2.5) +
  theme_bw() +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = 0.03, linetype = "dashed", size = 0.75) +
  xlab("Magic-Imputed Expression") + ylab("") +
  scale_fill_manual(values = c("CTR" = "gray", 
                               "Treated" = "orange")) +
  labs(fill = "Group") + 
  ggtitle("Trbv13-2") + 
  theme(plot.title = element_text(vjust = -8, hjust = 0.9, size = 18)) +
  xlim(-0.02,0.25)

ggplot(sc_harmony@meta.data,
       aes(x=sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",],
           y=orig.ident,
           fill=orig.ident)) +
  geom_density_ridges(alpha=0.5, scale = 2.5) +
  theme_bw() +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = 0.04, linetype = "dashed", size = 0.75) +
  xlab("Magic-Imputed Expression") + ylab("") +
  scale_fill_manual(values = c("CTR" = "gray", 
                               "Treated" = "orange")) +
  labs(fill = "Group") +
  ggtitle("Trbv13-3") + 
  theme(plot.title = element_text(vjust = -8, hjust = 0.9, size = 18)) +
  xlim(-0.02,0.25)


as.data.frame.matrix(table(sc_harmony$TRBV13_2and3_pos,
                           sc_harmony$orig.ident)) %>%
  mutate(Perc_CTR = paste(round((CTR/sum(CTR))*100,2), "%"),
         Perc_Treated = paste(round((Treated/sum(Treated))*100,2), "%"))

as.data.frame.matrix(table(sc_harmony$TRBV13_2_pos,
                           sc_harmony$orig.ident)) %>%
  mutate(Perc_CTR = paste(round((CTR/sum(CTR))*100,2), "%"),
         Perc_Treated = paste(round((Treated/sum(Treated))*100,2), "%"))

as.data.frame.matrix(table(sc_harmony$TRBV13_3_pos,
                           sc_harmony$orig.ident)) %>%
  mutate(Perc_CTR = paste(round((CTR/sum(CTR))*100,2), "%"),
         Perc_Treated = paste(round((Treated/sum(Treated))*100,2), "%"))

# check cells that are double positive
TRBV13_2and3_pos_test <- ifelse(sc_harmony$TRBV13_3_pos == "YES" & sc_harmony$TRBV13_2_pos == "YES", "DOUBLE", sc_harmony$TRBV13_2and3_pos)
as.data.frame.matrix(table(TRBV13_2and3_pos_test,
                           sc_harmony$orig.ident)) %>%
  mutate(Perc_CTR = CTR/sum(CTR),
         Perc_Treated = Treated/sum(Treated))


table(sc_harmony$TRBV13_3_pos,
      sc_harmony$TRBV13_2_pos)
table(sc_harmony$new_labels,
      sc_harmony$orig.ident)

# prep tables
sc_harmony@meta.data[,c("TRBV13_2and3_pos", "orig.ident", "new_labels")] %>%
  count(TRBV13_2and3_pos, orig.ident, new_labels) %>%
  filter(TRBV13_2and3_pos == "YES") %>%
  ggplot(., aes(x=orig.ident,
                y=n,
                fill=orig.ident)) +
  geom_bar(stat="identity", position="dodge", width=0.75, col="black") +
  theme_classic() +
  scale_fill_manual(values = c("CTR" = "lightgray", "Treated" = "orange")) +
  ylab(paste("Number of cells positive for", "TRBV13-2 or TRBV13-3", sep="\n")) + xlab("") +
  facet_wrap(~ new_labels, scales = "free_y", ncol=6) +
  labs(fill = "Group")

table(sc_harmony$new_labels,
      sc_harmony$TRBV13_2and3_pos)
table(sc_harmony$new_labels,
      sc_harmony$TRBV13_2_pos)
table(sc_harmony$new_labels,
      sc_harmony$TRBV13_3_pos)

sc_harmony@meta.data[,c("TRBV13_2_pos", "orig.ident", "new_labels")] %>%
  count(TRBV13_2_pos, orig.ident, new_labels) %>%
  filter(TRBV13_2_pos == "YES")

table(paste0(sc_harmony$new_labels, "_",sc_harmony$TRBV13_2and3_pos),
      sc_harmony$orig.ident)



# run DGE for the greous that Jonathan shared
# DEG analysis for:
#  deg_1: Vehicle group TRBV13+ vs Vehicle group TRBV13- (there should be no/minimal difference at all)
#  deg_2: Vehicle group TRBV13- vs Treatment group TRBV13- (also anticipate minimal difference)
#  deg_3: Treatment group TRBV13- vs Treatment group TRBV13+ 
#  deg_4: Entire vehicle group vs Treatment group TRBV13+

cell_pop <- unique(sc_harmony$new_labels)
for(ii in 1:length(cell_pop)){
  
  message(paste("Working with", cell_pop[ii]))
  temp.cells <- subset(sc_harmony, cells=which(sc_harmony$new_labels == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "RNA"

  # dge_6
  temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
                                   ifelse(temp.cells$orig.ident == "CTR" & temp.cells$TRBV13_2and3_pos == "YES", "CTR","remove"))
  
  table(temp.cells$temp_groups)
  temp.cells <- subset(temp.cells, cells=which(temp.cells$temp_groups %in% c("CTR", "Treated")))
  
  Idents(temp.cells) <- as.character(temp.cells$temp_groups)
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.2, 
                                   min.pct = 0.2),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  # temp.res$FDR <- p.adjust(temp.res$p_val)
  
  if(ii == 1){
    res_deg_trbvs <- list(temp.res)
    names(res_deg_trbvs) <- cell_pop[ii]
  } else {
    res_deg_trbvs[[ii]] <- temp.res
    names(res_deg_trbvs)[ii] <- cell_pop[ii]
  }
  message(paste(" >> Done for", cell_pop[ii]))
}

res_deg_1 <- res_deg_trbvs
res_deg_1 <- rbindlist(res_deg_1, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_1, "all_res_deg_1.txt", sep="\t", quote = F, row.names = F)
table(res_deg_1$.id,
      res_deg_1$Direction)
res_deg_2 <- res_deg_trbvs
res_deg_2 <- rbindlist(res_deg_2, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_2, "all_res_deg_2.txt", sep="\t", quote = F, row.names = F)
table(res_deg_2$.id,
      res_deg_2$Direction)
res_deg_3 <- res_deg_trbvs
res_deg_3 <- rbindlist(res_deg_3, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_3, "all_res_deg_3.txt", sep="\t", quote = F, row.names = F)
table(res_deg_3$.id,
      res_deg_3$Direction)
res_deg_4 <- res_deg_trbvs
res_deg_4 <- rbindlist(res_deg_4, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_4, "all_res_deg_4.txt", sep="\t", quote = F, row.names = F)
table(res_deg_4$.id,
      res_deg_4$Direction)

res_deg_6 <- res_deg_trbvs
res_deg_6 <- rbindlist(res_deg_6, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_6, "all_res_deg_6.txt", sep="\t", quote = F, row.names = F)
table(res_deg_6$.id,
      res_deg_6$Direction)


heat <- as.data.frame.matrix(table(res_deg_6$.id,
                                   res_deg_6$Direction)) %>%
  mutate(NS = NULL)

Heatmap(heat,
        cluster_rows = T,
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,10,50,150),
                         colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
        border = T,
        name = "DGE_6")
#make UMAP splitting by both treatment and cell type
Idents(sc_harmony)
DimPlot(sc_harmony,
        group.by = "new_labels",
        cols = c25)+
  facet_wrap(sc_harmony$new_labels ~ sc_harmony$orig.ident, ncol=6)

# check other TRBVs
trbvs <- rownames(sc_harmony)[grep("Trbv", rownames(sc_harmony), ignore.case = T)]
trbv_cor <- cor(t(sc_harmony@assays$MAGIC_RNA@data[trbvs,]), method = "pearson")
Heatmap(trbv_cor,
        border = T,
        name = paste("Pearson's","Corr coeff", sep="\n"))

trbv_2_plot <- sc_harmony@assays$MAGIC_RNA@data[c("Trbv13-2","Trbv13-3","Trbv12-2","Trbv24"),]
pairs(t(trbv_2_plot), 
      cex=0.25, pch=16, col=scales::alpha("blue", 0.5),
      upper.panel = NULL)

Idents(sc_harmony) <- sc_harmony$orig.ident
table(Idents(sc_harmony))
all_t_vs_c <- tryCatch(FindMarkers(sc_harmony, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.2, 
                                   min.pct = 0.2),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       }) 
all_t_vs_c <- all_t_vs_c %>%
  mutate(Gene = rownames(.)) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))

table(all_t_vs_c$Direction)
write.table(all_t_vs_c, "all_DGE_Treated_vs_CTR.txt", sep="\t", quote = F, row.names = F)
ggplot(sc_harmony@meta.data,
       aes(x=sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",],
           y=orig.ident,
           fill=orig.ident)) +
  geom_density_ridges(alpha=0.5, scale = 2.5) +
  theme_bw() +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = 0.04, linetype = "dashed", size = 0.75) +
  xlab("Magic-Imputed Expression") + ylab("") +
  scale_fill_manual(values = c("CTR" = "gray", 
                               "Treated" = "orange")) +
  labs(fill = "Group") +
  ggtitle("Trbv13-3") + 
  theme(plot.title = element_text(vjust = -8, hjust = 0.9, size = 18)) +
  xlim(-0.02,0.25)

# make heatmap as requested by Jonathan via email 03.10.2022 ####
Idents(sc_harmony) <- sc_harmony$new_labels
table(Idents(sc_harmony))
markers_per_cell_type <- FindAllMarkers(sc_harmony, 
                                        only.pos = T,
                                        assay = "SCT", 
                                        logfc.threshold = .25, 
                                        min.pct = 0.5)
head(markers_per_cell_type)
length(grep("^Gm[0-9]", markers_per_cell_type$gene))
top <- markers_per_cell_type %>%
  arrange(desc(avg_log2FC)) %>%
  filter(!duplicated(gene)) %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  filter(!grepl("^Gm[0-9]", gene)) %>%
  filter(!grepl("^AY[0-9]", gene)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(10))
sum(duplicated(top$gene))

heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% top$gene),])
heat <- aggregate(heat, by=list(sc_harmony$new_labels), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
heat <- apply(heat,2,scale)
rownames(heat) <- rws
heat <- t(heat[match(unique(top$cluster), rownames(heat)),match(top$gene, colnames(heat))])
heat[1:5,1:5]
h_all <- Heatmap(heat, 
                 cluster_columns = F,
                 cluster_rows = F,
                 name="Expression")

markers_per_cell_type_ctr <- FindAllMarkers(subset(sc_harmony, cells = colnames(sc_harmony)[which(sc_harmony$orig.ident == "CTR")]),
                                            assay = "SCT", 
                                            only.pos = T,
                                            logfc.threshold = .25, 
                                            min.pct = 0.5)

top <- markers_per_cell_type_ctr %>%
  filter(cluster != "Monocytes") %>%
  arrange(desc(avg_log2FC)) %>%
  filter(!duplicated(gene)) %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  filter(!grepl("^Gm[0-9]", gene)) %>%
  filter(!grepl("^AY[0-9]", gene)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(10))
sum(duplicated(top$gene))
top_ctr <- top

markers_per_cell_type_treat <- FindAllMarkers(subset(sc_harmony, cells = colnames(sc_harmony)[which(sc_harmony$orig.ident == "Treated")]),
                                              assay = "SCT",
                                              only.pos = T, 
                                              logfc.threshold = .25, 
                                              min.pct = 0.5)

top <- markers_per_cell_type_treat %>%
  filter(cluster != "Monocytes") %>%
  arrange(desc(avg_log2FC)) %>%
  filter(!duplicated(gene)) %>%
  filter(avg_log2FC > 0.5 & p_val_adj < 0.05) %>%
  filter(!grepl("^Gm[0-9]", gene)) %>%
  filter(!grepl("^AY[0-9]", gene)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(10))
sum(duplicated(top$gene))
top_treated <- top

# heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(top_ctr$gene,top_treated$gene))),which(sc_harmony$orig.ident == "Treated"  & sc_harmony$new_labels != "Monocytes")])
heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(res_deg_6$Gene[which(res_deg_6$p_val_adj < 0.05 & abs(res_deg_6$avg_log2FC) > 0.75)]))),
                                           which(sc_harmony$orig.ident == "Treated"  & 
                                                   sc_harmony$new_labels %in% unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")] &
                                                   # sc_harmony$new_labels != "Monocytes" & 
                                                   sc_harmony$TRBV13_2and3_pos == "YES")])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "Treated"  & 
                                                              sc_harmony$new_labels %in% unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")] & 
                                                              # sc_harmony$new_labels != "Monocytes"  & 
                                                              sc_harmony$TRBV13_2and3_pos == "YES")]), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
# heat <- apply(heat,2,scale)
rownames(heat) <- rws
# heat_treated <- t(heat[match(unique(top$cluster), rownames(heat)),
heat_treated <- t(heat[match(unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")],rownames(heat)),
                       match(unique(c(res_deg_6$Gene[which(res_deg_6$p_val_adj < 0.05 & abs(res_deg_6$avg_log2FC) > 0.75)])), colnames(heat))])
heat_treated[1:5,]

heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(res_deg_6$Gene[which(res_deg_6$p_val_adj < 0.05 & abs(res_deg_6$avg_log2FC) > 0.75)]))),
                                           which(sc_harmony$orig.ident == "CTR" & 
                                                   # sc_harmony$new_labels != "Monocytes" & 
                                                   sc_harmony$new_labels %in% unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")] &
                                                   sc_harmony$TRBV13_2and3_pos == "YES")])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "CTR" & 
                                                              sc_harmony$new_labels %in% unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")] & 
                                                              # sc_harmony$new_labels != "Monocytes"  & 
                                                              sc_harmony$TRBV13_2and3_pos == "YES")]), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
# heat <- apply(heat,2,scale)
rownames(heat) <- rws
head(t(heat))
heat_ctr <- t(heat[match(unique(res_deg_6$`.id`)[which(unique(res_deg_6$`.id`) != "Monocytes" & unique(res_deg_6$`.id`) != "CD8_Effector_Exhauted")],rownames(heat)),
                   match(unique(c(res_deg_6$Gene[which(res_deg_6$p_val_adj < 0.05 & abs(res_deg_6$avg_log2FC) > 0.75)])), colnames(heat))])
heat_ctr[1:5,]
heat <- t(apply(cbind(heat_ctr, heat_treated),1,scale))
colnames(heat) <- rep(colnames(heat_ctr),2)

h_ctr <- Heatmap(heat[,1:9][,c(7,4,1,2,6,5,3,9,8)], 
                 cluster_columns = F,
                 # cluster_rows = F,
                 name=paste("Expression", "CTR","(left)", sep="\n"),
                 col = colorRamp2(breaks = c(-3,0,3),
                                  colors = c("#4169E1", "white", "red")),
                 border=T)
h_treated <- Heatmap(heat[,10:18][,match(colnames(h_ctr@matrix), colnames(heat)[10:18])], 
                     cluster_columns = F,
                     # cluster_rows = F,
                     name=paste("Expression", "Treated","(right)", sep="\n"),
                     row_names_gp = gpar(fontsize=8),
                     col = colorRamp2(breaks = c(-3,0,3),
                                      colors = c("#4169E1", "white", "red")),
                     border=T)

pdf("Heatmap_ctr_treated_by_cell_type_deg6_trbv_pos_only.pdf", height = 15, width = 8)
h_ctr + h_treated
dev.off()
gplots::venn(list(ctr = top_ctr$gene,
                  treated = top_treated$gene))

# run DEG_5 as requested via email ####
Idents(sc_harmony) <- "orig.ident"

cell_pop <- unique(sc_harmony$new_labels)
for(ii in 1:length(cell_pop)){
  
  message(paste("Working with", cell_pop[ii]))
  temp.cells <- subset(sc_harmony, cells=which(sc_harmony$new_labels == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "RNA"
  
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.2, 
                                   min.pct = 0.2),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  # temp.res$FDR <- p.adjust(temp.res$p_val)
  
  if(ii == 1){
    res_deg_overall <- list(temp.res)
    names(res_deg_overall) <- cell_pop[ii]
  } else {
    res_deg_overall[[ii]] <- temp.res
    names(res_deg_overall)[ii] <- cell_pop[ii]
  }
  res_deg_overall[[ii]]$CellType <- cell_pop[ii]
  res_deg_overall[[ii]]$Significant <- ifelse(res_deg_overall[[ii]]$p_val_adj < 0.05, "YES", "NO")
  res_deg_overall[[ii]]$Direction <- ifelse(res_deg_overall[[ii]]$Significant == "YES" & res_deg_overall[[ii]]$avg_log2FC > 0, "UP", 
                                            ifelse(res_deg_overall[[ii]]$Significant == "YES" & res_deg_overall[[ii]]$avg_log2FC < 0, "DOWN", "NS"))
  message(paste(" >> Done for", cell_pop[ii]))
}

res_deg_overall_df <- rbindlist(res_deg_overall, idcol = F)
head(res_deg_overall_df)
table(res_deg_overall_df$CellType,
      res_deg_overall_df$Direction)
write.table(res_deg_overall_df, "DEG_5.txt", sep = "\t", quote = F, row.names = F)

# heat map using the gene list from DEG5, with expression on TRBV13+ cells (vehicle vs treatment)
# DEG5 breakdown by cell type.
heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(res_deg_overall_df$Gene[which(res_deg_overall_df$p_val_adj < 0.05 & abs(res_deg_overall_df$avg_log2FC) > 1.25)]))),
                                           which(sc_harmony$orig.ident == "Treated"  & 
                                                   sc_harmony$TRBV13_2and3_pos == "YES" &
                                                   sc_harmony$new_labels %in% unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")])])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "Treated"  &  
                                                              sc_harmony$TRBV13_2and3_pos == "YES" &
                                                              sc_harmony$new_labels %in% unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")])]), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
rownames(heat) <- rws
heat_treated <- t(heat[match(unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")],rownames(heat)),
                       match(unique(c(res_deg_overall_df$Gene[which(res_deg_overall_df$p_val_adj < 0.05 & abs(res_deg_overall_df$avg_log2FC) > 1.25)])), colnames(heat))])
heat_treated[1:5,]

heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(res_deg_overall_df$Gene[which(res_deg_overall_df$p_val_adj < 0.05 & abs(res_deg_overall_df$avg_log2FC) > 1.25)]))),
                                           which(sc_harmony$orig.ident == "CTR"  &
                                                   sc_harmony$TRBV13_2and3_pos == "YES" & 
                                                   sc_harmony$new_labels %in% unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")])])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "CTR"  &  
                                                              sc_harmony$TRBV13_2and3_pos == "YES" &
                                                              sc_harmony$new_labels %in% unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")])]), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
rownames(heat) <- rws
heat_ctr <- t(heat[match(unique(res_deg_overall_df$CellType)[which(unique(res_deg_overall_df$CellType) != "Monocytes")],rownames(heat)),
                   match(unique(c(res_deg_overall_df$Gene[which(res_deg_overall_df$p_val_adj < 0.05 & abs(res_deg_overall_df$avg_log2FC) > 1.25)])), colnames(heat))])
heat_ctr[1:5,]


heat <- t(apply(cbind(heat_ctr, heat_treated),1,scale))
colnames(heat) <- rep(colnames(heat_treated),2)

h_ctr <- Heatmap(heat[,1:10][,c(7,1,8,2,6,5,4,3,9,10)], 
                 cluster_columns = F,
                 # cluster_rows = F,
                 name=paste("Expression", "CTR","(left)", sep="\n"),
                 col = colorRamp2(breaks = c(-3,0,3),
                                  colors = c("#4169E1", "white", "red")),
                 border=T)
h_treated <- Heatmap(heat[,11:20][,match(colnames(h_ctr@matrix), colnames(heat)[11:20])], 
                     cluster_columns = F,
                     # cluster_rows = F,
                     name=paste("Expression", "Treated","(right)", sep="\n"),
                     row_names_gp = gpar(fontsize=8),
                     col = colorRamp2(breaks = c(-3,0,3),
                                      colors = c("#4169E1", "white", "red")),
                     border=T)

pdf("Heatmap_deg_5_trbv.pdf", height = 17, width = 8)
h_ctr + h_treated
dev.off()

# as email from Jonathan ####
# run merge just TRBV5 and TRBV14? And make it green? Like the above?
DefaultAssay(sc_harmony) <- "MAGIC_RNA"
rownames(sc_harmony)[grep("trbv13", rownames(sc_harmony), ignore.case = T)]
sc_harmony$Trvb5_14 <- sc_harmony@assays$MAGIC_RNA@data["Trbv5",] + sc_harmony@assays$MAGIC_RNA@data["Trbv14",]
sc_harmony$Trvb132_133 <- sc_harmony@assays$MAGIC_RNA@data["Trbv13-2",] + sc_harmony@assays$MAGIC_RNA@data["Trbv13-3",]
table(sc_harmony$orig.ident)
FeaturePlot(sc_harmony, 
            features = c("Trvb5_14"), 
            # split.by = "orig.ident",
            cols = c("lightgray", "green"))

ggplot(NULL,
       aes(x = sc_harmony$Trvb5_14,
           y = sc_harmony$Trvb132_133)) +
  geom_point()


# sc_harmony$diff <- sc_harmony$Trvb132_133 - sc_harmony$Trvb5_14
sc_harmony$diff <- scales::rescale(sc_harmony$Trvb132_133, to = c(0,1)) - scales::rescale(sc_harmony$Trvb5_14, to=c(0,1))
hist(sc_harmony$diff)
ggplot(sc_harmony@meta.data,
       aes(x = sc_harmony@reductions$umap@cell.embeddings[,1],
           y = sc_harmony@reductions$umap@cell.embeddings[,2],
           col=diff)) +
  theme_classic() +
  geom_point(size=0.3, col=ifelse(sc_harmony$diff < 0.21 & sc_harmony$diff > -0.05, "gray",
                                  ifelse(sc_harmony$diff >= 0.21, "red",
                                         ifelse(sc_harmony$diff <= -0.05, "green", "white")))) +
  # xlim(-10,-7) + ylim(8, 11) +
  # geom_point(size=1.25) +
  # scale_color_gradientn(colors = c("green", scales::alpha("green", 0.75), "gray", "gray", scales::alpha("tomato", 0.5), "red"),
  #                       values = c(0, 0.125, 0.35, 0.45,  0.5, 1),
  #                       breaks = c(-0.015, -0.0125, -0.01, 0.1, 0.15, 0.125)) +
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(col="") +
  facet_wrap(~ orig.ident)

Idents(sc_harmony) <- sc_harmony$new_labels
DimPlot(sc_harmony, 
        cols = c25,
        label = T, 
        label.color = "white",
        label.box = T)

ggplot(sc_harmony@meta.data,
       aes(x = new_labels,
           y = diff,
           fill = orig.ident)) +
  geom_boxplot(alpha=0.1) +
  geom_violin(alpha = 0.5, scale="width") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, size=14)) +
  labs(fill = "Group")

# write.table(sc_harmony@meta.data, "metadata_040822_w_trbvDiff_no_newCells.txt", sep = "\t", col.names = NA)

# do DEG 3 as per meeting request ####
# To be done within only the treatment group, and looking at TRBV13+ vs TRBV13- 
# read in DEG 3 
setwd("~/Google Drive/BridgeBfx_SB/MarengoTx/scRNAseq_test_2_samples")
deg_3 <- read.table("all_res_deg_3.txt", header = T, sep="\t")
colnames(deg_3)[1] <- "CellType"
head(deg_3)
heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(deg_3$Gene[which(deg_3$p_val_adj < 0.05 & abs(deg_3$avg_log2FC) > 1)]))),
                                           which(sc_harmony$orig.ident == "Treated"  &
                                                   sc_harmony$TRBV13_2and3_pos == "YES" &
                                                   sc_harmony$new_labels %in% unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")])])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "Treated"  &  
                                                              sc_harmony$TRBV13_2and3_pos == "YES" &
                                                              sc_harmony$new_labels %in% unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")])]), 
                  FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
rownames(heat) <- rws
heat_treated <- t(heat[match(unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")],rownames(heat)),
                       match(unique(c(deg_3$Gene[which(deg_3$p_val_adj < 0.05 & abs(deg_3$avg_log2FC) > 1)])), colnames(heat))])
heat_treated[1:5,]

heat <- t(sc_harmony@assays$MAGIC_RNA@data[which(rownames(sc_harmony@assays$MAGIC_RNA@data) %in% unique(c(deg_3$Gene[which(deg_3$p_val_adj < 0.05 & abs(deg_3$avg_log2FC) > 1)]))),
                                           which(sc_harmony$orig.ident == "CTR"  &
                                                   sc_harmony$TRBV13_2and3_pos == "YES" & 
                                                   sc_harmony$new_labels %in% unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")])])
heat <- aggregate(heat, by=list(sc_harmony$new_labels[which(sc_harmony$orig.ident == "CTR"  &  
                                                              sc_harmony$TRBV13_2and3_pos == "YES" &
                                                              sc_harmony$new_labels %in% unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")])]), FUN=mean)
rws <- heat$Group.1
heat$Group.1 <- NULL
rownames(heat) <- rws
heat_ctr <- t(heat[match(unique(deg_3$CellType)[which(unique(deg_3$CellType) != "Monocytes")],rownames(heat)),
                   match(unique(c(deg_3$Gene[which(deg_3$p_val_adj < 0.05 & abs(deg_3$avg_log2FC) > 1)])), colnames(heat))])
heat_ctr[1:5,]


heat <- t(apply(cbind(heat_ctr, heat_treated),1,scale))
colnames(heat) <- rep(colnames(heat_treated),2)

h_ctr <- Heatmap(heat[,1:10][,c(2,1,6,7,4,5,3,8,9,10)], 
                 cluster_columns = F,
                 # cluster_rows = F,
                 name=paste("Expression", "CTR","(left)", sep="\n"),
                 col = colorRamp2(breaks = c(-3,0,3),
                                  colors = c("#4169E1", "white", "red")),
                 border=T)
h_treated <- Heatmap(heat[,11:20][,match(colnames(h_ctr@matrix), colnames(heat)[11:20])], 
                     cluster_columns = F,
                     # cluster_rows = F,  
                     name=paste("Expression", "Treated","(right)", sep="\n"),
                     row_names_gp = gpar(fontsize=8),
                     col = colorRamp2(breaks = c(-3,0,3),
                                      colors = c("#4169E1", "white", "red")),
                     border=T)

# pdf("Heatmap_deg_5_trbv.pdf", height = 17, width = 8)
h_ctr + h_treated
# dev.off()


# address questions in email 06.02.2022 ####
# 3)	When looking at the DGEs for each of these assigned groups, 
# some have many DGEs while others have nearly none. 
# I was wondering if this might be due to cell numbers in the cell types from each of the two groups?
# three values to visualize:
## change in cell number treated vs. control
## cell type
## number of DEGs
DimPlot(sc_harmony)
as.data.frame.matrix(table(sc_harmony$new_labels,
                           sc_harmony$orig.ident)) %>%
  mutate(DeltaCellNumber = Treated - CTR,
         CellType = rownames(.)) %>%
  melt(., id.vars = c("CellType", "DeltaCellNumber")) %>%
  ggplot(., aes(x = variable,
                y = value,
                fill = variable)) +
  geom_boxplot()

# get number of DGE from the saved data
deg_files <- list.files(".", pattern = "all_res_deg_")
for(ii in 1:length(deg_files)){
  temp_res <- read.table(deg_files[ii], header = T, sep='\t') %>%
    filter(Direction != "NS")
  if(ii == 5){
    temp_res <- as.data.frame.matrix(table(temp_res$CellType, temp_res$Direction)) %>%
      mutate(CellType = rownames(.))
  } else {
    colnames(temp_res)[1] <- "CellType"
    temp_res <- as.data.frame.matrix(table(temp_res$CellType, temp_res$Direction)) %>%
      mutate(CellType = rownames(.))
  }
  colnames(temp_res)[1:2] <- c("DOWN", "UP")
  temp_res$DEG_Group <- gsub("all_res_","",gsub("\\.txt","",deg_files[ii]))
  # temp_res <- as.data.frame.matrix(table(temp_res$CellType, temp_res$Direction)) %>%
  #   mutate(CellType = rownames(.))
  
  
  if(ii == 1){
    table_all <- temp_res
  } else {
    table_all <- data.frame(rbind(table_all, 
                                  temp_res))
  }
}
head(temp_res)
tt_ncells <- as.data.frame.matrix(table(sc_harmony$new_labels,
                                        sc_harmony$orig.ident))
table_all <- merge(table_all, tt_ncells, by.x = "CellType", by.y="row.names") %>%
  mutate(DeltaNCells = Treated - CTR,
         IO_Effect = ifelse(CellType %in% c("CD8_Effector", "CD8_Effector_Memory", "CD4_T_helper"), "Pro",
                            ifelse(CellType %in% c("CD4_CM", "CD8_Naive_CM"), "Neutral","Suppressive")))
head(table_all)
ggplot(table_all, 
       aes(x = factor(IO_Effect, levels = c("Suppressive", "Neutral", "Pro")),
           y = DeltaNCells)) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(size = (UP + DOWN),
                  col = CellType,
                  shape = DEG_Group), 
              alpha=0.75) +
  labs(size = "(nDEG)") + xlab("IO Effect") +
  scale_color_manual(values = c25) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2),
         shape = guide_legend(override.aes = list(size=5), ncol = 2)) +
  theme(axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = list(c("Suppressive", "Neutral"),
                                                           c("Suppressive", "Pro"),
                                                           c("Neutral", "Pro")))

# it might be good to make some heatmaps for selected genes. 
# I’d like your input on how we could best do this. Attached is a gene list, 
# though we might for visual effect need to split it in two. 
# Suggestion would be to plot these genes across all assigned cell types?
genes2use <- c(as.character(data.frame(readxl::read_xlsx("../human_scRNASeq_VDJ_CSP/S565/Additional genes for plotting 06062022.xlsx"))[,1]), "FASL", "LY6C1", "LY6C2")
genes2use[-which(genes2use %in% toupper(rownames(sc_harmony)))]
row_ind <- which(toupper(rownames(sc_harmony@assays$RNA@data)) %in% genes2use)
df_meta <- data.frame(t(sc_harmony@assays$RNA@data[row_ind, ]),
                      Group = sc_harmony@meta.data[,c("orig.ident", "new_labels", "diff")])
head(df_meta)
heat <- t(apply(df_meta[,1:59], 1, scale))
colnames(heat) <- colnames(df_meta)[1:59]
df_heat <- df_meta[,60:62]
colnames(df_heat) <- c("Treatment", "CellType", "TRBV_Diff")
Heatmap(t(heat),
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(CellType=c("CD8_Precursor_Exhausted" = c25[1],"CD4_Treg" = c25[2],"CD8_Effector_Memory" = c25[3],
                                                                 "CD4_T_helper" = c25[4],"CD8_Memory_Exhausted" = c25[5],"Monocytes" = c25[6],
                                                                 "CD8_Effector" = c25[7],"CD4_CM" = c25[8],"CD8_Naive_CM" = c25[9],"CD8_Exhausted" = c25[10],
                                                                 "CD8_Effector_Exhauted" = c25[11]),
                                                      Treatment = c("CTR" = "lightgray", "Treated" = "orange"),
                                                      TRBV_Diff = colorRamp2(breaks = c(-1,0,1), colors = c("#DDA0DD", "#FFFAFA", "#D2691E")))),
        col = colorRamp2(breaks = c(-3,0,3), colors = c("#66CDAA", "#E0FFFF", "tomato")),
        show_row_names = T,
        show_column_names = F,
        use_raster = F,
        name = "Scaled Expression")


# make heatmap as requested in email 06082022 ####
table(sc_harmony$orig.ident)
heat_ctr <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$orig.ident == "CTR")), 
                              assays = "MAGIC_RNA",
                              features = str_to_sentence(genes2use),
                              return.seurat = F, 
                              group.by = 'new_labels')[[1]]
cls_ctr <- colnames(heat_ctr)

heat_treat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$orig.ident == "Treated")), 
                                assays = "MAGIC_RNA",
                                features = str_to_sentence(genes2use),
                                return.seurat = F, 
                                group.by = 'new_labels')[[1]]
cls_treat <- colnames(heat_treat)
heat <- t(apply(cbind(heat_ctr, heat_treat), 1, scale))
colnames(heat) <- c(cls_ctr, cls_treat)

heat_c <- Heatmap(heat[,1:11], 
                  cluster_columns = F,
                  # cluster_rows = F,
                  name=paste("Expression", "CTR","(left)", sep="\n"),
                  col = colorRamp2(breaks = c(-3,0,3),
                                   colors = c("#4169E1", "white", "red")),
                  border=T)
heat_t <- Heatmap(heat[,12:22][,match(colnames(heat_c@matrix), colnames(heat)[12:22])], 
                  cluster_columns = F,
                  name=paste("Expression", "Treated","(right)", sep="\n"),
                  row_names_gp = gpar(fontsize=14),
                  col = colorRamp2(breaks = c(-3,0,3),
                                   colors = c("#4169E1", "white", "red")),
                  border=T)

pdf("Heatmap_deg_5_trbv.pdf", height = 17, width = 8)
heat_c + heat_t
dev.off()
VlnPlot(sc_harmony, features = c("Il2ra", "Icos"), group.by = "new_labels", split.by = "orig.ident", assay = "MAGIC_RNA")
res_all %>% filter(Gene %in% c("Il2ra", "Icos", "Tcf7", "Gzmb") & Significant == "YES")
# address email 06172022 and notes from meeting ####
head(table_all)
ggplot(table_all %>%
         mutate(TotalDGE = UP+DOWN,
                log2Ncells = log2(Treated/CTR)), 
       aes(x = factor(IO_Effect, levels = c("Suppressive", "Neutral", "Pro")),
           y = TotalDGE)) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(size = log2Ncells,
                  col = CellType,
                  shape = DEG_Group), 
              alpha=0.75) +
  labs(size = "(log2Ncells)") + xlab("IO Effect") +
  scale_color_manual(values = c25) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2),
         shape = guide_legend(override.aes = list(size=5), ncol = 2)) +
  theme(axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = list(c("Suppressive", "Neutral"),
                                                           c("Suppressive", "Pro"),
                                                           c("Neutral", "Pro")))


ggplot(table_all %>%
         mutate(TotalDGE = UP+DOWN,
                log2Ncells = log2(Treated/CTR)), 
       aes(x = factor(IO_Effect, levels = c("Suppressive", "Neutral", "Pro")),
           y = log2Ncells)) +
  geom_boxplot() +
  theme_classic() +
  geom_jitter(aes(size = TotalDGE,
                  col = CellType,
                  shape = DEG_Group), 
              alpha=0.75) +
  labs(size = "(log2Ncells)") + xlab("IO Effect") +
  scale_color_manual(values = c25) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 2),
         shape = guide_legend(override.aes = list(size=5), ncol = 2)) +
  theme(axis.text = element_text(size = 12)) +
  stat_compare_means(method = "t.test", comparisons = list(c("Suppressive", "Neutral"),
                                                           c("Suppressive", "Pro"),
                                                           c("Neutral", "Pro")))


# -	Plot x = delta cells, y = percentage change, col = group
ggplot(table_all %>%
         mutate(TotalDGE = UP+DOWN,
                log2Ncells = log2(Treated/CTR)), 
       aes(x = DeltaNCells,
           y = log2Ncells,
           col = IO_Effect,
           size = TotalDGE)) +
  geom_point() +
  theme_classic() +
  # geom_jitter(aes(size = TotalDGE,
  #                 col = CellType,
  #                 shape = DEG_Group), 
  #             alpha=0.75) +
  labs(size = "(log2Ncells)") + xlab("DeltaNCells (Treat - CTR)") +
  scale_color_manual(values = c25) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 1)) +
  theme(axis.text = element_text(size = 12))

# Jonathan email 06172022 ####
# Comparing treated vs vehicle for each of the T-cell subset/labels
# pdf("vln_genes_byCellType_06172022.pdf", height = 15, width = 7)
sc_harmony$new_labels_treatment <- paste0(sc_harmony$new_labels, "-", sc_harmony$orig.ident)
heat <- AverageExpression(sc_harmony, 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "red",
                                                                    "Treated" = "green"),
                                                      CellType = c("CD4_CM"  = c25[1],
                                                                   "CD4_T_helper"  = c25[2],
                                                                   "CD4_Treg"  = c25[3],
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = c25[5],
                                                                   "CD8_Effector" = c25[6],
                                                                   "CD8_Exhausted"  = c25[7],
                                                                   "CD8_Memory_Exhausted"  = c25[8],
                                                                   "CD8_Naive_CM"  = c25[9],
                                                                   "CD8_Precursor_Exhausted"  = c25[10],
                                                                   "Monocytes" = c25[11]))),
        name = "Expression")

VlnPlot(sc_harmony,
        assay = "MAGIC_RNA",
        features = str_to_sentence(genes2use),
        stack = T, flip = T,
        pt.size = 0,
        group.by = "new_labels",
        split.by = "orig.ident")
# dev.off()

# address email 06212022 ####
# a) repeat the same heatmap but only keeping 
# CD8_effector, CD8_effector_memory, CD8_naive_CM, CD8_memory_exhausted_CD4_treg. 
# And removing all other cell types
sc_harmony$new_labels_treatment <- paste0(sc_harmony$new_labels, "-", sc_harmony$orig.ident)
table(sc_harmony$new_labels_treatment)
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                        "CD8_Effector_Memory",
                                                                                        "CD8_Naive_CM",
                                                                                        "CD8_Memory_Exhausted",
                                                                                        "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")

# b) Is it possible to generate the same heatmap gene expression 
# comparing vehicle samples, vs just the TRBV13+ population in treated samples?, i.e. DEG4.  
# And repeat it for as in (a)
sc_harmony$new_labels_treatment <- paste0(sc_harmony$new_labels, "-", sc_harmony$orig.ident)
table(sc_harmony$new_labels_treatment)
heat <- AverageExpression(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                                 "CD8_Effector_Memory",
                                                                                                 "CD8_Naive_CM",
                                                                                                 "CD8_Memory_Exhausted",
                                                                                                 "CD4_Treg") & sc_harmony$orig.ident == "CTR"),
                                                              which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                                 "CD8_Effector_Memory",
                                                                                                 "CD8_Naive_CM",
                                                                                                 "CD8_Memory_Exhausted",
                                                                                                 "CD4_Treg") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES"))),
                                 
), 
assays = "MAGIC_RNA",
features = str_to_sentence(genes2use), 
group.by = "new_labels_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")


# same as above but with DEG4
deg_4 <- read.table("all_res_deg_4.txt", header = T, sep="\t")
colnames(deg_4)[1] <- "CellType"
deg_4 <- deg_4[which(deg_4$CellType %in% df_heat$CellType),]
deg_4 <- deg_4[which(deg_4$Significant == "YES"),]
head(deg_4)

table(sc_harmony$new_labels_treatment)
heat <- AverageExpression(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                                 "CD8_Effector_Memory",
                                                                                                 "CD8_Naive_CM",
                                                                                                 "CD8_Memory_Exhausted",
                                                                                                 "CD4_Treg") & sc_harmony$orig.ident == "CTR"),
                                                              which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                                 "CD8_Effector_Memory",
                                                                                                 "CD8_Naive_CM",
                                                                                                 "CD8_Memory_Exhausted",
                                                                                                 "CD4_Treg") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES"))),
), 
assays = "MAGIC_RNA",
features = unique(deg_4$Gene), 
group.by = "new_labels_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
pdf("heat_deg4_large.pdf", width = 12, height = 120)
Heatmap(heat,
        show_column_names = F,
        show_row_names = T,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")
dev.off()

# c) d) So back to UMAP plots as we discussed regarding to colors, some of the key populations, 
# can you use for the following (see above)
Idents(sc_harmony) <- sc_harmony$new_labels
DimPlot(sc_harmony, 
        # label = T, 
        # label.color = "white", 
        # repel = T,
        # label.box = T,
        cols = c("CD4_Treg"  = "green",
                 "CD8_Effector_Exhauted"  = c25[4],
                 "CD8_Effector_Memory"  = "#008080",
                 "CD8_Effector" = "black",
                 "CD8_Memory_Exhausted"  = "red",
                 "CD8_Naive_CM"  = "orange",
                 "CD4_CM"  = c25[1],
                 "CD4_T_helper"  = "purple",
                 "CD8_Exhausted"  = c25[7],
                 "CD8_Precursor_Exhausted"  = c25[10],
                 "Monocytes" = c25[11]),
        split.by = "orig.ident")

# And finally, to zoom in on the heatmap (slide11), 
# can you generate violin plots for the following genes: 
# (can you use ctrl: color gray, and treated: color orange) 
genes2chck <- c("Ly6c1", "Ly6c2", "Il2ra", "Bcl2", "Gzma", "Gzmb", "Gzmk", "Plac8", "Fos", "Jun", "Junb", "Tcf7", "Klf2", "Tox", "Runx3", "Fasl", "Ccl5", "Cd69", "Nkg7", "Ifng", "Klrg1", "Tnfrsf4", "Tox", "Prdm1", "Tnfrsf9", "Ctla4", "Eomes", "Lag3", "Havcr2", "Pdcd1", "Cd28")
genes2chck[-which(genes2chck %in% rownames(sc_harmony))]

pdf("Violins_all_cell_types.pdf", width = 5, height = 4)
sapply(genes2chck, function(xx){
  vv <- VlnPlot(sc_harmony,
                assay = "MAGIC_RNA",
                features = xx,
                group.by = "new_labels",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange"))
  print(vv)
})
dev.off()

pdf("Violins_TRBV_veh_selected_cell_types.pdf", width = 5, height = 4)
sapply(genes2chck, function(xx){
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                       "CD8_Effector_Memory",
                                                                                       "CD8_Naive_CM",
                                                                                       "CD8_Memory_Exhausted",
                                                                                       "CD4_Treg") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                       "CD8_Effector_Memory",
                                                                                       "CD8_Naive_CM",
                                                                                       "CD8_Memory_Exhausted",
                                                                                       "CD4_Treg") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = xx,
                group.by = "new_labels",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange"))
  print(vv)
})
dev.off()

# CD8_effector memory cell type subsets and CD8_effector subsets (since these are the most important ones), 
# and can you use DEG4, i.e. gene expression comparison in those CD8s subsets, 
# between TRBV13+ cells (treated samples) vs all cells (vehicle)
deg_4 <- deg_4[which(deg_4$CellType %in% c("CD8_Effector", "CD8_Effector_Memory")),]
deg_4 <- deg_4[which(abs(deg_4$avg_log2FC) >= 0.5),]
pdf("Violins_DEG4_logFC1_TRBV_veh_selected_cell_types.pdf", width = 5, height = 4)
sapply(unique(deg_4$Gene), function(xx){
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                       "CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels %in% c("CD8_Effector",
                                                                                       "CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = xx,
                group.by = "new_labels",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange"))
  print(vv)
})
dev.off()

# violin plots, can you make them individually, not bundled together
pdf("Violins_all_cell_types_not_stacked_v2.pdf", width = 5, height = 4)
sapply(str_to_sentence(genes2use)[-which(str_to_sentence(genes2use) %in% c("Faslg", "Ly6h"))], function(xx){
  vv <- VlnPlot(sc_harmony,
                assay = "MAGIC_RNA",
                features = xx,
                pt.size = 0,
                group.by = "new_labels",
                split.by = "orig.ident",
                cols = c("gray", "orange"))
  print(vv)
})
dev.off()


# send the violin plots for just the CD8_effector_memory subsets ####
deg_4 <- deg_4[which(deg_4$CellType == "CD8_Effector_Memory"),]
pdf("Violins_CD8_effector_memory_DEG4.pdf", width = 3, height = 4)
sapply(str_to_sentence(genes2use)[-which(str_to_sentence(genes2use) %in% c("Faslg", "Ly6h"))], function(xx){
  vv <- VlnPlot(subset(sc_harmony, cells = which(sc_harmony$new_labels == "CD8_Effector_Memory")),
                assay = "MAGIC_RNA",
                features = xx,
                pt.size = 0,
                group.by = "new_labels",
                split.by = "orig.ident",
                cols = c("gray", "orange"))
  print(vv)
})
dev.off()


DefaultAssay(sc_harmony) <- "MAGIC_RNA"
FeatureScatter(sc_harmony,
               feature1 = "Cd8a",
               feature2 = "Cd4",
               group.by = "new_labels") +
  facet_wrap(~ sc_harmony$new_labels)
FeatureScatter(sc_harmony,
               feature1 = "Cd8a",
               feature2 = "Cd4",
               pt.size = 0.25) +
  geom_vline(xintercept = 0.4, size=0.5, linetype = "dashed") +
  geom_hline(yintercept = 0.4, size=0.5, linetype = "dashed")

sc_harmony$new_labels_edited <- ifelse(sc_harmony@assays$MAGIC_RNA@data["Cd4",] < 0.4 &
                                         sc_harmony@assays$MAGIC_RNA@data["Cd8a",] < 0.4, "Monocytes",
                                       ifelse(sc_harmony$new_labels == "CD8_Effector" & 
                                                sc_harmony@assays$MAGIC_RNA@data["Cd4",] >= 0.4 &
                                                sc_harmony@assays$MAGIC_RNA@data["Cd8a",] < 0.4, "CD4_T_helper",
                                              ifelse(sc_harmony$new_labels == "CD8_Effector_Memory" & 
                                                       sc_harmony@assays$MAGIC_RNA@data["Cd4",] >= 0.4 &
                                                       sc_harmony@assays$MAGIC_RNA@data["Cd8a",] < 0.4, "CD4_Effector_Memory",
                                                     ifelse(sc_harmony$new_labels == "CD8_Precursor_Exhausted" & 
                                                              sc_harmony@assays$MAGIC_RNA@data["Cd4",] >= 0.4 &
                                                              sc_harmony@assays$MAGIC_RNA@data["Cd8a",] < 0.4, "CD4_Precursor_Exhausted",
                                                            sc_harmony$new_labels))))
as.data.frame.matrix(table(sc_harmony$new_labels_edited,
                           sc_harmony$new_labels)) %>% +1 %>%
  log10() %>%
  Heatmap(.,
          name = "Log10(nCells + 1)",
          border = T, 
          row_title = "Edited annotation",
          row_title_side = "right",
          column_title = "Original annotation",
          column_title_side = "bottom",
          col = colorRamp2(breaks = c(0,2,4),
                           colors = c("white", "lightgray", "black")))

FeatureScatter(sc_harmony,
               feature1 = "Cd8a",
               feature2 = "Cd4",
               group.by = "new_labels_edited") +
  facet_wrap(~ sc_harmony$new_labels_edited)


# get DEG 4 for violins ####
deg_4 <- read.table("all_res_deg_4.txt", header = T, sep="\t")
colnames(deg_4)[1] <- "CellType"
deg_4 <- deg_4[which(deg_4$CellType == "CD8_Effector_Memory"),]
deg_4 <- deg_4[which(deg_4$Significant == "YES"),]
deg_4 <- deg_4[which(abs(deg_4$avg_log2FC) >= 0.5),]

pdf("Violins_CD8_effector_memory_DEG4_edited_labels.pdf", width = 3.4, height = 4)
sapply(unique(deg_4$Gene), function(xx){
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = xx,
                # group.by = "new_labels_edited",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange")) +
    theme(axis.text.x = element_blank())
    # stat_compare_means(method = "t.test", vjust = 0.75, step.increase = 0.5)
  print(vv)
})
dev.off()

pdf("Violins_CD8_effector_memory_DEG4_original_labels.pdf", width = 3.4, height = 4)
sapply(unique(deg_4$Gene), function(xx){
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = xx,
                # group.by = "new_labels_edited",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange")) +
    theme(axis.text.x = element_blank())
  # stat_compare_means(method = "t.test", vjust = 0.75, step.increase = 0.5)
  print(vv)
})
dev.off()

for(ii in 1:length(unique(deg_4$Gene))){
  temp_gene <- unique(deg_4$Gene)[ii]
  if(ii == 1)
    pdf("Violins_CD8_effector_memory_DEG4_edited_labels_pv.pdf", width = 3.4, height = 3)
  
  title2use <- paste(temp_gene, 
                     paste0("adj.pv = ",
                            strsplit(as.character(deg_4$p_val_adj[which(deg_4$CellType == "CD8_Effector_Memory" & deg_4$Gene == temp_gene)]), "\\.")[[1]][1],
                            ".",
                            substr(strsplit(as.character(deg_4$p_val_adj[which(deg_4$CellType == "CD8_Effector_Memory" & deg_4$Gene == temp_gene)]), "\\.")[[1]][2], 1, 3),
                            "e-",
                            strsplit(as.character(deg_4$p_val_adj[which(deg_4$CellType == "CD8_Effector_Memory" & deg_4$Gene == temp_gene)]), "-")[[1]][2]), 
                     sep="\n")
  if(ii == 42){
    title2use <- paste(temp_gene, 
                       paste0("adj.pv = 0"), 
                       sep="\n")
  }
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = temp_gene,
                # group.by = "new_labels_edited",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange")) +
    theme(axis.text.x = element_blank()) +
    ggtitle(title2use)
  # stat_compare_means(method = "t.test", vjust = 0.75, step.increase = 0.5)
  print(vv)
  message(ii)
}
dev.off()

VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                      "CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                            which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                      "CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
        assay = "MAGIC_RNA",
        features = c("Cd8a", "Cd4"),
        group.by = "new_labels_edited",
        pt.size = 0,
        split.by = "orig.ident",
        cols = c("gray", "orange"))

# address Jonathan email 06.29.2022 ####
# Red more closer to orangy/red? And green more to purple.
ggplot(sc_harmony@meta.data,
       aes(x = sc_harmony@reductions$umap@cell.embeddings[,1],
           y = sc_harmony@reductions$umap@cell.embeddings[,2],
           col=diff)) +
  theme_classic() +
  geom_point(size=0.3, col=ifelse(sc_harmony$diff < 0.21 & sc_harmony$diff > -0.05, scales::alpha("gray", 0.75),
                                  ifelse(sc_harmony$diff >= 0.21, scales::alpha("tomato", 0.75),
                                         ifelse(sc_harmony$diff <= -0.05, scales::alpha("#4C0099", 0.5), "white")))) +
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(col="") +
  facet_wrap(~ orig.ident)


# address Jonathan email from 06.30.2022 ####
head(sc_harmony)
table(sc_harmony$new_labels_edited)
Idents(sc_harmony) <- sc_harmony$new_labels_edited
DimPlot(sc_harmony, 
        # label = T,
        # label.color = "white",
        # repel = T,
        # label.box = T,
        cols = c("CD4_Treg"  = "green",
                 "CD8_Effector_Exhauted"  = c25[4],
                 "CD8_Effector_Memory"  = "#008080",
                 "CD8_Effector" = "black",
                 "CD8_Memory_Exhausted"  = "red",
                 "CD8_Naive_CM"  = "orange",
                 "CD4_CM"  = c25[1],
                 "CD4_T_helper"  = "purple",
                 "CD8_Exhausted"  = c25[7],
                 "CD8_Precursor_Exhausted"  = c25[10],
                 "Monocytes" = c25[11],
                 "CD4_Effector_Memory" = c25[15],
                 "CD4_Precursor_Exhausted" = c25[13]),
        split.by = "orig.ident")

# prep barplots
sc_harmony@meta.data[,c("new_labels_edited", "orig.ident")] %>%
  count(new_labels_edited, orig.ident) %>%
  ggplot(., aes(x=new_labels_edited,
                y=n,
                fill=factor(orig.ident, levels = c("Treated", "CTR")))) +
  geom_bar(stat="identity", position="dodge", width=0.75, col="black") +
  theme_classic() +
  scale_fill_manual(values = c("CTR" = "lightgray", "Treated" = "orange")) +
  ylab(paste("Number of cells")) + xlab("") +
  coord_flip() +
  # facet_wrap(~ new_labels_edited, scales = "free_y", ncol=6) +
  labs(fill = "Group")

# run new DGE analysis with updated labels ####
table(sc_harmony$new_labels_edited, sc_harmony$orig.ident, sc_harmony$TRBV13_2and3_pos)
# deg1 - Vehicle group TRBV13+ vs Vehicle group TRBV13- (there should be no difference at all)
# deg2 - Vehicle group TRBV13- vs Treatment group TRBV13- (also anticipate minimal difference)
# deg3 - Treatment group TRBV13- vs Treatment group TRBV13+ 
# deg4 - Entire vehicle group vs Treatment group TRBV13+
# deg5 - Entire vehicle vs. entire treated
# deg6 - Vehicle TRBV13+ vs Treated TRBV13+
cell_pop <- unique(sc_harmony$new_labels_edited)
for(ii in 1:length(cell_pop)){
  
  message(paste("Working with", cell_pop[ii]))
  temp.cells <- subset(sc_harmony, cells=which(sc_harmony$new_labels_edited == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "RNA"
  
  # dge_1
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "CTR" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
  #                                  ifelse(temp.cells$orig.ident == "CTR" & temp.cells$TRBV13_2and3_pos == "NO", "CTR", "remove"))
  
  
  # dge_2
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "NO", "Treated",
  #                                  ifelse(temp.cells$orig.ident == "CTR" & temp.cells$TRBV13_2and3_pos == "NO", "CTR", "remove"))
  # 
  # dge_3
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
  #                                  ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "NO", "CTR", "remove"))
  
  # dge_4
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
  #                                   ifelse(temp.cells$orig.ident == "CTR", "CTR", "remove"))

  # dge_5
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated", "Treated",
  #                                   ifelse(temp.cells$orig.ident == "CTR", "CTR", "remove"))
  
  # dge_6
  temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
                                   ifelse(temp.cells$orig.ident == "CTR" & temp.cells$TRBV13_2and3_pos == "YES", "CTR","remove"))
  
  table(temp.cells$temp_groups)
  temp.cells <- subset(temp.cells, cells=which(temp.cells$temp_groups %in% c("CTR", "Treated")))
  
  Idents(temp.cells) <- as.character(temp.cells$temp_groups)
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   # slot = "data",
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.1, 
                                   min.pct = 0.1),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  # temp.res$FDR <- p.adjust(temp.res$p_val)
  
  if(ii == 1){
    res_deg_trbvs <- list(temp.res)
    names(res_deg_trbvs) <- cell_pop[ii]
  } else {
    res_deg_trbvs[[ii]] <- temp.res
    names(res_deg_trbvs)[ii] <- cell_pop[ii]
  }
  message(paste(" >> Done for", cell_pop[ii]))
}

res_deg_1 <- res_deg_trbvs
res_deg_1 <- rbindlist(res_deg_1, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_1, "all_res_deg_1_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_1$.id,
      res_deg_1$Direction)
as.data.frame.matrix(table(res_deg_1$.id,
                           res_deg_1$Direction)) %>%
  mutate(NS = NULL) %>%
  subset(rownames(.) != "Myeloid Lineage") %>%
Heatmap(.,
        cluster_rows = T,
        cluster_columns = F,
        col = colorRamp2(breaks = c(0,10,50,150),
                         colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
        border = T, 
        name = "DGE_1")


res_deg_2 <- res_deg_trbvs
res_deg_2 <- rbindlist(res_deg_2, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_2, "all_res_deg_2_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_2$.id,
      res_deg_2$Direction)
as.data.frame.matrix(table(res_deg_2$.id,
                           res_deg_2$Direction)) %>%
  mutate(NS = NULL) %>%
  subset(rownames(.) != "Myeloid Lineage") %>%
  Heatmap(.,
          cluster_rows = T,
          cluster_columns = F,
          col = colorRamp2(breaks = c(0,10,50,150),
                           colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
          border = T, 
          name = "DGE_2")


res_deg_3 <- res_deg_trbvs
res_deg_3 <- rbindlist(res_deg_3, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_3, "all_res_deg_3_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_3$.id,
      res_deg_3$Direction)
as.data.frame.matrix(table(res_deg_3$.id,
                           res_deg_3$Direction)) %>%
  mutate(NS = NULL) %>%
  subset(rownames(.) != "Myeloid Lineage") %>%
  Heatmap(.,
          cluster_rows = T,
          cluster_columns = F,
          col = colorRamp2(breaks = c(0,10,50,150),
                           colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
          border = T, 
          name = "DGE_3")

res_deg_4 <- res_deg_trbvs
res_deg_4 <- rbindlist(res_deg_4, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_4, "all_res_deg_4_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_4$.id,
      res_deg_4$Direction)
as.data.frame.matrix(table(res_deg_4$.id,
                           res_deg_4$Direction)) %>%
  mutate(NS = NULL) %>%
  filter(rownames(.) != "Myeloid Lineage") %>%
  Heatmap(.,
          cluster_rows = T,
          cluster_columns = F,
          col = colorRamp2(breaks = c(0,10,50,150),
                           colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
          border = T, 
          name = "DGE_4")


res_deg_5 <- res_deg_trbvs
res_deg_5 <- rbindlist(res_deg_5, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_5, "all_res_deg_5_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_5$.id,
      res_deg_5$Direction)
as.data.frame.matrix(table(res_deg_5$.id,
                           res_deg_5$Direction)) %>%
  mutate(NS = NULL) %>%
  filter(rownames(.) != "Myeloid Lineage") %>%
  Heatmap(.,
          cluster_rows = T,
          cluster_columns = F,
          col = colorRamp2(breaks = c(0,10,50,150),
                           colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
          border = T, 
          name = "DGE_5")

res_deg_6 <- res_deg_trbvs
res_deg_6 <- rbindlist(res_deg_6, use.names = T, idcol = T) %>%
  mutate(Significant = ifelse(p_val_adj <= 0.05, "YES", "NO")) %>%
  mutate(Direction = ifelse(Significant == "YES" & avg_log2FC > 0, "UP",
                            ifelse(Significant == "YES" & avg_log2FC < 0, "DWN", "NS")))
write.table(res_deg_6, "all_res_deg_6_07172022.txt", sep="\t", quote = F, row.names = F)
table(res_deg_6$.id,
      res_deg_6$Direction)
as.data.frame.matrix(table(res_deg_6$.id,
                           res_deg_6$Direction)) %>%
  mutate(NS = NULL) %>%
  filter(rownames(.) != "Myeloid Lineage") %>%
  Heatmap(.,
          cluster_rows = T,
          cluster_columns = F,
          col = colorRamp2(breaks = c(0,10,50,150),
                           colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
          border = T, 
          name = "DGE_6")

# do tables for all DEG as one
res_tables_list <- list(res_deg_1, res_deg_2, res_deg_3, res_deg_4, res_deg_5, res_deg_6)
names(res_tables_list) <- c("deg_1", "deg_2", "deg_3", "deg_4", "deg_5", "deg_6")
for(ii in 1:length(res_tables_list)){
  temp_res <- data.frame(res_tables_list[[ii]])
  temp_tt <- as.data.frame.matrix(table(temp_res$.id,
                                        temp_res$Direction))
  colnames(temp_tt) <- paste0(colnames(temp_tt), "_", names(res_tables_list)[ii])
  if(ii == 1){
    all_res_tables <- temp_tt
  } else {
    all_res_tables <- data.frame(cbind(all_res_tables, temp_tt))
  }
}
all_res_tables[,-grep("NS", colnames(all_res_tables))]

# do violins
genes2use <- c(as.character(data.frame(readxl::read_xlsx("../human_scRNASeq_VDJ_CSP/S565/Additional genes for plotting 06062022.xlsx"))[,1]), "FASL", "LY6C1", "LY6C2")
genes2use[-which(str_to_sentence(genes2use) %in% rownames(sc_harmony))]
row_ind <- which(toupper(rownames(sc_harmony@assays$RNA@data)) %in% genes2use)

pdf("Violins_CD8_effector_memory_DEG4_07102022.pdf", width = 3, height = 3.5)
sapply(unique(res_deg_4$Gene), function(xx){
  if(toupper(xx) %in% genes2use){
  vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                    which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
                assay = "MAGIC_RNA",
                features = xx,
                # group.by = "new_labels_edited",
                pt.size = 0,
                split.by = "orig.ident",
                cols = c("gray", "orange")) +
    theme(axis.text.x = element_blank())
  print(vv)
  }
})
dev.off()

pdf("Violins_CD8_effector_memory_DEG5_07102022.pdf", width = 3, height = 3.5)
sapply(unique(res_deg_5$Gene), function(xx){
  if(toupper(xx) %in% genes2use){
    vv <- VlnPlot(subset(sc_harmony, cells = unique(c(which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                                      which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated")))),
                  assay = "MAGIC_RNA",
                  features = xx,
                  pt.size = 0,
                  split.by = "orig.ident",
                  cols = c("gray", "orange")) +
      theme(axis.text.x = element_blank())
    print(vv)
  }
})
dev.off()

# redo heatmap 
sc_harmony$new_labels_edited_treatment <- paste0(sc_harmony$new_labels_edited, "-", sc_harmony$orig.ident)
table(sc_harmony$new_labels_edited_treatment)
str_to_sentence(genes2use)[-which(str_to_sentence(genes2use) %in% rownames(sc_harmony))]
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                                         "CD8_Effector_Memory",
                                                                                                         "CD8_Naive_CM",
                                                                                                         "CD8_Memory_Exhausted",
                                                                                                         "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")

# check email Allart 07.11.2022 ####
# run DEG4 and change threshold logFC
for(ii in 1:length(cell_pop)){
  
  message(paste("Working with", cell_pop[ii]))
  temp.cells <- subset(sc_harmony, cells=which(sc_harmony$new_labels_edited == cell_pop[ii]))
  DefaultAssay(temp.cells) <- "RNA"
  
  # dge_4
  temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated" & temp.cells$TRBV13_2and3_pos == "YES", "Treated",
                                   ifelse(temp.cells$orig.ident == "CTR", "CTR", "remove"))
  
  # dge_5
  # temp.cells$temp_groups <- ifelse(temp.cells$orig.ident == "Treated", "Treated",
  #                                   ifelse(temp.cells$orig.ident == "CTR", "CTR", "remove"))
  
  table(temp.cells$temp_groups)
  temp.cells <- subset(temp.cells, cells=which(temp.cells$temp_groups %in% c("CTR", "Treated")))
  
  Idents(temp.cells) <- as.character(temp.cells$temp_groups)
  tt_temp <- table(Idents(temp.cells))
  print(tt_temp)
  
  # find DEGs
  message("Finding DEGs")
  temp.res <- tryCatch(FindMarkers(temp.cells, 
                                   min.cells.group = 5,
                                   ident.1 = "Treated", ident.2 = "CTR", 
                                   logfc.threshold = 0.1, 
                                   min.pct = 0.1),
                       error = function(xx){
                         message(xx)
                         dummy_df <- data.frame(p_val = rep(NA,2),
                                                avg_log2FC = rep(NA,2), 
                                                pct.1 = rep(NA,2),
                                                pct.2  = rep(NA,2),
                                                p_val_adj = rep(NA,2))
                         return(dummy_df)
                       })
  
  temp.res$Gene <- rownames(temp.res)
  # temp.res$FDR <- p.adjust(temp.res$p_val)
  
  if(ii == 1){
    res_deg_trbvs <- list(temp.res)
    names(res_deg_trbvs) <- cell_pop[ii]
  } else {
    res_deg_trbvs[[ii]] <- temp.res
    names(res_deg_trbvs)[ii] <- cell_pop[ii]
  }
  message(paste(" >> Done for", cell_pop[ii]))
}

# check Bcl2
rbindlist(res_deg_trbvs, use.names = T, idcol = T) %>% filter(Gene %in% c("Bcl2", "Ly6c2", "Prdm1") & `.id` == "CD8_Effector_Memory")

VlnPlot(subset(sc_harmony, 
               cells = unique(c(which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "CTR"),
                                +                                             which(sc_harmony$new_labels_edited %in% c("CD8_Effector_Memory") & sc_harmony$orig.ident == "Treated" & sc_harmony$TRBV13_2and3_pos == "YES")))),
        assay = "RNA", 
        # slot = "counts",
        # log = T,
        features = "Bcl2",
        pt.size = 0,
        split.by = "orig.ident",
        cols = c("gray", "orange")) +
  theme(axis.text.x = element_blank())

# address Jonathan email 07.11.2022 ####
# re-label CD4 T helper/CD8 Effector > CD4 Effector
sc_harmony$new_labels_edited <- ifelse(sc_harmony$new_labels_edited == "CD4_T_helper" & sc_harmony$new_labels == "CD8_Effector", "CD4_Effector", sc_harmony$new_labels_edited)
table(sc_harmony$new_labels, sc_harmony$new_labels_edited)
Idents(sc_harmony) <- sc_harmony$new_labels_edited
DimPlot(sc_harmony, 
        # group.by = "new_labels",
        # label = T,
        # label.color = "white",
        # repel = T,
        # label.box = T,
        cols = c("CD4_Treg"  = "green",
                 "CD8_Effector_Exhauted"  = c25[4],
                 "CD8_Effector_Memory"  = "#008080",
                 "CD8_Effector" = "black",
                 "CD8_Memory_Exhausted"  = "red",
                 "CD8_Naive_CM"  = "orange",
                 "CD4_CM"  = c25[1],
                 "CD4_T_helper"  = "purple",
                 "CD8_Exhausted"  = c25[7],
                 "CD8_Precursor_Exhausted"  = c25[10],
                 "Monocytes" = c25[11],
                 "CD4_Effector_Memory" = c25[15],
                 "CD4_Precursor_Exhausted" = c25[13],
                 "CD4_Effector" = c25[22])) +
        facet_wrap(~ sc_harmony$new_labels_edited, ncol=5)
  
VlnPlot(sc_harmony,
        features = c("Cd4", "Cd8a","Cd3e", "Cd63"),
        group.by = "new_labels_edited",
        stack = T, flip = T,
        pt.size = 0)

# get top markers for cells with new annotation
all_markers_new_labels <- FindAllMarkers(sc_harmony,
                              only.pos = T,
                              min.pct = 0.25,
                              logfc.threshold = 0.3)
head(all_markers_new_labels)
dim(all_markers_new_labels %>%
       filter(p_val_adj < 0.05))
all_markers_new_labels %>%
  filter(cluster == "Monocytes") %>%
  # filter(gene == "Il2ra") %>%
  arrange(desc(avg_log2FC)) %>%
  slice(seq_len(10)) %>%
  View
write.table(all_markers_new_labels, "cell_types_top_markers_up_07132022.txt", sep="\t", quote=F, row.names = F)


unique_pop <- unique(sc_harmony$new_labels_edited)
for(ii in 1:length(unique_pop)){
  temp_genes <- all_markers_new_labels %>%
    filter(cluster == unique_pop[ii]) %>%
    filter(p_val_adj <= 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    slice(seq_len(7))
  temp_genes <- temp_genes$gene
  temp_mat <- data.matrix(sc_harmony@assays$SCT@counts[which(rownames(sc_harmony) %in% temp_genes),])
  
  if(ii == 1){
    markers_mat <- temp_mat
  } else {
    markers_mat <- rbind(markers_mat, temp_mat)
  }
  message(paste("Done for", unique_pop[ii]))
}
head(markers_mat)
markers_mat <- markers_mat[-which(duplicated(rownames(markers_mat))),]
markers_mat[1:10,1:5]

cls <- colnames(markers_mat)
heat <- apply(markers_mat,1,scale)
rownames(heat) <- cls
head(heat)
colnames(sc_harmony@meta.data)
df_heat <- sc_harmony@meta.data[,c("orig.ident", "new_labels_edited")]
colnames(df_heat) <- c("Treatment", "CellType")
df_heat <- df_heat[order(df_heat$CellType, decreasing = T),]
head(df_heat)
heat <- merge(df_heat, heat, by="row.names")
heat$Treatment <- NULL

heat <- aggregate(heat[,-c(1,2)], by=list(heat$CellType), FUN=mean)
rownames(heat) <- heat$Group.1
heat$Group.1 <- NULL
heat <- t(heat)
Heatmap(heat,
        cluster_columns = F,
        cluster_rows = T,
        name = "Expression")


# address points in meeting 07132022 and email 07182022 ####
# Rename monocytes to “Myeloid lineage”
sc_harmony$new_labels_edited <- ifelse(sc_harmony$new_labels_edited == "Monocytes", "Myeloid Lineage", sc_harmony$new_labels_edited)
df_heat$CellType <- ifelse(df_heat$CellType == "Monocytes", "Myeloid Lineage", df_heat$CellType)

# Divide treated and untreated in slide 5
table(sc_harmony$orig.ident)
markers_mat[1:10,1:5]
heat <- apply(markers_mat,1,scale)
rownames(heat) <- cls
head(heat)

heat <- merge(df_heat, heat, by="row.names")
heat_ctr <- heat[which(heat$Treatment == "CTR"),]
heat_ctr$Treatment <- NULL

heat_ctr <- aggregate(heat_ctr[,-c(1,2)], by=list(heat_ctr$CellType), FUN=mean)
rownames(heat_ctr) <- heat_ctr$Group.1
heat_ctr$Group.1 <- NULL
heat_ctr <- t(heat_ctr)
h_ctr <- Heatmap(heat_ctr[,c(8,7,2,4,6,5,10,14,13,12,11,9,3,1)],
                 cluster_columns = F,
                 cluster_rows = T,
                 name = paste("Expression",
                              "CTR (left)", sep = "\n"))


heat_t <- heat[which(heat$Treatment == "Treated"),]
heat_t$Treatment <- NULL

heat_t <- aggregate(heat_t[,-c(1,2)], by=list(heat_t$CellType), FUN=mean)
rownames(heat_t) <- heat_t$Group.1
heat_t$Group.1 <- NULL
heat_t <- t(heat_t)
h_t <- Heatmap(heat_t[,c(8,7,2,4,6,5,10,14,13,12,11,9,3,1)],
                 cluster_columns = F,
                 cluster_rows = T,
                 name = paste("Expression",
                              "Treatment (right)", sep = "\n"))
h_ctr + h_t

# UMAP as requested in email 07182022 ####
# correct the UMAPs
# the DimPlot is not reporting color correctly, do it manually
ggplot(data.frame(UMAP1 = sc_harmony@reductions$umap@cell.embeddings[,1],
                  UMAP2 = sc_harmony@reductions$umap@cell.embeddings[,2],
                  CellType = sc_harmony$new_labels_edited,
                  Treatment = sc_harmony$orig.ident),
       aes(x = UMAP1,
           y = UMAP2,
           col = CellType)) +
  geom_point(size = 0.25) +
  theme_classic() +
  scale_color_manual(values = c("CD4_Treg"  = "green",
                                "CD8_Effector_Exhauted"  = c25[4],
                                "CD8_Effector_Memory"  = "#008080",
                                "CD8_Effector" = "black",
                                "CD8_Memory_Exhausted"  = "red",
                                "CD8_Naive_CM"  = "orange",
                                "CD4_CM"  = c25[1],
                                "CD4_T_helper"  = "purple",
                                "CD8_Exhausted"  = c25[7],
                                "CD8_Precursor_Exhausted"  = c25[10],
                                "Myeloid Lineage" = c25[11],
                                "CD4_Effector_Memory" = c25[15],
                                "CD4_Precursor_Exhausted" = c25[13],
                                "CD4_Effector" = c25[22])) +
  facet_wrap(~ Treatment) +
  guides(colour = guide_legend(override.aes = list(size=3), ncol = 1))


# heatmaps as requested in email 07182022 ####
# only selected cell types
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                               "CD8_Effector_Memory",
                                                                                               "CD8_Naive_CM",
                                                                                               "CD8_Memory_Exhausted",
                                                                                               "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")

# all cells
heat <- AverageExpression(sc_harmony, 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- gsub("Monocytes","Myeloid Lineage",colnames(heat))
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat[,c(17:18,13:14,11:12,15:16,19:20,7:8,1:2,3:4,5:6,9:10,21:22,23:24,25:26)],
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType  = c("CD4_Treg"  = "green",
                                                                    "CD8_Effector_Exhauted"  = c25[4],
                                                                    "CD8_Effector_Memory"  = "#008080",
                                                                    "CD8_Effector" = "black",
                                                                    "CD8_Memory_Exhausted"  = "red",
                                                                    "CD8_Naive_CM"  = "orange",
                                                                    "CD4_CM"  = c25[1],
                                                                    "CD4_T_helper"  = "purple",
                                                                    "CD8_Exhausted"  = c25[7],
                                                                    "CD8_Precursor_Exhausted"  = c25[10],
                                                                    "Myeloid Lineage" = c25[11],
                                                                    "CD4_Effector_Memory" = c25[15],
                                                                    "CD4_Precursor_Exhausted" = c25[13],
                                                                    "CD4_Effector" = c25[22]))),
        name = "Expression")

# redo the above with the top genes per each cell type
rownames(markers_mat)
heat <- AverageExpression(sc_harmony, 
                          assays = "MAGIC_RNA",
                          features = rownames(markers_mat), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- gsub("Monocytes","Myeloid Lineage",colnames(heat))
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat[,c(7:8,11:12,19:20,13:14,17:18,15:16,25:26,1:2,3:4,23:24,9:10,5:6,21:22)],
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType  = c("CD4_Treg"  = "green",
                                                                    "CD8_Effector_Exhauted"  = c25[4],
                                                                    "CD8_Effector_Memory"  = "#008080",
                                                                    "CD8_Effector" = "black",
                                                                    "CD8_Memory_Exhausted"  = "red",
                                                                    "CD8_Naive_CM"  = "orange",
                                                                    "CD4_CM"  = c25[1],
                                                                    "CD4_T_helper"  = "purple",
                                                                    "CD8_Exhausted"  = c25[7],
                                                                    "CD8_Precursor_Exhausted"  = c25[10],
                                                                    "Myeloid Lineage" = c25[11],
                                                                    "CD4_Effector_Memory" = c25[15],
                                                                    "CD4_Precursor_Exhausted" = c25[13],
                                                                    "CD4_Effector" = c25[22]))),
        name = "Expression")

# only for the selected cell types
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                               "CD8_Effector_Memory",
                                                                                               "CD8_Naive_CM",
                                                                                               "CD8_Memory_Exhausted",
                                                                                               "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = rownames(markers_mat), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = "Expression")

# count number of cells per cell type based on trbv status
table(sc_harmony$new_labels_edited,
      sc_harmony$TRBV13_2and3_pos,
      sc_harmony$orig.ident)

# address Jonathan email from 07212022 ####
# make heatmap with DEG4 and DEG5
# DEG4
genes2use <- read.table("all_res_deg_4_07172022.txt", sep="\t", header = T)
genes2use <- genes2use %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  filter(`.id` %in% c("CD8_Effector",
                      "CD8_Effector_Memory",
                      "CD8_Naive_CM",
                      "CD8_Memory_Exhausted",
                      "CD4_Treg"))
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                               "CD8_Effector_Memory",
                                                                                               "CD8_Naive_CM",
                                                                                               "CD8_Memory_Exhausted",
                                                                                               "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = unique(genes2use$Gene), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = paste("Expression",
                     "DEG4",
                     "Log2FC > 0.5", sep="\n"))


# DEG5
genes2use <- read.table("all_res_deg_5_07172022.txt", sep="\t", header = T)
genes2use <- genes2use %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
  filter(`.id` %in% c("CD8_Effector",
                      "CD8_Effector_Memory",
                      "CD8_Naive_CM",
                      "CD8_Memory_Exhausted",
                      "CD4_Treg"))
heat <- AverageExpression(subset(sc_harmony, cells = which(sc_harmony$new_labels_edited %in% c("CD8_Effector",
                                                                                               "CD8_Effector_Memory",
                                                                                               "CD8_Naive_CM",
                                                                                               "CD8_Memory_Exhausted",
                                                                                               "CD4_Treg"))), 
                          assays = "MAGIC_RNA",
                          features = unique(genes2use$Gene), 
                          group.by = "new_labels_edited_treatment")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1]))
Heatmap(heat,
        show_column_names = F,
        cluster_columns = F,
        top_annotation = HeatmapAnnotation(df = df_heat,
                                           col = list(Treatment = c("CTR" = "#6495ED",
                                                                    "Treated" = "#800000"),
                                                      CellType = c("CD4_Treg"  = "green",
                                                                   "CD8_Effector_Exhauted"  = c25[4],
                                                                   "CD8_Effector_Memory"  = "#008080",
                                                                   "CD8_Effector" = "black",
                                                                   "CD8_Memory_Exhausted"  = "red",
                                                                   "CD8_Naive_CM"  = "orange"))),
        name = paste("Expression",
                     "DEG5",
                     "Log2FC > 0.5", sep="\n"))




