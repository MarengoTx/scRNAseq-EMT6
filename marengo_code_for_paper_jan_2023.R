# load libraries
library(Seurat)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)

# load colors
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


# load data
# the code assumes that the data files ar put in your current working directory
sc_data <- readRDS("Marengo_newID_March242023.rds")

# label cells as positive for TRBV13
sc_data$DP <- sc_data@assays$MAGIC_RNA@data["Trbv13-2",] + sc_data@assays$MAGIC_RNA@data["Trbv13-3",]

# pick the values used to set the color scale
sc_data$Trvb5 <- sc_data@assays$MAGIC_RNA@data["Trbv5",]
sc_data$Trvb132_133 <- sc_data@assays$MAGIC_RNA@data["Trbv13-2",] + sc_data@assays$MAGIC_RNA@data["Trbv13-3",]


## Fig 5A
jpeg("UMAP_figure5A.jpeg",height=6,width=10,units="in",res=600)
ggplot(sc_data@meta.data,
       aes(x = sc_data@reductions$umap@cell.embeddings[,"UMAP_1"],
           y = sc_data@reductions$umap@cell.embeddings[,"UMAP_2"],
           col = new_names_mar2023)) +
  geom_point(size = 0.75) +
  theme_classic() +
  scale_color_manual(values = c("CD4_Treg"  = "green",
                                "CD8_Exhausted_2"  = c25[4],
                                "CD8_Effector_Memory"  = "#008080",
                                "CD8_Effector" = "black",
                                "CD8_Exhausted_1"  = "red",
                                "CD8_Naive"  = "orange",
                                "CD4_CM"  = c25[1],
                                "CD4_T_helper"  = "purple",
                                "CD8_Exhausted_3"  = c25[7],
                                "CD8_Precursor_Exhausted"  = c25[10],
                                "Other" = c25[11],
                                "CD4_Effector_Memory" = c25[15],
                                "CD4_Precursor_Exhausted" = c25[13],
                                "CD4_Effector" = c25[22])) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 1)) +
  labs(col = "CellType") +
  ylab("UMAP 2") + xlab("UMAP 1") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10)) +
  facet_wrap(~ orig.ident)
dev.off()


## Fig 5B
FeaturePlot(sc_data,
            features = "DP",
            col=c("lightgray", "red"),
            # order = T, 
            min.cutoff = "q15") +
  facet_wrap(~ sc_data$orig.ident) +
  ggtitle("Trbv13-2 + Trbv13-3")


## Fig 5C
sc_data@meta.data %>%
  dplyr::filter(grepl("^CD", new_labels_edited)) %>%
  group_by(new_labels_edited, orig.ident) %>%
  tally() %>%
  ungroup() %>%
  group_by(orig.ident) %>%
  mutate(perc_TILs = (n/sum(n))*100) %>% 
ggplot(.,
       aes(x = new_labels_edited,
           y = perc_TILs,
           fill = factor(orig.ident, levels = c("Treated", "CTR")))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, col = "black") +
  theme_classic() +
  # ylim(0, 30) +
  coord_flip() +
  scale_fill_manual(values = c("CTR" = "gray", 
                               "Treated" = "darkorange")) +
  labs(fill = "Group") + xlab("") + ylab("% of TILs") +
  theme(axis.text = element_text(size = 12))


## Fig 5D
# Not prepared with R

## Fig 5E
# read DEG data
deg <- read.table("all_res_deg_for_heat_updated_march2023.txt", sep="\t", header = T) 

jpeg("Figure5E_heatmap.jpeg",height=6,width=5,unit="in",res=600)
as.data.frame.matrix(table(deg$.id,
                            deg$Direction)) %>%
     mutate(NS = NULL) %>%
     filter(rownames(.) != "Myeloid Lineage") %>%
     Heatmap(.,
             cluster_rows = T,
             cluster_columns = F,
             col = colorRamp2(breaks = c(0,10,50,150),
                              colors = c("#F5F5F5", "#DCDCDC", "#696969", "black")),
             border = T, 
             name = "N")
dev.off()

## Fig 5F
# read genes
genes2use <- readxl::read_xlsx("genes_for_heatmap_fig5F.xlsx", sheet = 1, col_names = T) %>% pull(Gene)
row_ind <- which(rownames(sc_data@assays$RNA@data) %in% genes2use)

ind_samples_deg4 <- c(which(sc_data$new_names_mar2023 %in% c("CD4_T_helper",
                                                                "CD8_Exhausted_1",
                                                                "CD8_Effector",
                                                                "CD8_Effector_Memory") & sc_data$orig.ident == "CTR"),
                      which(sc_data$new_names_mar2023 %in% c("CD4_T_helper",
                                                                "CD8_Exhausted_1",
                                                                "CD8_Effector",
                                                                "CD8_Effector_Memory") & sc_data$orig.ident == "Treated" & sc_data$TRBV13_2and3_pos == "YES"))

heat <- AverageExpression(subset(sc_data, cells = ind_samples_deg4), 
                          assays = "MAGIC_RNA",
                          features = str_to_sentence(genes2use), 
                          group.by = "new_labels_edited_treatment_mar2023")[[1]]
cls <- colnames(heat)
heat <- t(apply(heat,1,scale))
colnames(heat) <- cls
df_heat <- data.frame(Treatment = sapply(cls, function(xx) strsplit(xx, "-")[[1]][2]),
                      CellType = sapply(cls, function(xx) strsplit(xx, "-")[[1]][1])) %>%
  dplyr::mutate(Treatment = ifelse(Treatment == "Treated", "Treated+TRBVpos", "CTR"))
ha <- HeatmapAnnotation(df = df_heat,
                        col = list(Treatment = c("CTR" = "#6495ED",
                                                 "Treated+TRBVpos" = "#800000"),
                                   CellType = c("CD8_Effector_Memory"  = "#008080",
                                                "CD8_Effector" = "black",
                                                "CD8_Exhausted_1"  = "red",
                                                "CD4_T_helper"  = "purple")),
                        annotation_legend_param = list(CellType = list(at = c("CD4_T_helper",
                                                                              "CD8_Exhausted_1",
                                                                              "CD8_Effector",
                                                                              "CD8_Effector_Memory"))))

ht <- Heatmap(heat,
              show_column_names = F,
              cluster_columns = F,
              top_annotation = ha,
              name = "Expression")
ht@column_order <- c(1:2,7:8,5:6,3:4)
ht

## Fig 5G
genes4violin <- c("Plac8","Gzma","Gzmb", "Ctla2a","Nr4a1","Il2ra", "Zeb2", "Ly6c1","Zfp36","Zfp36l1","Cd6","Tox","Klf4","Rgs16")
ind_samples <- c(which(sc_data$new_labels_edited %in% c("CD8_Effector_Memory") & sc_data$orig.ident == "CTR"),
                      which(sc_data$new_labels_edited %in% c("CD8_Effector_Memory") & sc_data$orig.ident == "Treated" & sc_data$TRBV13_2and3_pos == "YES"))
df_vln <- data.frame(t(sc_data@assays$MAGIC_RNA@data[which(rownames(sc_data@assays$MAGIC_RNA@data) %in% genes4violin),
                                                        ind_samples])) %>%
  dplyr::mutate(ID = rownames(.)) %>%
  dplyr::left_join(sc_data@meta.data %>%
                     dplyr::select(orig.ident, new_labels_edited, TRBV13_2and3_pos) %>%
                     dplyr::mutate(ID = rownames(.)), 
                   by = "ID") %>%
  tidyr::pivot_longer(cols = genes4violin)

for(ii in 1:length(genes4violin)){
  
  temp_gn <- genes4violin[ii]
  
  temp_df <- df_vln %>%
    dplyr::filter(name == temp_gn)
  
  gg <- ggplot(temp_df,
               aes(x = orig.ident,
                   y = value,
                   fill = orig.ident)) +
    geom_violin(scale = "width") +
    theme_classic() +
    scale_fill_manual(values = c("CTR" = "lightgray",
                                 "Treated" = "orange")) +
    ylab("") + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text = element_text(size = 12),
          legend.position = "none") +
    labs(fill = "Group") + ggtitle(temp_gn)
  print(gg)
  
}


## Fig 6A
ggplot(sc_data@meta.data,
       aes(x = sc_data@assays$MAGIC_RNA@data["Cd8a",],
           y = sc_data@assays$MAGIC_RNA@data["Cd4",],
           col = new_labels_edited)) +
  geom_point(size = 0.75) +
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
                                "Other" = c25[11],
                                "CD4_Effector_Memory" = c25[15],
                                "CD4_Precursor_Exhausted" = c25[13],
                                "CD4_Effector" = c25[22])) +
  guides(colour = guide_legend(override.aes = list(size=5), ncol = 1)) +
  labs(col = "CellType") +
  ylab("Cd4") + xlab("Cd8a") +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 10))

## Fig 6B
ggplot(sc_data@meta.data,
       aes(x=sc_data@assays$MAGIC_RNA@data["Trbv13-2",],
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

ggplot(sc_data@meta.data,
       aes(x=sc_data@assays$MAGIC_RNA@data["Trbv13-3",],
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

## Fig 6C

# calculate the delta between pos trbvs and negative trbvs
# scale the values for visualization purposes
sc_data$diff <- scales::rescale(sc_data$Trvb132_133, to = c(0,1)) - scales::rescale(sc_data$Trvb5, to=c(0,1))

ggplot(sc_data@meta.data,
       aes(x = sc_data@reductions$umap@cell.embeddings[,1],
           y = sc_data@reductions$umap@cell.embeddings[,2],
           col=diff)) +
  theme_classic() +
  geom_point(size=0.3) +
  scale_color_gradientn(colors = c("purple", scales::alpha("purple", 0.25), scales::alpha("lightgray", 0.55),scales::alpha("red", 0.25), "red"),
                        # n.breaks = 10) +
                        values = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(col="") +
  facet_wrap(~ orig.ident)


sc_data@meta.data %>%
  dplyr::filter(new_labels == "CD8_Effector_Memory" | new_labels_edited == "CD8_Effector_Memory") %>%
  group_by(new_labels, new_labels_edited) %>%
  tally()


## Supp Fig 7
# the number of labels shown will depend from the size of the window in which the plot is displayed

cell_types_for_volcano <- c("CD8_Effector_Memory", "CD4_Memory_Exhausted", "CD8_Effector", "CD4_T_helper")

for(ii in length(cell_types_for_volcano)){
  
  selectedCellType <- cell_types_for_volcano[ii]
  # selectedCellType <- "CD8_Effector_Memory"
  gg <- res_deg_4 %>%
    dplyr::mutate(Label = ifelse(Gene %in% genes2use, Gene, "")) %>%
    dplyr::rename(CellType = `.id`) %>%
    dplyr::filter(CellType %in% selectedCellType) %>%
    ggplot(.,
           aes(x=avg_log2FC, 
               y=-log10(p_val_adj), 
               col=Direction)) +
    geom_vline(xintercept = 0, size=0.75, col="gray", linetype="dotted") +
    geom_hline(yintercept = -log10(0.05), size=0.75, col="gray", linetype="dotted") +
    geom_point(alpha = 0.7) +
    xlim(-3.5,3.5) +
    theme_minimal() +
    # facet_wrap( ~ CellType, scales = "free") +
    scale_color_manual(values = c("DWN" = "cyan",
                                  "UP" = "tomato",
                                  "NS" = "lightgray")) +
    geom_label_repel(aes(label = Label), color="black",
                     box.padding   = 0.35, 
                     point.padding = 0.5, 
                     segment.color = 'grey50', 
                     max.overlaps = 100) +
    ggtitle(paste(selectedCellType, "DEG4", sep="\n"))
  
  print(gg)
}


#Supplementary Figure 8A:
jpeg("FeatureScatter_Supplementary_Figure8.jpeg",height=8,width=10,units="in",res=600)
FeatureScatter(sc_data,
                feature1 = "Cd8a",
                feature2 = "Cd4", cols = values, 
                group.by = "new_names_mar2023",pt.size=1.5) + guides(color = guide_legend(override.aes = list(size = 5))) + labs(fill = "Cell Type")
dev.off()
