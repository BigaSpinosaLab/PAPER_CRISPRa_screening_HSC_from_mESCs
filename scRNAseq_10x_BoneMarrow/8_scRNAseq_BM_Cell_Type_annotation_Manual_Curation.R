################################################################################
##        PAPER : Galan-Palma et al 2024. 
##  Description : This script is for manually curating the cell annotation
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
##                Eric Canton (ecanton@carrerasresearch.org) at BigSpinLab
################################################################################

# Cell Type annotation: MIX OF AUTOMATED (SingleR) AND MANUAL ANNOTATION (based 
# on gene markers)

# After final cell annotation: HEATMAP with gene markers expressio is generated
# (Fig3D of the manuscript)

################################################################################
## 1. Load packages
################################################################################

require(Seurat)
require(dplyr)
require(ComplexHeatmap)

################################################################################
## 2. Path and general colors
################################################################################

res.path = "Draft_figures_FINAL/"

# Define colors for conditions
cond.colors <- c("CD45_4_5" = "gray80",
                 "CD45_BFP_4_5" = "skyblue3")

general.palette <- c("B-cell" = "springgreen3",
                     "Erythroid" = "indianred3",
                     "Basophil" = "tan1",
                     "Macrophage" = "turquoise3",
                     "Myeloid" = "lightpink2",
                     "Neutrophil" = "khaki",
                     "Neutrophil prog" = "goldenrod3",
                     "T-cell" = "mediumpurple1", 
                     "ILC" = "hotpink", 
                     "HSPC" = "#113226", 
                     "Undetermined_Macro/Neu" = "gray")

################################################################################
## 3. Import data and define gene markers to be inspected
################################################################################

# Import (final) prior where SingleR labels are stored and final clusters are in 
# metadata
seu <- readRDS("/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_Annotated_Integrated.RDS")
# If just to explore
seu <- readRDS("scRNAseq_CRISPR_Screening_BM_Annotated_Integrated_FINAL_v2.RDS")

seu$Sample_ID <- factor(seu$Sample_ID, levels= c("CD45_4_5", "CD45_BFP_4_5"))

# Define known gene markers for different cell types
neutrophil_prog <- c("Mpo", "Elane", "Prtn3", "Ms4a3")
erythroid <- c("Hbb-bs", "Hba-a1", "Hbq1b", "Hbb-bt", "Hba-a2", "Klf1", "Sox6")
basophil <- c("Cd200r3", "Il4", "Itga2b", "Mcpt8", "Prss34")
hspc <- c("Hlf", "Mecom", "Procr", "Hoxa9", "Mycn")
macrophage <- c("Tnfsf13",  # Paracrine factor upregulated in M2 macrophages
                "Tnfrsf21", # M2
                "Adgre1", #F4/80
                "Ms4a6c",  # This is from MCA - Bone Marrow data
                "Ms4a4a", #Novel cell surface marker for M2 macrophages and plasma cells
                "Irf5","Irf9", "Irf7", "Irf8") # Several interferon (may not be specific for macroph)

ilc <- c("Il2rb", "Il18r1", 
         "Gata3", "Icos",
         "Rorc", "Ahr", "Csf2",
         "Id2", "Il7r")

tcell <- c("Ccr6", "Tbx21", "Il2ra", "Il23r", "Klrb1f", "Rorc")

bcell <- c("Cd80", "Cd38", "Cd86", "Cd19",
           "Ighd", "Ighg3", "Ighm")

neutrophil <- c("Ly6g", "Ncf4", "Camp")


# PHASE 1: CONCLUSIONS FROM AUTOMATED ANNOTATION
# ==============================================

# First pruning #########
# We only keep Myeloid and Neutrophil labels from the result of automated cell type annotation,
# the latter with the exception of Mpo high class since Mpo is clearly not expressed in certain
# groups of cells. Additionally, we restrict this labelling to identified 
# (final.clusters.integrated_snn) clusters: 1, 4, 8, 2 and 10 (mainly myeloid) and 0, 6, 9, 13 and 5
# (mainly neutrophil) which are mainly composed (some cases 100%) by myeloid and neu cells

seu@meta.data <- seu@meta.data %>%
  mutate(INTERMEDIATE_cell.type = case_when(
    cell.ident == "Myeloid cell" & (final.clusters.integrated_snn == "0" | 
                                      final.clusters.integrated_snn == "1" |  
                                      final.clusters.integrated_snn == "2" | 
                                      final.clusters.integrated_snn == "4" | 
                                      final.clusters.integrated_snn == "5" | 
                                      final.clusters.integrated_snn == "6" | 
                                      final.clusters.integrated_snn == "8" | 
                                      final.clusters.integrated_snn == "10" | 
                                      final.clusters.integrated_snn == "13")  ~   "Myeloid", 
    cell.ident == "Neutrophil_Chil3 high"  & (final.clusters.integrated_snn == "0" | 
                                                final.clusters.integrated_snn == "1" |  
                                                final.clusters.integrated_snn == "2" | 
                                                final.clusters.integrated_snn == "4" | 
                                                final.clusters.integrated_snn == "5" | 
                                                final.clusters.integrated_snn == "6" | 
                                                final.clusters.integrated_snn == "8" | 
                                                final.clusters.integrated_snn == "9" | 
                                                final.clusters.integrated_snn == "10" | 
                                                final.clusters.integrated_snn == "13")  ~   "Neutrophil", 
    cell.ident == "Neutrophil_Fcnb high" & (final.clusters.integrated_snn == "0" | 
                                                    final.clusters.integrated_snn == "1" |  
                                                    final.clusters.integrated_snn == "2" | 
                                                    final.clusters.integrated_snn == "4" | 
                                                    final.clusters.integrated_snn == "5" | 
                                                    final.clusters.integrated_snn == "6" | 
                                                    final.clusters.integrated_snn == "8" | 
                                                    final.clusters.integrated_snn == "9" |
                                                    final.clusters.integrated_snn == "10" | 
                                                    final.clusters.integrated_snn == "13")  ~   "Neutrophil", 
    cell.ident == "Neutrophil_Ltf high"  & (final.clusters.integrated_snn == "0" | 
                                                   final.clusters.integrated_snn == "1" |  
                                                   final.clusters.integrated_snn == "2" | 
                                                   final.clusters.integrated_snn == "4" | 
                                                   final.clusters.integrated_snn == "5" | 
                                                   final.clusters.integrated_snn == "6" | 
                                                   final.clusters.integrated_snn == "8" | 
                                                    final.clusters.integrated_snn == "9" |
                                                   final.clusters.integrated_snn == "10" | 
                                                   final.clusters.integrated_snn == "13")  ~   "Neutrophil",
    cell.ident == "Neutrophil_Mmp8 high"  & (final.clusters.integrated_snn == "0" | 
                                                final.clusters.integrated_snn == "1" |  
                                                final.clusters.integrated_snn == "2" | 
                                                final.clusters.integrated_snn == "4" | 
                                                final.clusters.integrated_snn == "5" | 
                                                final.clusters.integrated_snn == "6" | 
                                                final.clusters.integrated_snn == "8" | 
                                                final.clusters.integrated_snn == "9" |
                                                final.clusters.integrated_snn == "10" | 
                                                final.clusters.integrated_snn == "13")  ~   "Neutrophil", 
    #cell.ident == "Neutrophil_Mpo high"   ~   NA, 
    #cell.ident == "B cell"  ~   NA,
    #cell.ident == "Basophil" ~   NA,
    #cell.ident == "Erythroid cell"   ~   NA,
    #cell.ident == "Hematopoietic stem and progenitor cell"  ~   NA,
    #cell.ident == "Macrophage_Ctss high"  ~   NA,
    #cell.ident == "Macrophage_Fcna high"  ~   NA,
    #cell.ident == "Macrophage_Ms4a6c high"   ~   NA,
    #cell.ident == "Plasmacytoid dendritic cell"   ~   NA,
    #cell.ident == "Pre B cell"   ~   NA,
    #cell.ident == "T cell_Ccl5 high"  ~   NA,
    TRUE ~ NA,  # The rest, to be manually curated, are forced to NA for the final annot
  ))


# Second pruning #########
# For those cells located in previous clusters but not labelled as Myeloid or Neutrophil,
# (so with NA label from previous point)..
# we force them to the 'majority vote' label of each specific cluster. For instance, 
# cluster 1 which is mainly composed by Myeloid (1815 cells, also includes 29 neutrophils) 
# has cells labelled as Basophil (5) or T_cell_Ccl5_high according to SingleR with 
# MCA ref dataset. These cells are relabel to Myeloid

# table(seu$INTERMEDIATE_cell.type,seu$final.clusters.integrated_snn, useNA= "always")
#                   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21    3    4    5    6    7    8    9 <NA>
# Myeloid cell      31 1815  277    0    0   35    0    0    0    0    0    0 1062    0    0    0 1041    0  127    0  450    0    0
# Neutrophil cell 1780   29  169    0    0  126    0    0    0    0    0    0  315    0    0    0   94  556  599    0  127    0    0
# <NA>              25   14    0  114  198   32  281  189   89  100   42   26    9   55   86 1299    1  226    4  593    1  513    0

# This situation is given in clusters:
# 0 which is mainly Neutrophil
# 1 which is mainly Myeloid
# 13 which is mainly Neutrophil
# 2 which is mainly Myeloid
# 4 which is mainly Myeloid
# 5 which is exclusively Neutrophil
# 6 which is mainly Neutrophil
# 8 which is mainly Myeloid

# This second pruning affects a very limited number of cells except for the Neutrophil Mpo
# located in cluster 5 - which are OK (expressing Neutrophil progenitors). 

seu@meta.data <- seu@meta.data %>%
  mutate(INTERMEDIATE_cell.type_v2 = case_when(
    final.clusters.integrated_snn == "0" & is.na(INTERMEDIATE_cell.type) ~ "Neutrophil", 
    final.clusters.integrated_snn == "1" & is.na(INTERMEDIATE_cell.type) ~ "Myeloid", 
    final.clusters.integrated_snn == "13" & is.na(INTERMEDIATE_cell.type) ~ "Neutrophil",  
    final.clusters.integrated_snn == "2" & is.na(INTERMEDIATE_cell.type) ~ "Myeloid", 
    final.clusters.integrated_snn == "4" & is.na(INTERMEDIATE_cell.type) ~ "Myeloid", 
    final.clusters.integrated_snn == "5" & is.na(INTERMEDIATE_cell.type) ~ "Neutrophil", 
    final.clusters.integrated_snn == "6" & is.na(INTERMEDIATE_cell.type) ~ "Neutrophil", 
    final.clusters.integrated_snn == "8" & is.na(INTERMEDIATE_cell.type) ~ "Myeloid",
    TRUE ~ INTERMEDIATE_cell.type,  # No change for the rest
  ))


table(seu$INTERMEDIATE_cell.type_v2, seu$Sample_ID, useNA = "always")

# CD45_4_5 CD45_BFP_4_5 <NA>

# Myeloid        2953         1910    0
# Neutrophil     3238          941    0
# <NA>           1319         2169    0

table(seu$INTERMEDIATE_cell.type_v2, seu$final.clusters.integrated_snn, seu$Sample_ID, useNA = "always")
# , ,  = CD45_4_5
#               0    1   10   11   12   13   14   15   16   17   18   19    2   20   21    3    4    5    6    7    8    9 <NA>
# Myeloid       4 1089  173    0    0   24    0    0    0    0    0    0  683    0    0    0  629    0   75    0  276    0    0
# Neutrophil 1366   17  132    0    0   83    0    0    0    0    0    0  274    0    0    0   77  657  444    0  102   86    0
# <NA>          0    0    0  108   33    0    9   22   64    1   42   24    0    6   12  157    0    0    0  501    0  340    0
# 
# , ,  = CD45_BFP_4_5
# 
# 
#               0    1   10   11   12   13   14   15   16   17   18   19    2   20   21    3    4    5    6    7    8    9 <NA>
# Myeloid      27  740  104    0    0   11    0    0    0    0    0    0  388    0    0    0  413    0   52    0  175    0    0
# Neutrophil  439   12   37    0    0   75    0    0    0    0    0    0   41    0    0    0   17  125  159    0   25   11    0
# <NA>          0    0    0    6  165    0  272  167   25   99    0    2    0   49   74 1142    0    0    0   92    0   76    0

# ==============================================
# PHASE 2: MANUAL ANNOTATION based on lineage specific markers
# ==============================================

DefaultAssay(seu) <- 'MAGIC_SCT' # Visualized expression based on smoothed and imputed data

# We have to review and manually annotate following clusters (the ones that from automated - pruned
# labelling are NAs): 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 3, 7 and 9

# MARKERS FOR NEUTROPHIL PROGENITORS: Mpo, Elane, Prtn3, Ms4a3. --> Expression restricted to Cluster 9
# This conclusion aligns with SingleR automated labelling (Neutrophil Mpo label in this cluster)
FeaturePlot(object = seu, 
            reduction="umap", 
            features = neutrophil_prog)

FeaturePlot(object = seu, 
            reduction="umap", 
            features = neutrophil_prog)

# MARKERS FOR ERYTHROIDS: Klf1 and Sox6 (definitive and mature) and 
# Hbb-bs, Hba-a1, Hbq1b, Hbb-bt and Hba-a2 (hemoglobin genes) --> Expression restricted to Cluster 21
# This conclusion aligns with SingleR automated labelling (Erythroid cell)
FeaturePlot(object = seu, 
            reduction="umap", 
            features = erythroid)

Seurat::VlnPlot(seu, features = erythroid, 
                group.by = "final.clusters.integrated_snn")

# MARKERS FOR BASOPHILS: Mcpt8 and Prss34 (Basophil progenitors) and 
# Cd200r3, Il4 and Itga2b --> Expression restricted to Clusters 18 and 19
# This conclusion aligns with SingleR automated labelling for Cluster 18
FeaturePlot(object = seu, 
            reduction="umap", 
            features = basophil)

Seurat::VlnPlot(seu, features = basophil, 
                group.by = "final.clusters.integrated_snn")


# MARKERS FOR HSPCs: Hlf, Mecom, Procr, Hoxa9 and Mycn --> Expression restricted to Cluster 11
FeaturePlot(object = seu, 
            reduction="umap", 
            split.by = "Sample_ID",
            features = hspc)

Seurat::VlnPlot(seu, features = hspc, 
                split.by = "Sample_ID",
                group.by = "final.clusters.integrated_snn")

# MARKERS FOR Macrophages:  --> Expression restricted to Cluster 7 and 16
# These were already identified with automated labelling 
FeaturePlot(object = seu, 
            reduction="umap", 
            split.by  = "Sample_ID",
            features = macrophage)

Seurat::VlnPlot(seu, features = macrophage, 
                split.by="Sample_ID",
                group.by = "final.clusters.integrated_snn")

# MARKERS FOR ILCs:  --> Expression restricted to Cluster 14 - strong signal is 
# as well detected in cluster 17 (T cell) but in a lower level than for cluster 14
# This cell type was not present in the reference dataset
FeaturePlot(object = seu, 
            reduction="umap", 
            features = ilc)

Seurat::VlnPlot(seu, features = ilc, 
                group.by = "final.clusters.integrated_snn")


# MARKERS FOR T-cells:  --> Expression observed in clusters 14 and 17
# ILCs and T-cells are functionally alike with the key distinguishing factor of
# lacking antigen-specific T-cell receptor TRC in ILCs
# Here we can see that markers expression for T-cells is in both clusters but since
# we observed that in Cluster 17 there was a stronger signal for ILC markers, 
# we annotate Cluster 14 to ILCs and Cluster 17 to T-cell
# ILCs seen as the ultimate invariant T-cells

FeaturePlot(object = seu, 
            reduction="umap", 
            split.by="Sample_ID",
            features = tcell)

Seurat::VlnPlot(seu, features = tcell, 
                group.by = "final.clusters.integrated_snn")

# REMARK: There is  1 cell BFP negative included in this Cluster 17 (following 
# cluster 11 of HSPCs), this is visualized in the integrated and non integrated
# dataset. Moreover, this cell is showing HSPC markers expression -> redefined as
# HSPCs. To be edited at the end with the rest of annotations

# MARKERS FOR B-cells:  --> Expression observed in Clusters 3, 15 and 12
# Single R annotation mainly labeled clusters 3, 15 and 12 as B-cells
# Let's inspect some gene markers here
# Violin plots clearly show how Bcell markers are not expressed for CD45 neg cells
# These cell are wrongly clustered here, most probably, as a result of integration
# SingleR already labelled them as other thing different to Bcells

FeaturePlot(object = seu, 
            reduction="umap", 
            #split.by = "Sample_ID",
            features = bcell)

Seurat::VlnPlot(seu, features = bcell, 
                group.by = "final.clusters.integrated_snn", 
                split.by = "Sample_ID")


# After inspecting gene markers, conclusion is:
# a) Bcells for BFP positive are located in Clusters 3, 12 and 15
# b) For BFP negative, Cluster 3 cells are Macrophages, 
#    Cluster 12 cells Undetermined between Macro/Neu
#    and Cluster 15 Undetermined between Macro/Neu 


# MARKERS FOR NEUTROPHIL (just as a sanity check): Ly6g, Ncf4, Camp and Mmp25
FeaturePlot(object = seu, 
            reduction="umap", 
            #split.by = "Sample_ID",
            features = neutrophil)

Seurat::VlnPlot(seu, features = neutrophil, 
                group.by = "final.clusters.integrated_snn",
                split.by = "Sample_ID")


# REMARK: Cluster 20 and 12/BFP_neg and 15/BFP_neg are labelled as 'Undetermined_Macro/Neu'
# since there is no clear consistency/homogeneity in gene markers expression between 
# Macrophages and Neutrophils markers

# ==============================================
# FINAL: Include manual annotation in the metadata Seurat object
# ==============================================

seu@meta.data <- seu@meta.data %>%
  mutate(FINAL_cell.type = case_when(
    final.clusters.integrated_snn == "9" & is.na(INTERMEDIATE_cell.type) ~ "Neutrophil prog",
    final.clusters.integrated_snn == "21" ~ "Erythroid", 
    final.clusters.integrated_snn == "18" ~ "Basophil",  
    final.clusters.integrated_snn == "19" ~ "Basophil", 
    final.clusters.integrated_snn == "11" ~ "HSPC", 
    final.clusters.integrated_snn == "7"  ~ "Macrophage", 
    final.clusters.integrated_snn == "16" ~ "Macrophage", 
    final.clusters.integrated_snn == "14"  ~ "ILC",
    final.clusters.integrated_snn == "17" & (Sample_ID == "CD45_BFP_4_5") ~ "T-cell",
    final.clusters.integrated_snn == "17" & (Sample_ID == "CD45_4_5") ~ "HSPC", # This is 1 cell 
    final.clusters.integrated_snn == "20"  ~ "Undetermined_Macro/Neu",
    final.clusters.integrated_snn == "3" & (Sample_ID == "CD45_BFP_4_5") ~ "B-cell",
    final.clusters.integrated_snn == "12" & (Sample_ID == "CD45_BFP_4_5") ~ "B-cell",
    final.clusters.integrated_snn == "15" & (Sample_ID == "CD45_BFP_4_5") ~ "B-cell",
    final.clusters.integrated_snn == "3" & (Sample_ID == "CD45_4_5") ~ "Macrophage",
    final.clusters.integrated_snn == "12" & (Sample_ID == "CD45_4_5") ~ "Undetermined_Macro/Neu",
    final.clusters.integrated_snn == "15" & (Sample_ID == "CD45_4_5") ~ "Undetermined_Macro/Neu",
    TRUE ~ INTERMEDIATE_cell.type_v2,  # No change for the rest
  ))

# Resume results - remove intermediate for final Seurat object
colnames(seu@meta.data)[6] <- "SingleR_cell.type"
seu@meta.data <- seu@meta.data[,-c(7,8,12,13)]


# STORE Final results: this is the final RDS to be shared
saveRDS("scRNAseq_CRISPR_Screening_BM_Annotated_Integrated_FINAL_v2.RDS")


# ==============================================
# HEATMAP showing expression of marker genes along all manually curated clusters:
# 3, 7, 9 (only neutrophil prog!), 11, 12, 14, 15, 16, 17, 18, 19, 20 and 21
# ==============================================

# Prepare list of markers to be plot in the heatmap
all_markers <- data.frame(Markers = c(hspc, bcell, tcell, ilc, basophil,
                                      erythroid,neutrophil_prog, macrophage),
                 CellType = c(
                   rep("HSPC", length(hspc)), 
                   rep("B-cell", length(bcell)), 
                   rep("T-cell", length(tcell)),
                   rep("ILC", length(ilc)),
                   rep("Basophil",  length(basophil)),
                   rep("Erythroid", length(erythroid)),
                   rep("Neutrophil prog", length(neutrophil_prog)),
                   rep("Macrophage", length(macrophage))))

all_markers$CellType <- factor(all_markers$CellType, 
                               levels= unique(all_markers$CellType))

# Obtain exprs matrix from proper assay
exprs <- as.data.frame(GetAssayData(object = seu, 
                                    assay = "MAGIC_SCT"))

metadata <- seu@meta.data
metadata$FINAL_cell.type <- factor(metadata$FINAL_cell.type , 
                               levels= c("HSPC", "B-cell",  "T-cell", "ILC", "Basophil", 
                                         "Erythroid", "Neutrophil prog", "Macrophage",
                                         "Undetermined_Macro/Neu"))

# Cells are order first per cell type and then per condition (BFP negative o BFP positive)
# This will ease the heatmap visualization

metadata.reorder <- arrange(metadata, 
                            FINAL_cell.type, 
                            group_by = Sample_ID)
exprs.reorder <- exprs[,match(rownames(metadata.reorder), colnames(exprs))]

# Identify gene matches: we are only showing selected genes (markers)
matches <- intersect(all_markers$Markers, rownames(exprs.reorder))
exprs.subset <- exprs.reorder[which(rownames(exprs.reorder) %in% matches),]

# In case we are only interested in selected clusters, columns/rows  subset
# SKIP this if all clusters are want to be plotted
# NOTE: Here, we are removing those clusters referring to Neutrophil and Myeloid (most of the samples)
exprs.subset <- exprs.subset[ , which(metadata.reorder$FINAL_cell.type %in% 
                                        c("HSPC", "B-cell",  "T-cell", "ILC", "Basophil", 
                                          "Erythroid", "Neutrophil prog", "Macrophage",
                                          "Undetermined_Macro/Neu"))]
metadata.reorder <- metadata.reorder[which(metadata.reorder$FINAL_cell.type %in% 
                                             c("HSPC", "B-cell",  "T-cell", "ILC", "Basophil", 
                                               "Erythroid", "Neutrophil prog", "Macrophage",
                                               "Undetermined_Macro/Neu")) ,]

# Scale the rows (genes)
exprs.subset <- t(scale(t(exprs.subset)))
# Sort them in the same order
all_markers <- all_markers[match(rownames(exprs.subset), all_markers$Markers),]


# Get the final annotated cell type and condition (sample origin)
cell_type <- metadata.reorder$FINAL_cell.type
condition <- metadata.reorder$Sample_ID

## Color ramp: make the white color map to 0. the red map to highest and the blue map to the lowest
col_fun = circlize::colorRamp2(c(quantile(exprs.subset,0.01), 0, quantile(exprs.subset,0.99)), # To remove extreme cases (default behaviour in ComplexHeatmap)
                               c("skyblue3", "white", "#67001F"))

htmp <- Heatmap(exprs.subset, 
                name = "Scaled Exprs (MAGIC_SCT)",  
                column_split = cell_type,
                row_split = all_markers$CellType,  
                row_gap = unit(1, "mm"),
                cluster_columns = FALSE,
                cluster_column_slices = TRUE,
                column_title_gp = gpar(fontsize = 3, fontface="bold"),
                column_gap = unit(0.5, "mm"),
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                clustering_distance_rows = "pearson",
                row_names_side = "left",
                col = col_fun,
                row_names_gp = gpar(fontsize = 4, fontface="bold"),
                row_title_gp = gpar(fontsize = 3, fontface="bold", fontcolor="grey"), 
                top_annotation = HeatmapAnnotation(CELL_TYPE= cell_type, ORIGIN=condition,
                                                   col = list(ORIGIN = cond.colors, CELL_TYPE = general.palette),
                                                   simple_anno_size=unit(0.2,"cm"),
                                                   show_annotation_name = c(CELL_TYPE=FALSE, ORIGIN=FALSE)),
                show_column_names = FALSE,
                use_raster = TRUE,
                heatmap_legend_param = list(legend_direction = "horizontal"),
                raster_quality = 4)


pdf(file = paste0(res.path,"Fig3D_Heatmap_Gene_Markers_without_Neutrophil_Myeloid.pdf"), 
    width=8,  height=4) # Values in inches

draw(htmp, heatmap_legend_side="bottom", column_title_gp=grid::gpar(fontsize=16))

dev.off()

