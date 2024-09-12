#------------------------------------------------------------------------------# 
##        PAPER : G. Palma et al 2024. "..."
##  Description : This script  utilizes SingleR to perform unbiased cell type
##                recognition from the BM samples using the Mouse Cell Atlas 
##                dataset as reference
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
##                Eric Canton (ecanton@carrerasresearch.org) at BigSpinLab
#------------------------------------------------------------------------------# 

##### Import data #####

mca <- readRDS("/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_Mouse_Cell_Atlas.RDS")

seu <- readRDS("/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/final/scRNAseq_CRISPR_Screening_BM_Annotated_Integrated_FINAL.RDS")

##### Create SingleCellExperiment #####

###### Mouse Cell Atlas #####
mca <- NormalizeData(mca, normalization.method = "LogNormalize", scale.factor = 10000)
mca <- ScaleData(mca)

mca_sce <- as.SingleCellExperiment(mca)
rowData(mca_sce)$feature_symbol <- rownames(mca_sce)
mca_sce <- selectFeatures(mca_sce)

###### BM CRISPR Screening data ######
# split the dataset into a list of two seurat objects (stim and CTRL)
BM.list <- SplitObject(seu, split.by = "orig.ident")

CD45_4_5 <- BM.list[["CD45_4_5"]]
CD45_4_5 <- NormalizeData(CD45_4_5)
CD45_4_5 <- ScaleData(CD45_4_5)
CD45_4_5 <- as.SingleCellExperiment(CD45_4_5)
rowData(CD45_4_5)$feature_symbol <- rownames(CD45_4_5)
CD45_4_5 <- selectFeatures(CD45_4_5)


CD45_BFP_4_5 <- BM.list[["CD45_BFP_4_5"]]
CD45_BFP_4_5 <- NormalizeData(CD45_BFP_4_5, normalization.method = "LogNormalize", scale.factor = 10000)
CD45_BFP_4_5 <- ScaleData(CD45_BFP_4_5)
CD45_BFP_4_5 <- as.SingleCellExperiment(CD45_BFP_4_5)
rowData(CD45_BFP_4_5)$feature_symbol <- rownames(CD45_BFP_4_5)
CD45_BFP_4_5 <- selectFeatures(CD45_BFP_4_5)


##### Cell Type Recognition using SingleR #####

predictions_CD45_4_5 <- SingleR::SingleR(test=CD45_4_5, ref=mca_sce, labels=mca_sce$cell.ident)
predictions_CD45_BFP_4_5 <- SingleR::SingleR(test=CD45_BFP_4_5, ref=mca_sce, labels=mca_sce$cell.ident)

CD45_4_5 <- BM.list[["CD45_4_5"]]
CD45_4_5$cell.ident <- predictions_CD45_4_5$labels

CD45_BFP_4_5 <- BM.list[["CD45_BFP_4_5"]]
CD45_BFP_4_5$cell.ident <- predictions_CD45_BFP_4_5$labels


seu <- merge(CD45_4_5, CD45_BFP_4_5, project = "Combined_BM_CRISPR_Screening")

###### Reduced cell identity ######
seu$red.cell.ident <- seu$cell.ident

seu@meta.data <- seu@meta.data %>%
  mutate(red.cell.ident = case_when(cell.ident == "Myeloid cell" ~ "Myeloid cell",
                                    cell.ident == "Neutrophil_Mmp8 high" | cell.ident == "Neutrophil_Mpo high" | cell.ident == "Neutrophil_Ltf high" | 
                                      cell.ident == "Neutrophil_Fcnb high" | cell.ident == "Neutrophil_Chil3 high" ~ "Neutrophil",
                                    cell.ident ==  "Macrophage_Ms4a6c high" | cell.ident ==  "Macrophage_Fcna high" | cell.ident == "Macrophage_Ctss high" ~ "Macrophage",
                                    cell.ident ==  "Basophil" ~ "Basophil",
                                    cell.ident ==  "Plasmacytoid dendritic cell" ~ "Plasmacytoid dendritic cell",
                                    cell.ident ==  "Erythroid cell" ~ "Erythroid cell",
                                    cell.ident ==  "T cell_Ccl5 high" ~ "T cell_Ccl5 high",
                                    cell.ident ==  "B cell" ~ "B cell",
                                    cell.ident ==  "Pre B cell" ~ "Pre B cell",
                                    cell.ident ==  "Hematopoietic stem and progenitor cell" ~ "Hematopoietic stem and progenitor cell"))


## Save RDS object
saveRDS(seu,"/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_Annotated_Only2Samples.RDS")

