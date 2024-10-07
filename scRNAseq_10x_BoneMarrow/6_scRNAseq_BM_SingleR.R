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

# This RDS was generated from script 4 (final step from pre-processing)
seu <- readRDS("/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_After_QC_and_dbl_removal_QC2_Only2Samples.RDS")

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
CD45_4_5$SingleR_cell.type <- predictions_CD45_4_5$labels

CD45_BFP_4_5 <- BM.list[["CD45_BFP_4_5"]]
CD45_BFP_4_5$SingleR_cell.type <- predictions_CD45_BFP_4_5$labels


seu <- merge(CD45_4_5, CD45_BFP_4_5, project = "Combined_BM_CRISPR_Screening")

# REMARK: These cell annotations will be further inspected (manual curation) after samples integration and clustering
# identification is carried out


## Save RDS object
saveRDS(seu,"/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_SingleR_Annotated.RDS")

