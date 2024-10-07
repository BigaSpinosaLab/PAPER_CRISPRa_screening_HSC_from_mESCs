#------------------------------------------------------------------------------# 
##        PAPER : G. Palma et al 2024. "..."
##  Description : This script is for performing a scRNAseq analysis.
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
##                Eric Canton (ecanton@carrerasresearch.org) at BigSpinLab
#------------------------------------------------------------------------------# 

##### Import libraries #####

# Included any additional libraries required i.e. openxlsx
libraries <- c("ggplot2", "data.table", "dplyr", "stringr", "ggpubr", 
               "Seurat", "glmGamPoi", "clustree", "multtest","metap", 
               "RColorBrewer","ggpubr", "openxlsx", 
               "readxl", "Rmagic", "SeuratWrappers", 
               "patchwork") 
check.libraries <- is.element(libraries, installed.packages()[, 1])==FALSE
libraries.to.install <- libraries[check.libraries]
if (length(libraries.to.install!=0)) {
  require(libraries.to.install)
}

success <- sapply(libraries,require, quietly = FALSE,  character.only = TRUE)
if(length(success) != length(libraries)) {stop("A package failed to return a success in require() function.")}


##### Import Seurat Object #####
# REMARK: This Seurat object already contains an automated cell type labelling with SingleE
seu <- readRDS(file ="/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_SingleR_Annotated.RDS")


##### Normalization, dimensionality reduction (PCA, UMAP) and clustering over WT sample #####

#Prior to normalization, dataset is split into a list of three objects, one per condition:

# split the dataset into a list of two seurat objects (stim and CTRL)
BM.list <- SplitObject(seu, split.by = "Sample_ID")

CD45_4_5 <- BM.list[["CD45_4_5"]]
CD45_BFP_4_5 <- BM.list[["CD45_BFP_4_5"]]


# Normalization aims to remove technical factors such as the library size that 
# may confound the real biological heterogeneity. In this case, during 
# normalization, mitochondrial content is also removed as a  possible source 
# of variation in our samples. SCTransform v2 regularization is used.

# **This is firstly done over the WT sample**

# Normalize
CD45_4_5 <- SCTransform(CD45_4_5,
                        vars.to.regress = "percent.mito",
                        vst.flavor = "v2", # Default value in v5
                        verbose = FALSE)

CD45_4_5 <- RunPCA(CD45_4_5, verbose = FALSE)

pct <- CD45_4_5[["pca"]]@stdev / sum(CD45_4_5[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

n_dim =35 # Practically gets 85% (84.97%)
# This number of dimensions will be used for all the samples


##### Prepare Integration #####

CD45_4_5 <- SCTransform(CD45_4_5, vars.to.regress = "percent.mito", vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs =n_dim, verbose = FALSE)

CD45_BFP_4_5 <- SCTransform(CD45_BFP_4_5, vars.to.regress = "percent.mito", vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs =n_dim, verbose = FALSE)

BM.list <- list(CD45_4_5 = CD45_4_5, CD45_BFP_4_5 = CD45_BFP_4_5)

features <- SelectIntegrationFeatures(object.list = BM.list, nfeatures = 3000)
BM.list <- PrepSCTIntegration(object.list = BM.list, anchor.features = features)

BM.anchors <- FindIntegrationAnchors(object.list = BM.list, normalization.method = "SCT",
                                     anchor.features = features)


##### Integrate Samples #####

BM.combined.sct <- IntegrateData(anchorset = BM.anchors, normalization.method = "SCT")

BM.combined.sct <- RunPCA(BM.combined.sct, verbose = FALSE)

pct <- BM.combined.sct[["pca"]]@stdev / sum(BM.combined.sct[["pca"]]@stdev) * 100

# # Calculate cumulative percents for each PC
cumu <- cumsum(pct)

n_dim = 34 # 85.25% captured variance
BM.combined.sct <- RunUMAP(BM.combined.sct, reduction = "pca",
                           dims = 1:n_dim, verbose = FALSE)

BM.combined.sct <- FindNeighbors(BM.combined.sct, reduction = "pca",
                                 dims = 1:n_dim)

BM.combined.sct <- FindClusters(BM.combined.sct, resolution = seq(0.1, 1, by=0.1))

BM.combined.sct <- Seurat::SetIdent(BM.combined.sct, value = BM.combined.sct$integrated_snn_res.0.4)


##### Subclustering #####
# Apply subclustering in those clusters where more cell type heterogeneity (based on automated labelling) has been detected
# This subclustering will be useful for manual curation of cell types annotation

BM.combined.sct <- Seurat::SetIdent(BM.combined.sct, value = BM.combined.sct$integrated_snn_res.0.4)

BM.combined.sct <- FindSubCluster(BM.combined.sct, 11, "integrated_snn", subcluster.name = "integrated_subclust11_snn_res.0.4",  resolution = 1,   algorithm = 1)
BM.combined.sct <- FindSubCluster(BM.combined.sct, 12, "integrated_snn", subcluster.name = "integrated_subclust12_snn_res.0.4",  resolution = 0.5, algorithm = 1)
BM.combined.sct <- FindSubCluster(BM.combined.sct, 13, "integrated_snn", subcluster.name = "integrated_subclust13_snn_res.0.4",  resolution = 1,   algorithm = 1)
BM.combined.sct <- FindSubCluster(BM.combined.sct, 15, "integrated_snn", subcluster.name = "integrated_subclust15_snn_res.0.4",  resolution = 0.1, algorithm = 1)

##### Data Imputation (MAGIC) #####
## Smoothing and imputation is applied over normalized data - not integrated
BM.combined.sct <- magic(BM.combined.sct, assay = "SCT")

##### Combined subclustering annotation #####

# Each subcluster is found in a different column, so we create a new one to store all subclusters
# We create a dummy variable to store metadata from the Seurat object.
meta <- BM.combined.sct@meta.data

meta$subcluster.integrated_snn_res.0.4 <- as.character(meta$integrated_snn_res.0.4)
meta$subcluster.integrated_snn_res.0.4 <- ifelse((meta$integrated_snn_res.0.4 == "11"), as.character(meta$integrated_subclust11_snn_res.0.4), as.character(meta$subcluster.integrated_snn_res.0.4))
meta$subcluster.integrated_snn_res.0.4 <- ifelse((meta$integrated_snn_res.0.4 == "12"), as.character(meta$integrated_subclust12_snn_res.0.4), as.character(meta$subcluster.integrated_snn_res.0.4))
meta$subcluster.integrated_snn_res.0.4 <- ifelse((meta$integrated_snn_res.0.4 == "13"), as.character(meta$integrated_subclust13_snn_res.0.4), as.character(meta$subcluster.integrated_snn_res.0.4))
meta$subcluster.integrated_snn_res.0.4 <- ifelse((meta$integrated_snn_res.0.4 == "15"), as.character(meta$integrated_subclust15_snn_res.0.4), as.character(meta$subcluster.integrated_snn_res.0.4))
meta <- meta[,c(1:9, 15, 10)]

BM.combined.sct@meta.data <- meta


##### Final Clusters #####

# Manual rearrangement of subclusters - this is based on automated annotation and visualization of specific gene markers expression
# We create a dummy variable to store metadata from the Seurat object.
meta <- BM.combined.sct@meta.data

meta$final.clusters.integrated_snn <- as.character(meta$integrated_snn_res.0.4)
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "13_0" | meta$subcluster.integrated_snn_res.0.4 == "13_1" | meta$subcluster.integrated_snn_res.0.4 == "13_3"), "1", as.character(meta$final.clusters.integrated_snn))
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "11_1" | meta$subcluster.integrated_snn_res.0.4 == "11_5"), "17", as.character(meta$final.clusters.integrated_snn))
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "11_3"), "18", as.character(meta$final.clusters.integrated_snn))
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "11_6"), "19", as.character(meta$final.clusters.integrated_snn))
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "13_2" | meta$subcluster.integrated_snn_res.0.4 == "12_0"), "13", meta$final.clusters.integrated_snn)
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "11_2"), "20", as.character(meta$final.clusters.integrated_snn))
meta$final.clusters.integrated_snn <- ifelse((meta$subcluster.integrated_snn_res.0.4 == "15_2" | meta$subcluster.integrated_snn_res.0.4 == "15_1"), "21", meta$final.clusters.integrated_snn)

meta <- meta[,c(4,2,3,6,7,12,13,27,19,32)]

BM.combined.sct@meta.data <- meta

## Save RDS object
saveRDS(BM.combined.sct,"/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_Annotated_Integrated.RDS")
