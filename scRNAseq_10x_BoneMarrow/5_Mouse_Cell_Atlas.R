#------------------------------------------------------------------------------# 
##        PAPER : G. Palma et al 2024. "..."
##  Description : This script is for importing the Mouse Cell Atlas files and 
##                creating the SingleCellExperiment object
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
##                Eric Canton (ecanton@carrerasresearch.org) at BigSpinLab
#------------------------------------------------------------------------------# 

##### Import libraries #####

libraries <- c("SingleCellExperiment", "dplyr", 
               "stringr", "ggpubr", "Seurat") 
check.libraries <- is.element(libraries, installed.packages()[, 1])==FALSE
libraries.to.install <- libraries[check.libraries]
if (length(libraries.to.install!=0)) {
  require(libraries.to.install)
}

success <- sapply(libraries,require, quietly = FALSE,  character.only = TRUE)
if(length(success) != length(libraries)) {stop("A package failed to return a success in require() function.")}


##### Import data #####

# REMARK: Data was directly downloaded (July'24) from https://bis.zju.edu.cn/MCA/gallery.html?tissue=Bone-Marrow
# and imported into R environment. From MCA website > Gallery > Bone-Marrow > Bone-Marrow (adult, with 9049 cells)


# Data frame with expression data (First row contains row numbers, so we do not need to import it)
bm_atlas <- read.csv("/home/user/Documents/Files/Projects/BM_scRNAseq/Bone-Marrow/Bone-Marrow_dge.csv")[-1]

# Genes corresponding to each row
bm_genes <- read.csv("/home/user/Documents/Files/Projects/BM_scRNAseq/Bone-Marrow/Bone-Marrow_gene.csv")

# Cluster for each sample
bm_clusters <- read.csv("/home/user/Documents/Files/Projects/BM_scRNAseq/Bone-Marrow/Bone-Marrow_barcodes_anno.csv")



##### Restructure data #####

rownames(bm_atlas) <- bm_genes$x

colnames(bm_atlas) <- gsub("\\.", "_", colnames(bm_atlas))
bm_clusters$X <- gsub("\\.", "_", bm_clusters$X)

bm_ident <- bm_clusters[,c(1,2)]
colnames(bm_ident) <- c("cellNames","clusters")
bm_ident$cells <- paste(str_sub(bm_ident$cellNames, 0, 12))

samples <- unique(bm_ident$cells)



##### Create Seurat Object from data #####
seurat_CellAtlas = CreateSeuratObject(counts = bm_atlas, 
                                      project="BoneMarrow")

seurat_CellAtlas$cell.ident <- bm_clusters$Idents.pbmc.
seurat_CellAtlas$cell.num <- bm_clusters$pbmc.cluster


# Save Seurat Object
saveRDS(object = seurat_CellAtlas, file = "/home/user/Documents/Files/Projects/BM_scRNAseq/RDS/scRNAseq_CRISPR_Screening_BM_Mouse_Cell_Atlas.RDS")



