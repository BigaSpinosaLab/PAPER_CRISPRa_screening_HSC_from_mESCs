################################################################################
##        Title : scRNAseq data analysis CRISPR Screening BM
##  Description : This script is for creating corresponding Seurat objects and 
##                merge them (no integration included)
##   Researcher : Researcher: Luís Galan
##       Author : María Maqueda
##         Date : 19th June 2024 
################################################################################

# Remark: No samples integration will be conducted here. Integration will be 
#         assessed once clusters are visualized. In this case, samples are merged
#         into one single object. 

# Remark 2: Doublet detection previosly conducted is added to Seurat objects to be
# created here. Final decision about this detection is included here.

################################################################################
## 1. Load packages
################################################################################

require(Seurat)
require(purrr)

################################################################################
## 2. Define paths and samples
################################################################################

# Project path in IMIM cluster where CellRanger outputs is located
cluster.path <- "/run/user/1001/gvfs/smb-share:server=172.20.4.46,share=projectscomput/cancer/CRISPR_screening_mESCs_ABigas/scRNAseq_Sorted_BM/cellranger_out/"

# Define all directories where data is stored
samples <- c("CD45_4_5", "CD45_BFP_3_5", "CD45_BFP_4_5")

# Define all directories
datadirs <- sapply(samples, function(s) paste0(cluster.path, s, "/outs/filtered_feature_bc_matrix"))

# Load results from scDbl FINAL
dbl <- readRDS(file = "RDS/scDblFinder_results_FINAL.RDS")

################################################################################
## 3. Load sparse data matrix and create individual Seurat objects
################################################################################

seurat_list <- purrr::map2(datadirs, samples, function(path,sample){
  data <- Seurat::Read10X(data.dir = path)
  seurat_object = CreateSeuratObject(counts = data, 
                                     project=sample) # Project name is required for further merging
  
  # Add sample ID into metadata
  seurat_object$Sample_ID <- sample
  
  # Add results from scDblFinder (final conclusion) into metadata slot
  scDblFinder.sample <- dbl[which(dbl$Sample %in% sample),]
  seurat_object$scDblFinder.class <- scDblFinder.sample$Final
  
  # Return full seurat object
  return(seurat_object)
})

seurat_list

# $CD45_4_5
# An object of class Seurat 
# 32285 features across 8461 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $CD45_BFP_3_5
# An object of class Seurat 
# 32285 features across 2821 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $CD45_BFP_4_5
# An object of class Seurat 
# 32285 features across 5635 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts


################################################################################
## 4. Merge individual objects into one (THIS IS NOT INTEGRATION!)
##    and store it into a RDS
################################################################################

all.merged.samples <- merge(x = seurat_list[[1]],
                            y = seurat_list[2:length(seurat_list)], 
                            add.cell.ids=samples, # To differentiate barcodes for specific cells
                            project = "CRISPR_Screening_BM")

all.merged.samples

# An object of class Seurat 
# 32285 features across 16917 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 3 layers present: counts.CD45_4_5, counts.CD45_BFP_3_5, counts.CD45_BFP_4_5

################################################################################
## 5. Store merged object
################################################################################

saveRDS(object = all.merged.samples, file = "RDS/scRNAseq_CRISPR_Screening_BM_RAW_original_data.RDS")
