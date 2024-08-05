################################################################################
##        Title : scRNAseq data analysis (EBs 7g)
##  Description : This script is for creating corresponding Seurat objects and 
##                merge them (no integration included)
##   Researcher : Luis Galán Palma
##       Author : María Maqueda
##         Date : 14th May 2024
################################################################################

# Remark: No samples integration will be conducted here. Integration will be 
#         assessed once clusters are visualized. In this case, samples are merged
#         into one single object. 

################################################################################
## 1. Load packages
################################################################################

require(Seurat)
require(purrr)

################################################################################
## 2. Define paths and samples
################################################################################

# Project path in IMIM cluster where CellRanger outputs is located
cluster.path <- "/Volumes/cancer/CRISPR_screening_mESCs_ABigas/scRNAseq_EBs_7g_vs_WT/cellranger_out/"

# Define all directories where data is stored
samples <- c("WT_s1", "WT_s2", "WT_s3", "g7_s1", "g7_s2", "g7_s3")

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

# $WT_s1
# An object of class Seurat 
# 32285 features across 7757 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $WT_s2
# An object of class Seurat 
# 32285 features across 7214 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $WT_s3
# An object of class Seurat 
# 32285 features across 8711 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $g7_s1
# An object of class Seurat 
# 32285 features across 8785 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $g7_s2
# An object of class Seurat 
# 32285 features across 7601 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts
# 
# $g7_s3
# An object of class Seurat 
# 32285 features across 7109 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 1 layer present: counts

################################################################################
## 4. Merge individual objects into one (THIS IS NOT INTEGRATION!)
##    and store it into a RDS
################################################################################

all.merged.samples <- merge(x = seurat_list[[1]],
                            y = seurat_list[2:length(seurat_list)], 
                            add.cell.ids=samples, # To differentiate barcodes for specific cells
                            project = "EBs_CRISPRa")

all.merged.samples

# An object of class Seurat 
# 32285 features across 47177 samples within 1 assay 
# Active assay: RNA (32285 features, 0 variable features)
# 6 layers present: counts.WT_s1, counts.WT_s2, counts.WT_s3, counts.g7_s1, counts.g7_s2, counts.g7_s3

# Finally, add a new metadata variable indicating the cell condition:WT or g7
all.merged.samples$Group <- sapply(all.merged.samples$orig.ident, function(s) unlist(strsplit(s, split="_"))[1])

################################################################################
## 5. Store merged object
################################################################################

saveRDS(object = all.merged.samples, file = "RDS/All_samples_original_data.RDS")
