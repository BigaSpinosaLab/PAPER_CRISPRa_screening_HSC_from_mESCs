################################################################################
##        Title : scRNAseq data analysis (EBs 7g)
##  Description : This script is for annotating scRNAseq samples to Embryo Mouse
##                Atlas from Pijuan-Sala et al. 2019 Nature
##   Researcher : Luis Galán Palma
##       Author : María Maqueda
##         Date : 25th Juyly 2024
################################################################################

# Regarding Mouse embryo atlas: 350 Whole embryos were dissociated at timepoints 
# between embryonic days (E) 6.5 and 8.5 of development.
# Libraries were generated using the 10x Genomics Chromium platform (v1 chemistry)

# All samples were clustered and annotated into 37 major cell populations
# However, after inspection, some cell types are relabel as Under-represented based
# on the much reduced number of cells detected:
# Cell types labelled as underrepresented: those appearing in < 0.1% of all cells
# (34,950 cells): this means < 35 cells among all samples

# Cell types affected by these are: Caudal epiblast/mesoderm/neurectoderm, ExE
# endoderm, intermediate endoderm, NMP, notochord, parietal endoderm, somitic
# mesoderm and spinal cord

################################################################################
## 1. Load packages
################################################################################

library(MouseGastrulationData) # This is for retrieving the atlas
require(Seurat)
require(scuttle) # To proceed with the normalization as defined by Marioni's
require(SingleR)
require(scater) # OPTIONAL - to plot some heatmaps as annot diagnostics as in 
# SingleR tutorial
require(openxlsx)
require(ggplot2)

################################################################################
## 2. Reference preparation: Mouse Embryo Atlas - prepare SingleCellExperiment
################################################################################

# Atlas data format: Metadata information for all samples included in the dataset
AtlasSampleMetadata

# Following code already executed, load the corresponding RDS object
#sce <- EmbryoAtlasData(samples = c(1:10,12:20, 23:37)) # Mixed gastrulation excluded and sample 11 (authors)

 # For normalization we use scuttle package to use the proper size factors provided by the authors
# sce <- scuttle::logNormCounts(x = sce, 
#                               size.factors= sizeFactors(sce),
#                               transform="log")

 # Associated symbol to rows
# rownames(sce) <- rowData(sce)$SYMBOL # In order to match with our dataset
 # Remove those cells with NA in cell type
# sce.red <- sce[,!is.na(sce$celltype)]

sce.red <- readRDS(file = "RDS/EmbryoAtlas_Reference.RDS")

################################################################################
## 3. Own single cel data - prepare SingleCellExperiment
################################################################################

# Import data - minimum dataset (for memory consumption reasons) - this does not
# include dim reduction or transformed data - although it could also be used

ebs_seu <- readRDS(file = "RDS/scRNAseq_EBs_WT_7g_After_QC_Filt_woRibo.RDS")
EB.list <- SplitObject(ebs_seu, split.by = "orig.ident")

# Per sample, convert to SCE and logNorm as with the reference

EB.list.SCE <- lapply(EB.list, function(sample) {

  sample.sce <- as.SingleCellExperiment(sample)
  
  sample.sce <- scuttle::logNormCounts(x = sample.sce, 
                                      transform="log")
  
  return(sample.sce)
})

rm(ebs_seu, EB.list) # To free RAM
gc()

################################################################################
## 4. SingleR execution: cell type annotations
################################################################################

# Execute SingleR

preds.singleR <- lapply(EB.list.SCE, function(sample)
  {
  pred <- SingleR::SingleR(test=sample, 
                   de.method = "wilcox",  
                   ref=sce.red, 
                   labels=sce.red$celltype)
  return(pred)
})
  
#saveRDS(object = preds.singleR, file = "RDS/scRNAseq_EBs_SingleR_Preds_EmbryoAtlas.RDS")

################################################################################
## 5a. Explore results: focus on the already pruned labels!!
################################################################################

results <- lapply(preds.singleR, function(pred) {
  #  Get the table of predictions
  res <- as.data.frame(table(pred$pruned.labels, useNA="always"))
  colnames(res) <- c("Cell.Type", "Number.Cells")
  return(res)
})

# Not all samples have representation of all 37 cell types in reference!
# WT_s1 WT_s2 WT_s3 g7_s1 g7_s2 g7_s3 
# 37    36    34    34    32    34 

all <- purrr::reduce(.x = results, merge, by = c('Cell.Type'), all = T)
colnames(all) <- c("Cell.Type", paste0("Number.Cells_" , names(results))) 

all$Cell.Type <- as.character(all$Cell.Type)
all$Cell.Type[which(is.na(all$Cell.Type))] <- "UNLABELLED"

################################################################################
## 5b. Simplification of cell types
## There are cell types under-represented in this dataset, to improve visualization
## cell types that are in that situation are labelled as Under-represented. 
## Additionally, some types are aggregated (Blood progenitors 1/2 and Erythroid 1/2/3)
################################################################################

# Cell types labelled as underrepresented: those appearing in < 0.1% of all cells
# (34,950 cells): this means < 35 cells among all samples
cutoff <- ceiling(sum(all[,grep("Number.", colnames(all))],na.rm = TRUE) *0.001) # 35 cells

# Total cells per cell type across all samples
all$Total.cells <- rowSums(all[,grep("Number.", colnames(all))], na.rm = TRUE)

# Add a new column as 'Final.Cell.Annotation' adding  Underrepresented_type 
all$FINAL.CELL.ANNOTATION <- ifelse(all$Total.cells <cutoff, 
                                    'Underrepresented Type', 
                                    all$Cell.Type)

# Simplify Blood progenitors and Erythroid
all$FINAL.CELL.ANNOTATION[grep("Blood", all$Cell.Type)] <- "Blood progenitors"
all$FINAL.CELL.ANNOTATION[grep("Erythroid", all$Cell.Type)] <- "Erythroid"


# How many unlabelled cells we have?
all$Total.cells[which(all$FINAL.CELL.ANNOTATION == "UNLABELLED")] # 466 cells

# How many unlabelled cells we have?
sum(all$Total.cells[which(all$FINAL.CELL.ANNOTATION == "Underrepresented Type")]) # 125 cells

# Store this as an excel file
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "orchid4")
write.xlsx(x = all, 
             file= "PAPER_RELATED/SingleR_related/CellAnnotations_EmbryoAtlas_SingleR_results.xlsx", 
             headerStyle=hs)


# Annotation diagnostics per sample
# In these diagnostics plots: the original labels from SingleR are shown

for(i in names(preds.singleR))
    {
      pdf(file = paste0("PAPER_RELATED/SingleR_related/Diagnostics_plots/",i,"_SingleR_Cell_Labels.pdf"), 
          width=12, height=12)
      SingleR::plotScoreHeatmap(preds.SingleR[[i]],show.pruned = TRUE)
      dev.off()
    }

################################################################################
## 6. Aggregate results to complete dataset RDS
################################################################################


# First, put all SingleR labels into a same dataframe together with cell_barcode

cell.annotations <- as.data.frame(do.call(rbind, preds.singleR))
cell.annotations$Barcode <- rownames(cell.annotations)
cell.annotations <- cell.annotations[,c("Barcode", "pruned.labels")]

# Second, import the complete Seurat object 
ebs_complete <- readRDS(file = "RDS/scRNAseq_EBs_FINAL_v2.RDS")

  # Does the number of cells match? (They should!!!)
nrow(cell.annotations) == ncol(ebs_complete)
#[1] TRUE


# Third, redefine cell type annotations as previously done by adding 
# 'Under-represented type' and aggregating Blood Prog and Erythroid

# Underrepressented Caudal epiblast/mesoderm/neurectoderm, ExE
# endoderm, intermediate endoderm, NMP, notochord, parietal endoderm, somitic
# mesoderm and spinal cord

cell.annotations$pruned.labels.final  <- cell.annotations$pruned.labels

pattern <- "Caudal epiblast|Caudal Mesoderm|Caudal neurectoderm|ExE endoderm|Intermediate mesoderm|NMP|Notochord|Parietal endoderm|Somitic mesoderm|Spinal cord"
cell.annotations$pruned.labels.final[grep(pattern = pattern, x = cell.annotations$pruned.labels.final)] <- "Underrepresented Type" 

cell.annotations$pruned.labels.final[grep("Blood", cell.annotations$pruned.labels.final)] <- "Blood progenitors"
cell.annotations$pruned.labels.final[grep("Erythroid", cell.annotations$pruned.labels.final)] <- "Erythroid" 

# Add cell annotations label into metadata slot: Original from SingleR and the FINAL 
# previous considerations
ebs_complete$SingleR_EmbryoAtlas <- cell.annotations$pruned.labels
ebs_complete$SingleR_FINAL_EmbryoAtlas <- cell.annotations$pruned.labels.final

saveRDS(object = ebs_complete, file = "RDS/scRNAseq_EBs_FINAL_v2.RDS")

################################################################################
## 7. Plot UMAP with these labels
################################################################################

# REMARK: Neural crest cells are hard to visualize

# Define colors for cell types
cell.colors = c("Allantois" = "slateblue2",
                   "Anterior Primitive Streak" = "burlywood1",
                   "Blood progenitors" = "rosybrown4", 
                   "Cardiomyocytes" = "orchid3",
                   "Def. endoderm" = "blue",
                   "Endothelium" = "darkorange2",
                   "Epiblast" = "burlywood4",
                   "Primitive Streak" = "burlywood2",
                   "Erythroid" = "red3", 
                   "ExE ectoderm" = "magenta", #"orangered",
                   "ExE mesoderm" = "aquamarine3",
                   "Forebrain/Midbrain/Hindbrain" = "black", #"darkolivegreen1",
                   "Gut" = "lightpink2",
                   "Haematoendothelial progenitors" = "rosybrown1",
                   "Mesenchyme" = "tan3",
                   "Mixed mesoderm" = "plum1",
                   "Nascent mesoderm" = "plum4",
                   "Neural crest" = "magenta",
                   "Paraxial mesoderm" = "paleturquoise4",
                   "PGC" = "gold",
                   "Pharyngeal mesoderm" = "paleturquoise",
                   "Rostral neurectoderm" = "green",
                   "Surface ectoderm" = "orangered4",
                   "Visceral endoderm" = "navy",
                   "Underrepresented Type" = "gray90")


# Plot UMAP with prediction labels - All cells together

pdf(file = "PAPER_RELATED/SingleR_related/UMAP_Complete_SingleR_Cell_Annotations_FINAL.pdf", height=6, width=10)

Seurat::DimPlot(ebs_complete, 
                group.by="SingleR_FINAL_EmbryoAtlas",
                reduction =  "umap", 
                repel = TRUE,
                cols = cell.colors, 
                label = TRUE,
                na.value = "gray90",
                label.size = 3,
                alpha=0.5,
                pt.size = 0.7) &  
  theme(legend.text = element_text(size=10, face = "bold"),
        axis.text.x=element_text( hjust=1, size=15),
        axis.text.y=element_text( hjust=1, size=15),
        axis.title = element_text(size=20,face = "bold")) &
  
  theme_void() 

dev.off()

# Plot UMAP with prediction labels - Two UMAPs - split per condition
pdf(file = "PAPER_RELATED/SingleR_related/UMAP_Split_by_Sample_SingleR_Cell_Annotations_FINAL.pdf", height=7, width=18)

Seurat::DimPlot(ebs_complete, 
                group.by="SingleR_FINAL_EmbryoAtlas", 
                split.by = "orig.ident",
                ncol = 3,
                reduction =  "umap", 
                cols = cell.colors, 
                label = FALSE,
                label.size = 1.5,
                na.value = "gray90",
                alpha=0.5,
                pt.size = 0.15) &  
  theme(legend.text = element_text(size=10, face = "bold"),
        axis.text.x=element_text( hjust=1, size=15),
        axis.text.y=element_text( hjust=1, size=15),
        axis.title = element_text(size=20,face = "bold")) &
  
  theme_void() 

dev.off()
