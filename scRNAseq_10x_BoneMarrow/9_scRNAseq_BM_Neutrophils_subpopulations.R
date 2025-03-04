################################################################################
##        PAPER : Galan-Palma et al 2024. 
##  Description : This script is for exploring Neutrophils population
##       Author : María Maqueda (mmaqueda@researchmar.net) at BigSpinLab
##                Eric Canton (ecanton@carrerasresearch.org) at BigSpinLab
################################################################################

# In this script, we will annotate different neutrophils subpopulations already 
# described in the literature. After that, a heatmap with markers for each population
# will be generated with cells split by pop.

# Papers considered for exploring neutrophils subpopulation:
# 1) Xie et al Nat Imm 2020 -> here they described the 8 populations (G0-G4 and G5's)
# 2) Kim, Klu and Benayoun Scientific Data 2022 -> Dataset used for annotating 
#    populations. The only drawback for this dataset is that they do not have 
#    G0 or G1 cells (GMP and pre-neutrophils respectively)

################################################################################
## 1. Load packages
################################################################################

require(Seurat)
require(SingleR)
require(scuttle) # To perform normalization and log-transformation as indicated by SingleR
require(scater) # OPTIONAL - to plot some heatmaps as annot diagnostics as in 
# SingleR tutorial
require(openxlsx)
require(ggplot2)
require(dplyr)
require(ComplexHeatmap)

################################################################################
## 2. Reference preparation: 10x Genomics neutrophils dataset
################################################################################

# Neutrophils dataset: Final annotated Seurat object retrieved from (11thOct'24):
# https://figshare.com/articles/dataset/Annotated_Seurat/19623978?file=34856508

load("~/Documents/Projects/Luis_Galán/CRISPR_Screening/scRNAseq_BM/BM_data_ref/Neutrophils_specific/2022-04-12_10x_BM_Ntph_USC_Xie_SingleCellNet_predictions_Seurat_object.RData")
neut.paper <- ntph.singlet.cl

# Prepare dataset as SingleCellExperiments (logNorm and replace gene ids)
neut.sce <- SingleCellExperiment(assays= list(counts= neut.paper@assays$RNA@counts),
                                 colData = neut.paper@meta.data)
neut.sce <- scuttle::logNormCounts(neut.sce, transform="log")

# This dataset contains 6,025 cells and 11,199 genes (3mo mice)

################################################################################
## 3. Own single cell data - prepare SingleCellExperiment
################################################################################

# Import data - minimum dataset (for memory consumption reasons) - this does not
# include dim reduction or transformed data - although it could also be used

# This is a previous RDS version to the one at GEO
own.bm <- readRDS(file = "scRNAseq_CRISPR_Screening_BM_Annotated_Integrated_FINAL_v3.RDS")
own.bm.list <- SplitObject(own.bm, split.by = "Sample_ID")

# Per sample, convert to SCE and logNorm as with the reference

# CD45 BFP negative
bfp_neg <- SingleCellExperiment(assays= list(counts= own.bm.list$CD45_4_5[["RNA"]]$counts.CD45_4_5.1.1.1),  # These are the counts
                                colData = own.bm.list$CD45_4_5@meta.data)
bfp_neg <- scuttle::logNormCounts(x = bfp_neg, 
                                  transform="log")

# CD45 BFP positive
bfp_pos <- SingleCellExperiment(assays= list(counts= own.bm.list$CD45_BFP_4_5[["RNA"]]$counts.CD45_BFP_4_5.2.2.2),  # These are the counts
                                colData = own.bm.list$CD45_BFP_4_5@meta.data)
bfp_pos <- scuttle::logNormCounts(x = bfp_pos, 
                                  transform="log")

own.BM.list.sce <- list("CD45_4_5"  = bfp_neg,
                        "CD45_BFP_4_5" = bfp_pos)

rm(own.bm.list, bfp_neg, bfp_pos) # To free RAM
gc()

# Create a subset just focusing on neutrophils populations
own.bm.neut <- subset(x = own.bm, subset = FINAL_cell.type == c("Neutrophil", "Neutrophil prog"))

################################################################################
## 4. SingleR execution: cell type annotations
################################################################################

# SingleR annotated by cell independently - all cells are considered here but the
# interest is only with Neutrophils and Neutrophils progenitors 

# Execute SingleR: BM neutrophils 3mo mice

preds.singleR.neut <- lapply(own.BM.list.sce, function(sample)
{
  pred <- SingleR::SingleR(test=sample, 
                           de.method = "wilcox",  
                           ref=neut.sce, 
                           labels=neut.sce$SingleCellNet_Xie)
  return(pred)
})

################################################################################
## 5a. Explore results: focus on the already pruned labels!!
################################################################################

preds <- preds.singleR.neut

results <- lapply(preds, function(pred) {
  #  Get the table of predictions
  res <- as.data.frame(table(pred$pruned.labels, useNA="always"))
  colnames(res) <- c("Cell.Type", "Number.Cells")
  return(res)
})


all <- purrr::reduce(.x = results, merge, by = c('Cell.Type'), all = T)
colnames(all) <- c("Cell.Type", paste0("Number.Cells_" , names(results))) 

all$Cell.Type <- as.character(all$Cell.Type)
all$Cell.Type[which(is.na(all$Cell.Type))] <- "UNLABELLED"

################################################################################
## 5b. Annotation diagnostics per sample
################################################################################

# Annotation diagnostics per sample
# In these diagnostics plots: the original labels from SingleR are shown
# Only plot those cells that refer to Neutrophils and Neutrophils progenitors populations

for(i in names(preds))
{
  toplot <- preds[[i]]
  toplot <- toplot[which(rownames(toplot) %in% colnames(own.bm.neut)) , ]
  
  pdf(file = paste0("SingleR_preds_Diagnostics_plots/",i,"_SingleR_Cell_Labels.pdf"), 
      width=12, height=12)
  SingleR::plotScoreHeatmap(toplot,show.pruned = TRUE)
  dev.off()
}

################################################################################
## 6. Aggregate results to complete dataset RDS
################################################################################

# First, put all SingleR labels into a same dataframe together with cell_barcode

cell.annotations <- as.data.frame(do.call(rbind, preds))
cell.annotations$Barcode <- rownames(cell.annotations)
cell.annotations <- cell.annotations[,c("Barcode", "pruned.labels")]

# Does the number of cells match? (They should!!!)
nrow(cell.annotations) == ncol(own.bm)
#[1] TRUE

# Third, add cell annotations label into metadata slot: 
# Add NAs to those cells outside Neutrophils and Neut prog. Additionally,
# Neut Progenitors should be labelled as G0/G1. 

#own.bm$SingleR_Neutrophils <- cell.annotations$pruned.labels

own.bm@meta.data <- seu@meta.data %>%
  mutate(SingleR_Neutrophils = case_when(
    FINAL_cell.type == "Neutrophil" ~ cell.annotations$pruned.labels, 
    FINAL_cell.type == "Neutrophil prog" ~ 'G0/G1', 
    TRUE ~ NA,  # No change for the rest
  ))


# Store results: This is the final RDS version available at GEO
# =============

#saveRDS(object = own.bm, file = "scRNAseq_CRISPR_Screening_BM_Annotated_Integrated_FINAL_v4.RDS")

################################################################################
## 7. Plot UMAP with these labels (Potential Suppl Figures)
################################################################################

# table(own.bm$SingleR_Neutrophils)

# G0/G1    G2    G3    G4   G5a   G5b   G5c 
# 416     957  3221  1704  2336   135   528

# Define colors for neutrophils subpopulations
neut.palette <- c("G0/G1" = "coral1",
                     "G2" = "darkseagreen2",
                     "G3" = "lightpink1",
                     "G4" = "turquoise2",
                     "G5a" = "khaki",
                     "G5b" = "firebrick",
                     "G5c" = "forestgreen")

# Define colors for conditions
cond.colors <- c("CD45_4_5" = "gray80",
                 "CD45_BFP_4_5" = "skyblue3")


# Subset agains to get the label for neutrophils population
own.bm.neut <- subset(x = own.bm, subset = FINAL_cell.type == c("Neutrophil", "Neutrophil prog"))
df.metadata <- own.bm.neut@meta.data
df.metadata$UMAP1_INT <- Embeddings(own.bm.neut, reduction = "umap")[,1]
df.metadata$UMAP2_INT <- Embeddings(own.bm.neut, reduction = "umap")[,2]


# Plot UMAP (integrated) with prediction labels - 
# All cells together - Subset for Neutrophils population

pdf(file = "Suppl_FigS3_UMAP_Integrated_coloured_Neutrophil_type.pdf", 
    width=6,  height=5) # Values in inches

ggplot(df.metadata,aes(x=UMAP1_INT,y=UMAP2_INT, 
                       color=Sample_ID
)) +
  # facet_wrap(~Sample_ID) +
  geom_point(color= "gray20", size=0.5) +
  geom_point(aes(color=SingleR_Neutrophils),size=0.2,alpha = 0.7) +
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "right",  
        legend.spacing.x = unit(0, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=10, face = "bold"),
  ) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(values = neut.palette) 

dev.off()


# Plot UMAP (Non-integrated) with prediction labels - split by Sample
# G3 is the neut population that clearly differs

pdf(file = "Suppl_FigS3_UMAP_NON_Integrated_coloured_Neutrophil_type.pdf", 
    width=10,  height=5) # Values in inches

ggplot(df.metadata,aes(x=UMAP1_NonIn,y=UMAP2_NonIn, 
                       color=SingleR_Neutrophils
)) +
  facet_wrap(~Sample_ID) +
  geom_point(color= "gray20", size=0.5) +
  geom_point(aes(color=SingleR_Neutrophils),size=0.2,alpha = 0.7) +
  # geom_point(aes(color=Sample_ID),size=0.2,alpha = 0.7) +
  theme_void() +  # Empty theme (no axis at all)
  theme(legend.position = "right",  
        legend.spacing.x = unit(0, 'cm'),
        legend.title = element_blank(),
        legend.text = element_text(size=10, face = "bold"),
  ) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  scale_color_manual(values = neut.palette) 

dev.off()

################################################################################
## ANNEX: Plot UMAP with neutrophil markers from paper
################################################################################

# ==================================
g0_markers <- c("Cd34", "Sox4")
FeaturePlot(own.bm.subset, 
            features = g0_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g1_markers <- c("Npm1", "Elane", "Prtn3", "Mpo", "Ctsg") # Can also be expressed by G0

FeaturePlot(own.bm.subset, 
            features = g1_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g2_markers <- c("Fcnb", "Tuba1b")
FeaturePlot(own.bm.subset, 
            features = g2_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

FeaturePlot(own.bm.subset, 
            features = g2_markers, split.by="Sample_ID") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 
# ==================================
g2_g3_markers <- c("Chil3", "Camp", "Ngp","Ltf")
FeaturePlot(own.bm.subset, 
            features = g2_g3_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g3_g4_markers <- c("Mmp8", "Retnlg")
FeaturePlot(own.bm.subset, 
            features = g3_g4_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

VlnPlot(own.bm.subset, 
        features = g3_g4_markers, split.by="Sample_ID") &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g4_all_g5_markers <- c("Cxcl2","Ccl6", "Stfa2l1")
FeaturePlot(own.bm.subset, 
            features = g4_all_g5_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g5c_markers <- c("Gngt2", "Gm2a", "Fgl2")
FeaturePlot(own.bm.subset, 
            features = g5c_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 

# ==================================
g5b_markers <- c("Ifit3", "Rsad2", "Isg15")
FeaturePlot(own.bm.subset, 
            features = g5b_markers) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu"))) 





