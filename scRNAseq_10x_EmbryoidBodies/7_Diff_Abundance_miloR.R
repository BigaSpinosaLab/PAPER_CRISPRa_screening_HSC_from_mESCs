
################################################################################
##        Title : scRNAseq data analysis (EBs 7g)
##  Description : This script is for conducting a differential abundance testing 
##                with Milo
##   Researcher : Luis Galán Palma
##       Author : María Maqueda
##         Date : 28th July 2024
################################################################################

# https://marionilab.github.io/miloR/articles/milo_demo.html
# Milo is a tool for analysis of complex single cell datasets generated from 
# replicated multi-condition experiments, which detects changes in composition #
# between conditions. Cluster-free approach

################################################################################
## 1. Load packages
################################################################################

require(miloR)
require(Seurat)
require(SingleCellExperiment)
require(dplyr)
require(ggplot2)
require(scater) # For visualization results (UMAP embedding + DA results)
require(UpSetR)
require(openxlsx)

################################################################################
## 2. Load data
################################################################################

ebs_complete <- readRDS(file = "RDS/scRNAseq_EBs_FINAL_v2.RDS")

# Convert first to SCE and then to Milo object
ebs_sce <- as.SingleCellExperiment(ebs_complete)
ebs_milo <- Milo(ebs_sce)

# To free RM
rm(ebs_complete, ebs_sce)

################################################################################
## 3. Construct KNN graph
################################################################################

ebs_milo <- miloR::buildGraph(ebs_milo, 
                              k = 30,  # Default, number of NN to consider for graph 
                                      # Recommended to use the same as in UMAP (==30)
                              d = 38) # Number of dimensions to use (Default is 50) - We used 38 for UMAP

# Defining representative neighbourhoods

ebs_milo <- miloR::makeNhoods(ebs_milo, 
                              prop = 0.1, # Default value - Proportion of cells to randomly sample to start with
                                          # 0.1 recommended for datasets < 30k (we have 34k - so we leave this)
                              k = 30,  # The same k used to construct the input graph
                              d= 38,  # Default value
                              refined = TRUE) # Run sampling refinement algorithm

# Take a look at how big the neighbourhoods are (how many cells in each neighbor)
# Empyrically, distribution peaking between 50-100 is fine
# I was expecting maybe less since authors indicate a rule of thum of 5x N_samples avg dist (=30 in our case)

pdf(file = "PAPER_RELATED/miloR_related/miloR_distr_neighbourhood_size.pdf", height=4, width=4)
plotNhoodSizeHist(ebs_milo)
dev.off()

# Counting cells in neighbourhoods from each sample

ebs_milo <- countCells(ebs_milo, meta.data = as.data.frame(colData(ebs_milo)), sample="orig.ident")
head(nhoodCounts(ebs_milo))

# 6 x 6 sparse Matrix of class "dgCMatrix"
# WT_s1 WT_s2 WT_s3 g7_s1 g7_s2 g7_s3
# 1    22    13    10     .     2    10
# 2    14    11    14    11    11    12
# 3     4     8     7    12    12     5
# 4     7     3    12    13    13     6
# 5     7    14    14     7    14    11
# 6    10    11    11    13    10    17

# Defining the experimental design: we want to detect DA between g7-activated and WT conditiions
# Stored in the Group column
# No technical covariate to be added

ebs_design <- data.frame(colData(ebs_milo))[,c("orig.ident", "Group")]

ebs_design <- dplyr::distinct(ebs_design) # Otherwise, we have all list of cells
rownames(ebs_design) <- ebs_design$orig.ident

ebs_design

# orig.ident Group
# WT_s1      WT_s1    WT
# WT_s2      WT_s2    WT
# WT_s3      WT_s3    WT
# g7_s1      g7_s1    g7
# g7_s2      g7_s2    g7
# g7_s3      g7_s3    g7

# Computing neighbourhood connectivity: to correct p-values accounting for the amount
# of overlap between neighbourhoods. 
# IMPORTANT!!!!! This execution takes some time (30min approx)
ebs_milo <- miloR::calcNhoodDistance(ebs_milo, 
                              d=38)

# Testing
ebs_design$Group <- relevel(as.factor(ebs_design$Group),ref="WT") #To ensure test g7 vs WT IMPORTANT!!

da_results <- testNhoods(ebs_milo, 
                         design = ~  Group, 
                         design.df = ebs_design)

head(da_results) # Important metrics to look at:
# logFC indicates logFC in cell numbers between samples from g7 or WT
# SpatialFDR reports p-val corrected for multiple testing accounting for overlap between neigh

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# logFC   logCPM        F       PValue          FDR Nhood   SpatialFDR
# 557   2.727966 9.195126 45.98342 4.778562e-07 0.0004428134   557 0.0003661660
# 1329  2.608781 9.329211 48.14738 3.281017e-07 0.0004428134  1329 0.0003661660
# 2368  3.090191 9.101200 48.57267 3.051450e-07 0.0004428134  2368 0.0003661660
# 1292 -2.149224 9.410302 42.57036 8.864329e-07 0.0006058539  1292 0.0005121001
# 461  -3.234157 8.864387 39.65682 1.542167e-06 0.0006058539   461 0.0005479352
# 608   3.424775 8.549937 37.90402 2.179330e-06 0.0006058539   608 0.0005479352

################################################################################
## 4. Inspecting DA results
################################################################################

# Distribution of uncorrected P-values
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50) # Looks fine

# Volcano plot (each point is a neighbourhood)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)

# Visualize DA results relating them to the embedding of single cells
ebs_milo <- buildNhoodGraph(ebs_milo)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(ebs_milo, 
                          dimred = "UMAP", 
                          colour_by="Group", 
                          text_by = "SingleR_FINAL_EmbryoAtlas", 
                          text_size = 3, 
                          point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(ebs_milo, 
                                da_results, 
                                layout="UMAP",
                                alpha=0.1) 

pdf(file = "PAPER_RELATED/miloR_related/miloR_summary_UMAP_DA_results.pdf", height=7, width=15)
umap_pl + nh_graph_pl +
  patchwork::plot_layout(guides="collect")
dev.off()

# Assign cell types to neighbours
da_results <- annotateNhoods(ebs_milo, 
                             da_results, 
                             coldata_col = "SingleR_FINAL_EmbryoAtlas")

da_results %>%
  arrange(SpatialFDR) %>%
  head() 

# logFC   logCPM        F       PValue          FDR Nhood   SpatialFDR SingleR_FINAL_EmbryoAtlas SingleR_FINAL_EmbryoAtlas_fraction
# 557  -2.727966 9.195126 45.98342 4.778562e-07 0.0004428134   557 0.0003661660                 Allantois                          0.7634409
# 1329 -2.608781 9.329211 48.14738 3.281017e-07 0.0004428134  1329 0.0003661660                 Allantois                          0.6699029
# 2368 -3.090191 9.101200 48.57267 3.051450e-07 0.0004428134  2368 0.0003661660                 Allantois                          0.6162791
# 1292  2.149224 9.410302 42.57036 8.864329e-07 0.0006058539  1292 0.0005121001            Mixed mesoderm                          0.5523810
# 461   3.234157 8.864387 39.65682 1.542167e-06 0.0006058539   461 0.0005479352                  Epiblast                          0.4492754
# 608  -3.424775 8.549937 37.90402 2.179330e-06 0.0006058539   608 0.0005479352                 Allantois                          0.7818182

ggplot(da_results, aes(SingleR_FINAL_EmbryoAtlas_fraction)) + geom_histogram(bins=50)

# Let's exclude those neighbourhoods that are a mix of cells
da_results$SingleR_FINAL_EmbryoAtlas <- ifelse(da_results$SingleR_FINAL_EmbryoAtlas_fraction < 0.7, 
                                         "Mixed", da_results$SingleR_FINAL_EmbryoAtlas)

# Let's define the cell labels as factor to put them in a specific order

pdf(file = "PAPER_RELATED/miloR_related/miloR_summary_logFC_cell_types.pdf", height=7, width=8)
plotDAbeeswarm(da_results, group.by = "SingleR_FINAL_EmbryoAtlas",
               alpha = 0.1) + # Default for Spatial FDR. Colours in grey unsignificant neigh
  theme_bw() + 
  geom_point(size=0.2) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  annotate(geom="text", x=18,y=-2, label="Enriched in WT",
           color="brown3") +
  annotate(geom="text", x=18,y=2, label="Enriched in g7",
           color="slateblue")

dev.off()

##### Store results in RDS' files in case they have to be reloaded
# saveRDS(object = ebs_milo, file = "RDS/scRNAseq_EBs_milo_object.RDS")
# saveRDS(object = da_results, file = "RDS/scRNAseq_EBs_miloR_DA_results.RDS")

################################################################################
## 5. Compute the number of enriched neighboors per cell type (SpatialFDR < 0.1)
## Enriched in WT (logFC <0) or in g7 (logFC >0). This is a summary of da_results
################################################################################

neigh.sign_wt_enriched <- da_results %>% 
                            filter(SpatialFDR <0.1, logFC<0) %>%
                            arrange(SpatialFDR) 

neigh.sign_g7_enriched <- da_results %>% 
                            filter(SpatialFDR <0.1, logFC>0) %>%
                            arrange(SpatialFDR) 

# Put summary table all together
summary_da <- data.frame("Cell.Type" = unique(da_results$SingleR_FINAL_EmbryoAtlas))
summary_da <- merge.data.frame(x = summary_da, y= as.data.frame(table(neigh.sign_wt_enriched$SingleR_FINAL_EmbryoAtlas)),
                 by.x = "Cell.Type", by.y="Var1", all=TRUE)
colnames(summary_da)[2] <- "Nneigh.Enriched.WT"

summary_da <- merge.data.frame(x = summary_da, y= as.data.frame(table(neigh.sign_g7_enriched$SingleR_FINAL_EmbryoAtlas)),
                               by.x = "Cell.Type", by.y="Var1", all=TRUE)
colnames(summary_da)[3] <-"Nneigh.Enriched.g7"

################################################################################
## 6. AUTOMATIC GROUPING OF NEIGHBOURHOODS -> Selection of groups of interest to be tested 
## In this section, neighbourhood groups will be defined based on their enrichment
## in WT or g7. 
################################################################################

# ============
# Visualization of potential groups

# To automatically group neighbours into groups based on:
# (1) Number of shared cells between 2 neighbours,
# (2) Direction of logFC for DA
# (3) Difference in logFC

# Find groups
da_results <- groupNhoods(ebs_milo, 
                          da_results, 
                          overlap = 3, # Min overlap between adjacent neigh for merging (Default value:1)
                          da.fdr = 0.1, # Only aggregting this (default value: 0.1)
                          max.lfc.delta = 2) # Abs diff in lfc below which neighbourhoods should not be considered adjacent (NULL:default)
#head(da_results)

pdf(file = "PAPER_RELATED/miloR_related/Graph_milo_Complete_Nhoods_Groups_of_Interest.pdf", height=6, width=8)
plotNhoodGroups(ebs_milo, 
                da_results, 
                layout="UMAP",show_groups = c(11, # Pharyngeal mesoderm
                                              13, #Mixed mesoderm + nascent mesoderm
                                              14, # Cardiomyocytes
                                              15, # Erythroid
                                              16,  # Allantois
                                              4)) + #Primitive + Anterior Prim Streak
  scale_fill_manual(values=c("paleturquoise", 
                             "plum2", 
                             "orchid4", 
                             "red3", 
                             "slateblue2",
                             "burlywood2"), 
                    na.value= "ivory") # Rest of groups, not of interest
    
dev.off()


pdf(file = "PAPER_RELATED/miloR_related/miloR_summary_logFC_Nhood_Groups_ALL.pdf", height=7, width=8)

plotDAbeeswarm(da_results, "NhoodGroup") + 
  theme_bw() + 
  geom_point(size=0.2) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  annotate(geom="text", x=18,y=-2, label="Enriched in WT",
           color="brown3") +
  annotate(geom="text", x=18,y=2, label="Enriched in g7",
           color="slateblue")

dev.off()


# Only with the groups of interest
da_results_subset <- da_results[which(da_results$NhoodGroup %in% c(11,13,14,15,16,4)),] 


pdf(file = "PAPER_RELATED/miloR_related/miloR_summary_logFC_Nhood_Groups_HIGHLIGHTED.pdf", height=4, width=5)

plotDAbeeswarm(da_results_subset, "NhoodGroup") + 
  theme_bw() + 
  geom_point(size=0.2) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  annotate(geom="text", x=7,y=-2, label="Enriched in WT",
           color="brown3") +
  annotate(geom="text", x=7,y=2, label="Enriched in g7",
           color="slateblue")
dev.off()



################################################################################
## 6b.Store list of cells belonging to specific groups for further analysis, if
## required.
## With this list of cell barcodes we can then subset the complete seurat object
################################################################################

# FIRST SELECTION: All cells annotated to any of the neighbourhoods included in
# the groups of interest INDEPENDENTLY if these neighbourhoods are significantly
# enriched or not (in terms of abundance!)

# Define groups of interest
Groups = c("11", "13", "14", "15", "16" , "4")

# Identify neighbourhoods included in those groups
nhoods <- sapply(Groups, function(Gr){
  return(da_results$Nhood[which(da_results$NhoodGroup %in% Gr)])
})

cell_barcodes_NHoodGroups <- lapply(nhoods, function(nh)
  {
  
  interm <- sapply(nh, function(neigh)
    {
    return(colnames(ebs_milo)[ebs_milo@nhoods[,neigh]==1])
    })
  
  return(unique(unlist(interm)))
})


# A quick check about cells overlapping among these groups

pdf(file = "PAPER_RELATED/miloR_related/UPSET_Cell_Barcodes_NHood_Groups_of_interest.pdf", 
    width=5.5,  height=5) # Values in inches

upset(fromList(cell_barcodes_NHoodGroups),
      sets = names(cell_barcodes_NHoodGroups),
      keep.order=TRUE,
      mainbar.y.label = "Common Cell Barcodes",
      sets.x.label = "Cell Nhood Group size",
      point.size = 1.5,
      line.size = 0.5,
      text.scale = 1,
      set_size.show = TRUE,
      set_size.scale_max= 5500,
      nsets = length(cell_barcodes_NHoodGroups),
      sets.bar.color = c("paleturquoise", "plum2", "orchid4", "red3", "slateblue2","burlywood2"), ,
      main.bar.color = c("burlywood2", "paleturquoise", "plum2", "slateblue2", "red3", "orchid4", rep("gray",9)),
      nintersects = NA)

dev.off()

saveRDS(object = cell_barcodes_NHoodGroups, file = "RDS/scRNAseq_EBs_milo_CELLS_NhoodGroups_interest.RDS")

################################################################################
## FINAL. Store results of Diff Abundance Analysis
################################################################################

# Store this as an excel file
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "orchid4")

write.xlsx(x = list("Summary_Enriched_NHoods" = summary_da, "DAA" = da_results), 
           file= "PAPER_RELATED/miloR_related/scRNAseq_EBs_DiffAbundanceAnalysis_miloR.xlsx", 
           headerStyle=hs)
