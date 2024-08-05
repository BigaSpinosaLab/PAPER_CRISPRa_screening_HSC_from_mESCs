## The following functions 

PCA_custom <- function(Norm_Matrix,
                       TopGenes = NULL,
                       filename,
                       title,
                       ComponentX = "PC1",
                       ComponentY = "PC2",
                       colors,
                       colorby = NULL,
                       shapeby = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       res_path) 
{
  
  # Manual PCA computation with prcomp
  mat <- assay(Norm_Matrix) 
  
  # Remove the zero-variance columns: otherwise we cannot scale the data (optional)
  mat <- mat[which(apply(mat, 1, var) != 0),]
  
  # Set ntop value
  if (!is.null(TopGenes)) { 
    ntop = TopGenes 
  } else {
    ntop = nrow(mat)
  }
  
  # Optional: if you want to select the most variant genes
  rv = rowVars(mat)
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  pca = prcomp(t(mat[select,]), scale=TRUE, center=TRUE)
  PCs <- as.data.frame(pca$x)
  PCs$Sample <- rownames(PCs)
  
  if (!is.null(colorby)) { 
    PCs[[colorby]] <- as.character(colData(Norm_Matrix)[[colorby]]) 
  }
  
  if (!is.null(shapeby)) { 
    PCs[[shapeby]] <- as.character(colData(Norm_Matrix)[[shapeby]]) 
  }
  
  Variance <- round(summary(pca)$importance[2,]*100, digits=1)
  
  CompletePC <- PCs
  
  colnames(PCs)[colnames(PCs) == ComponentX] <- 'ComponentX'
  colnames(PCs)[colnames(PCs) == ComponentY] <- 'ComponentY'
  
  firstcomp <- as.numeric(substring(ComponentX, 3, 3))
  secondcomp <- as.numeric(substring(ComponentY, 3, 3))
  
  if (is.null(colorby)) {
    if (is.null(shapeby)) {
      plot1 <- ggplot(PCs, aes(ComponentX, 
                               ComponentY, 
                               label=Sample))
    } else {
      plot1 <- ggplot(PCs, aes(ComponentX, 
                               ComponentY,
                               shape=PCs[[shapeby]], 
                               label=Sample))
    }
    
  } else {
    if (is.null(shapeby)) {
      plot1 <- ggplot(PCs, aes(ComponentX, 
                               ComponentY,
                               color=PCs[[colorby]], 
                               label=Sample))
    } else {
      plot1 <- ggplot(PCs, aes(ComponentX, 
                               ComponentY,
                               color=PCs[[colorby]],
                               shape=PCs[[shapeby]], 
                               label=Sample))
    }
  }
  
  plot1 <- plot1 +
    geom_point(color="black",size=5,alpha=0.4) +
    geom_point(size=4,alpha=0.8) +
    scale_color_manual(values = colors) +
    geom_text_repel(size=3, color="black") +
    xlab(paste0(ComponentX, ": ", Variance[firstcomp], "% variance")) +
    ylab(paste0(ComponentY, ": ", Variance[secondcomp], "% variance")) +
    labs(color = colorby, shape = shapeby)+
    theme_bw()+
    theme(legend.title = element_text(size = 12,face="italic"),
          legend.text = element_text(size = 12),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(size=14, hjust = 1),
          axis.text.y = element_text(size=14, hjust = 1))
  
  if (is.null(xlim)) {
    xlim_vector = c(1.2*min(PCs$ComponentX), 1.2*max(PCs$ComponentX))
  } else {
    xlim_vector = xlim
  }
  
  if (is.null(ylim)) {
    ylim_vector = c(1.2*min(PCs$ComponentY), 1.2*max(PCs$ComponentY))
  } else {
    ylim_vector = ylim
  }
  
  plot1 <- plot1 +
    xlim(xlim_vector) +
    ylim(ylim_vector)  +
    coord_fixed() +
    ggtitle(title)
  
  ggsave(file = file.path(paste0(res_path, "/",filename,".pdf")), width=6, height=6,
         plot= last_plot(),
         device = pdf,
         path=res_path)  
  
  print(plot1)
  
  return(CompletePC)
}



#' This function is a wrap-up of getting results from DAA and generate_final_output
#' So, just one function calling is enough from the main script
#' Find parameters descriptions below

retrieve_final_results <- function(DESeqDataSet_res, 
                                   contrast.or.coef = NULL,
                                   test = NULL,
                                   comparison.name,
                                   lFC_shrinkage=TRUE,
                                   alpha = NULL,
                                   lFC = NULL,
                                   res_path,
                                   gene_annotations)
{
  
  results <- get_results_from_DEA(DESeqDataSet_res = DESeqDataSet_res, 
                                  contrast.or.coef = contrast.or.coef,
                                  test = test,
                                  comparison.name = comparison.name,
                                  lFC_shrinkage=lFC_shrinkage,
                                  alpha = alpha,
                                  lFC = lFC,
                                  res_path = res_path) 
  
  final.results <- generate_final_output(results.df = results, 
                                         gene_annotations = gene_annotations,
                                         comparison.name = comparison.name,
                                         alpha = alpha,
                                         lFC = lFC,
                                         res_path = res_path)
  
  return(final.results)
  
}


#' This function retrieves the results from DESseq2 given a specific contrast
#' 
#' @param DESeqDataSet_res DESeqDataSet object from runnig DESeq
#' @param contrast.or.coef The test you are specifying if a contrast or a coef name
#' Specify 'contrast' or 'coef' directly
#' @param test Comparison you want to retrieve. It could be a contrast or a coef name
#' as specified in previous parameter
#' @param comparison.name Name of the comparison under test
#' @param lFC_shrinkage Shall we conduct shrinkage? Not recommended for designs with interaction terms
#' @param alpha  Statistical criteria to be applied to adj pval. Default 5%
#' @param lFC log Fold Change criteria. Default 1 (abs value). If shrinkage applied
#' this criteria to be applied over the shrunken logFC
#' @param res_path Absolute path to store all tables and figures
#' @returns Data frame with DEA results


get_results_from_DEA <- function(DESeqDataSet_res, 
                                 contrast.or.coef = NULL,
                                 test,
                                 comparison.name,
                                 lFC_shrinkage=TRUE,
                                 alpha = 0.05,
                                 lFC = 1,
                                 res_path)
  
{
  # Show the comparison under evaluation
  print(paste0("The comparison to be made is: ",comparison.name))
  
  # Sanity check for the contrast equivalent coefficient
  if(lFC_shrinkage & (contrast.or.coef == "contrast"))
  {
    print("The equivalent model coefficient for your contrast will be inferred. Advisable to check if this inference is correct. In case of doubt, disable shrinkage")
    
  }
  
  if(!lFC_shrinkage){
    print("No log FC shrinkage will be conducted")
  }
  
  # a) Get the results from a previous run of DESeq2
  if(contrast.or.coef=="contrast")
  {
    results <- results(object = DESeqDataSet_res,
                       contrast= test,
                       alpha=alpha)
    
  }else if (contrast.or.coef=="coef")
  {
    results <- results(object = DESeqDataSet_res,
                       name= test,
                       alpha=alpha)
    
  }
  
  # Print a summary of results
  print("Results summary from DESeq2 analysis (no logFC criteria applied)")
  print(summary(results))
  
  # b) LogFC shrinkgae computation and aggregate results
  if(lFC_shrinkage)
  {
    if(contrast.or.coef=="contrast")
    {
      # Infer the coef from the contrast
      values <- unlist(strsplit(results@elementMetadata$description[2], split=": "))[2]
      values <- as.logical(as.numeric(unlist(strsplit(values, split=","))))
      
      # Conduct shrinkage 
      res.shr <- lfcShrink(dds = DESeqDataSet_res,
                           coef= resultsNames(DESeqDataSet_res)[values],
                           type = "apeglm")
      
    }else if(contrast.or.coef=="coef")
    {
      # Coef number for the test
      coef.number = which(resultsNames(DESeqDataSet_res) %in% test)
      
      # Conduct shrinkage 
      res.shr <- lfcShrink(dds = DESeqDataSet_res,
                           coef = coef.number,
                           type = "apeglm")
    }
    
    
    results.df <- cbind(as.data.frame(results),
                        "Shrunken_lFC" = as.data.frame(res.shr)$log2FoldChange,
                        "Shrunken_lFCSE" = as.data.frame(res.shr)$lfcSE)
  }
  
  results.df <- results.df[order(results.df$padj),]
  
  # Plot MAplot before and after the logFC shrinkage (if applicable)
  pdf(file = file.path(res_path,paste0("figures/MAplot_",comparison.name,".pdf")), width=16, height=10)
  if(!lFC_shrinkage)
  {
    plotMA(results,ylim=c(-3,3),main=paste0("TEST: ", comparison.name))
    
  }else{
    par(mfrow=c(1,2))    
    plotMA(results,ylim=c(-3,3),main=paste0("TEST: ", comparison.name))
    plotMA(res.shr, ylim=c(-3,3), main=paste0("TEST: ", comparison.name," - Shrinkage lFC"))
  }
  dev.off()
  
  # c) Print number of DEGs if applying abs(log2FoldChange) >1 - shrunken logFC
  print("The number of DEGs is (adj pval and log2FC criteria):")
  if(!lFC_shrinkage)
  {
    print(paste("UP:", nrow(results.df %>%
                              filter(padj < alpha, log2FoldChange > lFC)), sep=""))
    print(paste("DOWN:", nrow(results.df %>%
                                filter(padj < alpha, log2FoldChange < -(lFC))), sep=""))
  }else{
    print(paste("UP:", nrow(results.df %>%
                              filter(padj < alpha, Shrunken_lFC > lFC)), sep=""))
    print(paste("DOWN:", nrow(results.df %>%
                                filter(padj < alpha, Shrunken_lFC < -(lFC))), sep=""))
  }
  
  # Return results for further analysis
  return(results.df)
}


#' This function annotates results from DEa and generates proper outputs: excel, bed and figures
#' 
#' @param results_df DESeqDataSet object from runnig DESeq
#' @returns Data frame with DAA results


generate_final_output <- function(results.df, 
                                  gene_annotations,
                                  comparison.name,
                                  alpha = 0.05,
                                  lFC = 1,
                                  res_path)
  
{
  # Include gene annotations in final results
  # Match the order of annotations and results
  if(!is.null(gene_annotations))
  {
    gene_annotations <- gene_annotations[match(rownames(results.df), gene_annotations$gene_id),]
    results.df <- cbind(results.df, gene_annotations[,c("gene_name", "Entrez", "gene_id", "chrom", "start", "end")])
  }else{
    results.df$gene_name <- rownames(results.df)
  }
  
  # Store results in an excel highlighting DEGs
  formatting_write_results_dea(final_results = results.df,
                           adjpval_criteria = alpha,
                           excel_path = res_path,
                           comparison.name = comparison.name)
  
  # Check if there is log FC shrinkage or not to consider it for Volcanos
  if(length(which(colnames(results.df) %in% "Shrunken_lFC"))>0){shrink=TRUE}else{shrink=FALSE}
  
  # Generate Volcano Plot
  #===============================================
  if(shrink)
  {
    up.n = nrow(results.df %>% filter(padj < alpha, Shrunken_lFC>0))
    down.n = nrow(results.df %>% filter(padj < alpha, Shrunken_lFC<0))
  }else{
    up.n = nrow(results.df %>% filter(padj < alpha, log2FoldChange>0))
    down.n = nrow(results.df %>% filter(padj < alpha, log2FoldChange<0))
    }


  p <- EnhancedVolcano(toptable = results.df,
                       lab= results.df$gene_name,
                       x =ifelse(shrink, 'Shrunken_lFC','log2FoldChange'),
                       y= 'padj',
                       xlab= ifelse(shrink, 'Shrunken log2 FC','Log2 Fold change'),
                       ylab = "-Log10(adj pval)",
                       title = comparison.name,
                       subtitle= paste0("RNAseq data analysis (", nrow(results.df),") genes"),
                       #xlim=c(-5,5),
                       pCutoff = alpha,
                       FCcutoff = lFC,
                       caption = paste0("FC cutoff: ", lFC,"; adj p-val cutoff: ", alpha),
                       pointSize = 1.5,
                       labSize = 3,
                       max.overlaps = 25,
                       drawConnectors = TRUE,
                       legendLabels=c("Not sign.", "log2FC", "adj pval","adj pval & log2FC"),
                       legendPosition = "right",
                       legendIconSize = 5.0,
                       legendLabSize = 10)
  
  grob <- grobTree(textGrob(paste0("UP-REG (adj pval) N=", up.n), x=0.7,  y=0.5, hjust=0,
                            gp=gpar(col="gray43", 
                                    fontsize=13, 
                                    fontface="bold")))
  
  grob2 <- grobTree(textGrob(paste0("DOWN-REG (adj pval) N=", down.n), x=0.1,  y=0.5, hjust=0,
                             gp=gpar(col="gray43", 
                                     fontsize=13, 
                                     fontface="bold")))
  
  p <- p + annotation_custom(grob) + annotation_custom(grob2) 
  
  
  pdf(paste0(res_path, "/figures/VOLCANO_", comparison.name,".pdf"), width=16, height=10)
  print(p)
  dev.off()
  
  # Return the final results
  return(results.df)
  
}  




formatting_write_results_dea <- function(final_results,
                                     adjpval_criteria,
                                     excel_path,
                                     comparison.name)
{
  # Prepare workbook
  wb <- createWorkbook()
  addWorksheet(wb, comparison.name)
  writeData(wb, comparison.name, final_results)
  
  # Prepare Conditional Formatting
  highlight <- createStyle(fontColour = "#000000",
                           bgFill = "lightblue") #Colour rows that meet adj pval criteria
  
  conditionalFormatting(wb,
                        comparison.name,
                        cols = 1:length(final_results),
                        rows = 1:nrow(final_results),
                        rule = paste0("$F1<",adjpval_criteria),  # pval column is column F
                        style = highlight)
  
  
  ## create and add a style to the column headers
  headerStyle1 <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "orchid3")
  
  addStyle(wb, 
           sheet = 1, headerStyle1, 
           rows = 1, 
           cols=1:ncol(final_results),  # There are 11 columns in total
           gridExpand = TRUE)
  
  setColWidths(wb, sheet = 1, cols = 1:ncol(final_results), widths = "auto")
  
  # Store results in an excel file
  saveWorkbook(wb, 
               file=paste0(res_path, "/tables/", comparison.name,"_Diff_Exprs_Analysis.xlsx"), 
               overwrite = TRUE)
}


# This function has been adapted from the corresponding function from clusterProfiler package

simplify_custom <- function(res, 
                            cutoff, # Similarity cutoff
                            by, 
                            select_fun=min,  # Minimum pval 
                            measure="Wang", # This is the default method 
                            ontology ="BP",
                            semData=NULL) 
{
  
  if (measure == "Wang") {
    semData <- GOSemSim::godata(ont = ontology)
  } else {
    stop("godata should be provided for IC-based methods...")
  }
  
  sim <- GOSemSim::mgoSim(res$ID, res$ID,
                          semData = semData,
                          measure=measure,
                          combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- tidyr::gather(sim.df, go2, similarity, -go1)
  
  sim.df <- sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)
  
  ID <- res$ID
  
  GO_to_remove <- character()
  
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset <- sim.df[ii,]
    
    jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  
  return(res[!res$ID %in% GO_to_remove, ])
} # end function

