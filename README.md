# PAPER_CRISPRa_screening_HSC_from_mESCs
Collection of scripts used to perform the analysis of Galán-Palma L. et al. 

This repository includes scripts required for the scRNAseq datasets analysis included in Galán-Palme et al. 
'In vivo genome-wide CRISPRa screen identifies novel genes that generate hematopoietic stem and progenitor cells from mouse embryonic stem cells' 

All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

### scRNAseq 10x Embryoid Bodies (EBs)

Scripts required to reproduce the complete scRNAseq EBs dataset (3x SADEiGEN and 3x WT samples) analysis, specifically:

- Data preprocessing: scripts from 0 to 4.
    - Script 0 refers to CellRanger execution to obtain raw features, barcodes and expression matrices. 
    - Script 1 refers to exploration (basic statistics) of Cell Ranger output.
    - Script 2 refers to doublets identification.
    - Script 3 refers to generation of main Seurat object with raw data from all samples. This object is used for downstream analysis.
    - Script 4 refers to cells and genes filtering (Quality Control).

- Downstream analysis: scripts from 5 to 8.
    - Script 5 refers to dimensionality reduction and clustering process.
    - Script 6 refers to cell type annotation with Single R and based on the EmbryoMouseAtlas from Pijuan-Sala et al. 2019 Nature. The same script includes the retrieval of the reference dataset.
    - Script 7 refers to the differential abundance analysis with miloR for testing SADEiGEN vs WT condition.
    - Script 8 refers to the pseudobulk differential expression analysis for the subset of Kdr positive cells.

Raw FASTQ files and CellRanger output are available through GEO. Additionally, an RDS file including a Seurat object is also available in the same accession number. This RDS file includes up to Script 6 execution.

NOTE: Corresponding HTML reports generated from scripts 1, 2, 4 and 5 are also available in a subfolder.

### scRNAseq 10x Bone Marrow (BM)

Scripts required to reproduce the complete scRNAseq BM dataset (1x SADEiGEN and 1x WT samples) analysis, specifically:

- Data preprocessing: scripts from 0 to 4.

    - Script 0 refers to CellRanger execution to obtain raw features, barcodes and expression matrices.
    - Script 1 refers to exploration (basic statistics) of Cell Ranger output.
    - Script 2 refers to doublets identification.
    - Script 3 refers to generation of main Seurat object with raw data from all samples. This object is used for downstream analysis.
    - Script 4 refers to cells and genes filtering (Quality Control).

- Downstream analysis: scripts from 5 to 8.

    - Script 5 refers to the retrieval of the reference dataset.
    - Script 6 refers to cell type annotation with Single R and based on the Mouse Cell Atlas from Wang et al. 2023. 
    - Script 7 refers to dimensionality reduction, data integration and clustering process.
    - Script 8 refers to the results visualization.

Raw FASTQ files and CellRanger output are available through GEO. Additionally, an RDS file including a Seurat object is also available in the same accession number. This RDS file includes up to Script 7 execution.

NOTE: Corresponding HTML reports generated from scripts 1, 2, 4 and 8 are also available in a subfolder.
