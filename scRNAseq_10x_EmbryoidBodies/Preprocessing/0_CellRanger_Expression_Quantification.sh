#!/bin/bash

#SBATCH --job-name=CellRanger
#SBATCH --partition=fast
#SBATCH --cpus-per-task=6
#SBATCH --nodes=1 
#SBATCH --mem MaxMemPerNode
#SBATCH -t 40:00:00 
#SBATCH -o logs/CellRanger.out
#SBATCH -e logs/CellRanger.err

# REMARK 1: This script is not considering Feature Barcode analysis. Code should be
# adapted if this is required

# REMARK 2: This script is considering that one same library has been sequenced in
# different flow cells

# REMARK 3: No array execution can be done since there is only one 'fast' node

#=========================
# User defined parameters: relevant paths
#=========================
# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory
WKD=$ROOTDIR'/CRISPR_screening_mESCs_ABigas/scRNAseq_EBs_7g_vs_WT'

#=========================
# General configuration
#=========================
START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

################################################################################
##       Obtain raw gene expression information from scRNAseq 10X Genomics data
##                          using Cell Ranger software 
################################################################################

# Check the several tutorials for CellRanger in 10X Genomics website
# https://www.10xgenomics.com/support/software/cell-ranger/latest

###########################
## 1. Other relevant paths for I/O
###########################

# Specify path to raw data
RAWDATA=$WKD'/raw_data'

# Specify the samplesheet with info about sample name, flow cells and fastq dir
SAMPLESHEET=$RAWDATA'/Samplesheet_x_CellRanger.txt'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify singularity image to be used
CELLRANGER="cellranger_v7.2.0.sif"

# Specify the path to the transcriptome to be used
REF=$DB_PATH'/CellRanger/refdata-gex-mm10-2020-A'

#################################################
## 3. Run CellRanger
#################################################

cd $WKD'/cellranger_out'

while IFS=";" read -r sample_id flowcells fastqs; 
do
  singularity exec $IMAGES_PATH/$CELLRANGER cellranger count \
                --include-introns true \
                --id=$sample_id \
                --sample=$flowcells \
                --fastqs=$RAWDATA/$fastqs \
                --transcriptome=$REF

done < $SAMPLESHEET

#######################
## 4. End Preprocessing
#######################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Processing Time: $DIFF seconds"
