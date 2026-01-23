# Bulk RNA-seq Analysis (R)

This repository contains a reproducible bulk RNA-seq analysis pipeline in R.

## Workflow
- Read alignment using Rsubread
- Gene quantification with featureCounts
- Differential expression analysis using DESeq2
- Visualization (PCA, volcano plot, heatmap)

## Requirements
- R (>= 4.2)
- Rsubread, DESeq2, ggplot2, pheatmap, EnhancedVolcano

## Usage
1. Place paired-end FASTQ files in `data/`
2. Place reference genome and GTF in `reference/`
3. Run the script in `scripts/`

