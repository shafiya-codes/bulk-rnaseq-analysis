# 🧬 Bulk RNA-seq Analysis Pipeline (R)

**Author:** Shafiya Sakina  
**Focus:** RNA-seq | Differential Gene Expression | Transcriptomics  

---

## 📌 Overview

This repository contains an **end-to-end Bulk RNA-seq analysis workflow** implemented in **R**, covering all major steps from raw sequencing data to biological interpretation.  
The pipeline follows best practices used in real-world bioinformatics and research environments.

This project demonstrates practical understanding of:
- RNA-seq data processing
- Statistical modeling of gene expression
- Batch effect correction
- Visualization of high-dimensional biological data

---

## 🧪 Datasets Used

### 1. Public RNA-seq FASTQ files
- Source: ENA / SRA  
- Organism: *Saccharomyces cerevisiae*  
- Paired-end sequencing data

### 2. Cancer gene expression dataset
- Source: GEO (HTSeq counts)  
- Biological comparison: **Tumor Stage T1 vs T2**

⚠️ Raw sequencing data is not stored in this repository.  
Scripts automatically download data from public repositories.

---

## 🔬 Analysis Workflow

1. Quality control of raw FASTQ files (FastQC)
2. Reference genome and annotation preparation
3. Read alignment using Rsubread
4. Gene-level quantification using featureCounts
5. Exploratory data analysis (RLE, PCA)
6. Normalization and batch correction using RUVSeq
7. Differential expression analysis using DESeq2
8. Visualization (volcano plots, heatmaps)

---

## 🛠 Tools & Packages

- Rsubread  
- DESeq2  
- RUVSeq  
- fastqcr  
- ggplot2  
- pheatmap  
- EnhancedVolcano  

---

## 📁 Repository Structure


---

## 🎯 Learning Outcomes

Through this project, I gained hands-on experience in:
- Handling raw RNA-seq data
- Implementing statistically robust RNA-seq pipelines
- Correcting technical noise and batch effects
- Interpreting gene expression changes in a biological context

---

## 🚀 Future Improvements

- Adapter trimming and read filtering
- Functional enrichment analysis (GO / KEGG)
- Workflow automation using Snakemake / Nextflow
- Containerization with Docker

---

## 📬 Contact

**Shafiya Sakina**  
GitHub: https://github.com/shafiya-codes  
LinkedIn: (add your LinkedIn link here)
