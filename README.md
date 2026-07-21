# 🧬 Bulk RNA-seq Analysis Pipeline (R)

![R](https://img.shields.io/badge/R-%3E%3D4.2-276DC3?style=flat&logo=r&logoColor=white)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.17-brightgreen?style=flat)
![License](https://img.shields.io/badge/license-MIT-blue?style=flat)
![Status](https://img.shields.io/badge/status-active-success?style=flat)

A fully reproducible, end-to-end bulk RNA-seq analysis pipeline in R — from raw FASTQ files to publication-ready differential expression results. Covers QC, alignment, normalization, batch correction, DEG analysis, and rich visualizations.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Pipeline Workflow](#-pipeline-workflow)
- [Repository Structure](#-repository-structure)
- [Requirements](#-requirements)
- [Quick Start](#-quick-start)
- [Output & Visualizations](#-output--visualizations)
- [Dataset](#-dataset)
- [Key Parameters](#-key-parameters)
- [Citation](#-citation)

---

## 🔬 Overview

This pipeline implements a standard bulk RNA-seq workflow for differential expression analysis. It is designed with reproducibility and clarity in mind — each step is modular, well-commented, and produces interpretable outputs at every stage.

**Key features:**
- Paired-end FASTQ alignment with **Rsubread / STAR**
- Gene-level quantification with **featureCounts**
- Normalization and DEG analysis with **DESeq2**
- Batch effect detection and correction with **sva / limma**
- Publication-ready plots: PCA, volcano, MA, heatmap

---

## 🔄 Pipeline Workflow

```
Raw FASTQ files
      │
      ▼
┌─────────────────────┐
│  QC (FastQC/MultiQC)│  ← Per-sample quality metrics
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Alignment          │  ← Rsubread / STAR → BAM files
│  (genome index)     │
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  featureCounts      │  ← Gene-level read counts
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Normalization      │  ← DESeq2 size factors / vst
│  + Batch Correction │  ← ComBat-seq / removeBatchEffect
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Differential       │  ← DESeq2 Wald test
│  Expression (DESeq2)│  ← LFC shrinkage (apeglm)
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Visualization      │  ← PCA, Volcano, Heatmap, MA
└─────────────────────┘
      │
      ▼
┌─────────────────────┐
│  Enrichment         │  ← GO / KEGG (clusterProfiler)
└─────────────────────┘
```

---

## 📁 Repository Structure

```
bulk-rnaseq-analysis/
├── scripts/
│   ├── 01_qc.R                  # FastQC wrapper, MultiQC summary
│   ├── 02_alignment.R           # Rsubread index build + align
│   ├── 03_quantification.R      # featureCounts gene-level counts
│   ├── 04_normalization.R       # DESeq2 normalization, VST
│   ├── 05_batch_correction.R    # ComBat-seq batch correction
│   ├── 06_deseq2_analysis.R     # Differential expression
│   └── 07_visualization.R       # PCA, volcano, heatmap, MA plot
├── data/                        # Place raw FASTQ files here (gitignored)
├── reference/                   # Reference genome + GTF (gitignored)
├── results/                     # Output tables and plots (gitignored)
├── .gitignore
├── LICENSE
└── README.md
```

---

## ⚙️ Requirements

**R version:** ≥ 4.2

Install all dependencies in one go:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Rsubread",
  "DESeq2",
  "EnhancedVolcano",
  "pheatmap",
  "clusterProfiler",
  "org.Hs.eg.db",
  "sva",
  "limma",
  "apeglm"
))

install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "ggrepel"))
```

**System tools (optional but recommended):**
- FastQC ≥ 0.11
- MultiQC ≥ 1.14
- Trim Galore ≥ 0.6 (for adapter trimming)

---

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/shafiya-codes/bulk-rnaseq-analysis.git
cd bulk-rnaseq-analysis
```

### 2. Set up your data

```
data/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz

reference/
├── hg38.fa          # Reference genome (GRCh38)
└── hg38.gtf         # Gene annotation (Ensembl)
```

### 3. Run the pipeline step by step

```r
# In R or RStudio:
source("scripts/01_qc.R")
source("scripts/02_alignment.R")
source("scripts/03_quantification.R")
source("scripts/04_normalization.R")
source("scripts/05_batch_correction.R")
source("scripts/06_deseq2_analysis.R")
source("scripts/07_visualization.R")
```

### 4. Check results

Outputs are written to `results/`:
- `results/counts/` — raw and normalized count matrices
- `results/DEGs/` — DESeq2 results tables (CSV)
- `results/plots/` — all figures (PNG/PDF)

---

## 📊 Output & Visualizations

The pipeline produces the following plots in `results/plots/`:

| Plot | Description |
|------|-------------|
| `pca_plot.png` | PCA of VST-normalized samples, colored by condition/batch |
| `volcano_plot.png` | EnhancedVolcano — log2FC vs −log10(padj), top DEGs labeled |
| `heatmap_top50.png` | Heatmap of top 50 DEGs (pheatmap, z-score scaled) |
| `ma_plot.png` | MA plot with LFC shrinkage (apeglm) |
| `dispersion_plot.png` | DESeq2 dispersion estimates |
| `go_dotplot.png` | GO biological process enrichment (clusterProfiler) |

> **Note:** To add your own output screenshots, place PNG files in `results/plots/` and they will render here automatically when you update this README.

---

## 📦 Dataset

This pipeline was developed and tested using publicly available data from the **NCBI GEO** database.

**Recommended test dataset:**
- **GEO accession:** [GSE96870](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96870) — Mouse RNA-seq, airway smooth muscle cells
- **Why this dataset:** Small size (~500MB), well-annotated, commonly used for DESeq2 benchmarking

To download:
```bash
# Using SRA Toolkit
prefetch SRR5442248 SRR5442249 SRR5442250
fastq-dump --split-files --gzip SRR5442248
```

Or use the Bioconductor shortcut for a count matrix directly:
```r
library(airway)
data(airway)  # Pre-built SummarizedExperiment — great for quick testing
```

---

## 🔧 Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `padj_cutoff` | 0.05 | Adjusted p-value threshold for DEGs |
| `lfc_cutoff` | 1.0 | Log2 fold-change cutoff |
| `min_count` | 10 | Minimum read count filter (across all samples) |
| `lfc_shrinkage` | `apeglm` | LFC shrinkage estimator for volcano/MA plots |
| `batch_method` | `ComBat-seq` | Batch correction method (counts-aware) |

---

## 📚 Citation

If you use this pipeline in your work, please cite the core tools:

- **DESeq2:** Love MI, Huber W, Anders S. *Genome Biology* 2014. doi:10.1186/s13059-014-0550-8
- **Rsubread / featureCounts:** Liao Y et al. *Nucleic Acids Research* 2019. doi:10.1093/nar/gkz114
- **clusterProfiler:** Wu T et al. *The Innovation* 2021. doi:10.1016/j.xinn.2021.100141
- **EnhancedVolcano:** Blighe K et al. Bioconductor 2018.

---

## 👩‍💻 Author

**SHAFIYA SAKINA**
Bioinformatician | Computational Genomics | RNA-seq · Multi-omics · ML for Precision Biology

[![ORCID](https://img.shields.io/badge/ORCID-0009--0007--8341--5285-A6CE39?style=flat&logo=orcid)](https://orcid.org/0009-0007-8341-5285)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-shafiya--sakina-0077B5?style=flat&logo=linkedin)](https://linkedin.com/in/shafiya-sakina)
[![GitHub](https://img.shields.io/badge/GitHub-shafiya--codes-181717?style=flat&logo=github)](https://github.com/shafiya-codes)

---

*MIT License · Feel free to fork, adapt, and build on this pipeline.*
