############################################################
# Bulk RNA-seq Analysis Pipeline
# Author: Shafiya Sakina
#
# Description:
# End-to-end bulk RNA-seq workflow in R including:
# - Quality control
# - Alignment
# - Quantification
# - Normalization & batch correction
# - Differential expression analysis
# - Visualization
############################################################

############################
# 1. Load Required Packages
############################
suppressPackageStartupMessages({
  library(Rsubread)
  library(data.table)
  library(RUVSeq)
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(fastqcr)
  library(R.utils)
})

############################
# 2. Define Project Paths
############################
project_dir <- getwd()

data_dir    <- file.path(project_dir, "data")
ref_dir     <- file.path(project_dir, "reference")
results_dir <- file.path(project_dir, "results")
qc_dir      <- file.path(results_dir, "fastqc")

dir.create(data_dir, showWarnings = FALSE)
dir.create(ref_dir, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)
dir.create(qc_dir, showWarnings = FALSE)

############################
# 3. Download RNA-seq FASTQ Files
############################
fastq_urls <- c(
  "ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_1.fastq.gz",
  "ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/006/SRR5924196/SRR5924196_2.fastq.gz",
  "ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_1.fastq.gz",
  "ftp.sra.ebi.ac.uk/vol1/fastq/SRR592/008/SRR5924198/SRR5924198_2.fastq.gz"
)

for (url in fastq_urls) {
  dest <- file.path(data_dir, basename(url))
  if (!file.exists(dest)) {
    download.file(url, dest)
  }
}

############################
# 4. Download Reference Genome & Annotation
############################
genome_url <- "ftp://ftp.ensembl.org/pub/release-96/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz"
gtf_url    <- "ftp://ftp.ensembl.org/pub/release-96/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.96.gtf.gz"

genome_gz <- file.path(ref_dir, basename(genome_url))
gtf_gz    <- file.path(ref_dir, basename(gtf_url))

download.file(genome_url, genome_gz)
download.file(gtf_url, gtf_gz)

gunzip(genome_gz, overwrite = TRUE)
gunzip(gtf_gz, overwrite = TRUE)

genome_fa <- sub(".gz", "", genome_gz)
gtf_file  <- sub(".gz", "", gtf_gz)

############################
# 5. Quality Control (FastQC)
############################
fastqc_install()
fastqc(dir = data_dir, outdir = qc_dir)
qc_report <- qc_aggregate(qc_dir)
print(qc_report)

############################
# 6. Build Genome Index
############################
index_base <- file.path(ref_dir, "yeast_index")

if (!file.exists(paste0(index_base, ".00.b.array"))) {
  buildindex(
    basename = index_base,
    reference = genome_fa
  )
}

############################
# 7. Align Reads
############################
reads1 <- list.files(data_dir, pattern = "_1.fastq.gz$", full.names = TRUE)
reads2 <- list.files(data_dir, pattern = "_2.fastq.gz$", full.names = TRUE)

stopifnot(length(reads1) == length(reads2))

align(
  index = index_base,
  readfile1 = reads1,
  readfile2 = reads2,
  input_format = "gzFASTQ",
  output_format = "BAM",
  nthreads = 8
)

bam_files <- list.files(data_dir, pattern = ".BAM$", full.names = TRUE)

############################
# 8. Mapping Statistics
############################
mapping_stats <- propmapped(bam_files)
print(mapping_stats)

############################
# 9. Gene Quantification
############################
fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE
)

count_matrix <- fc$counts
write.csv(count_matrix,
          file = file.path(results_dir, "gene_counts.csv"))

############################################################
# 10. Differential Expression (Cancer Dataset)
############################################################
counts <- read.table(
  "GSE143630_RCC_htseq_counts.txt",
  header = TRUE,
  row.names = 1
)

# Filter low-count genes
keep <- apply(counts, 1, function(x) sum(x > 0) >= 2)
counts <- counts[keep, ]

# Define conditions
t1_n <- length(grep("^T1", colnames(counts)))
t2_n <- length(grep("^T2", colnames(counts)))

condition <- factor(c(rep("T1", t1_n), rep("T2", t2_n)))

set <- newSeqExpressionSet(
  as.matrix(counts),
  phenoData = data.frame(condition,
                         row.names = colnames(counts))
)

############################
# 11. Exploratory Analysis
############################
colors <- brewer.pal(3, "Set2")
plotRLE(set, col = colors[condition], outline = FALSE)
plotPCA(set, col = colors[condition])

############################
# 12. Normalization & Batch Correction
############################
set <- betweenLaneNormalization(set, which = "upper")
differences <- makeGroups(condition)
set_ruv <- RUVs(set, rownames(counts), k = 1, differences)

############################
# 13. DESeq2 Analysis
############################
dds <- DESeqDataSetFromMatrix(
  countData = counts(set_ruv),
  colData   = pData(set_ruv),
  design    = ~ W_1 + condition
)

dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res),
          file = file.path(results_dir, "DEG_results.csv"))

############################
# 14. Visualization
############################
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Differential Expression: T1 vs T2"
)

top_genes <- head(rownames(res[order(res$padj), ]), 20)
pheatmap(counts[top_genes, ],
         cluster_cols = TRUE,
         scale = "row")

############################################################
# End of Script
############################################################
