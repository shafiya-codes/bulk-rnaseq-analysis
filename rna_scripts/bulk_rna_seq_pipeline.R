############################################################
# Bulk RNA-seq Analysis Pipeline (Reproducible)
# Author: Shafiya Sakina
#
# Description:
# End-to-end bulk RNA-seq workflow using paired-end data:
#  - Alignment (Rsubread)
#  - Gene quantification (featureCounts)
#  - Differential expression analysis (DESeq2)
#  - Visualization (PCA, Volcano, Heatmap)
#
# Organism: Saccharomyces cerevisiae (yeast)
############################################################

############################
# 0. Global Settings
############################
set.seed(123)

############################
# 1. Load Required Packages
############################
suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(RColorBrewer)
})

############################
# 2. Project Structure
############################
project_dir <- getwd()

dirs <- list(
  data      = file.path(project_dir, "data"),
  reference = file.path(project_dir, "reference"),
  results   = file.path(project_dir, "results")
)

lapply(dirs, dir.create, showWarnings = FALSE)

############################
# 3. Input FASTQ Files
############################
reads1 <- list.files(dirs$data, pattern = "_1.fastq.gz$", full.names = TRUE)
reads2 <- list.files(dirs$data, pattern = "_2.fastq.gz$", full.names = TRUE)

if (length(reads1) == 0 || length(reads2) == 0) {
  stop("❌ FASTQ files not found. Please place *_1.fastq.gz and *_2.fastq.gz in /data")
}

stopifnot(length(reads1) == length(reads2))

############################
# 4. Reference Genome & Annotation
############################
genome_fa <- file.path(
  dirs$reference,
  "Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa"
)

gtf_file <- file.path(
  dirs$reference,
  "Saccharomyces_cerevisiae.R64-1-1.96.gtf"
)

if (!file.exists(genome_fa) || !file.exists(gtf_file)) {
  stop("❌ Reference genome or GTF file missing in /reference")
}

############################
# 5. Build Genome Index
############################
index_base <- file.path(dirs$reference, "yeast_index")

if (!file.exists(paste0(index_base, ".00.b.array"))) {
  message("🔧 Building genome index...")
  buildindex(
    basename = index_base,
    reference = genome_fa
  )
}

############################
# 6. Alignment
############################
message("🚀 Aligning reads...")

align(
  index = index_base,
  readfile1 = reads1,
  readfile2 = reads2,
  input_format = "gzFASTQ",
  output_format = "BAM",
  nthreads = 6
)

bam_files <- list.files(dirs$data, pattern = "\\.BAM$", full.names = TRUE)

if (length(bam_files) == 0) {
  stop("❌ Alignment failed: BAM files not generated")
}

############################
# 7. Gene Quantification
############################
message("🧮 Quantifying genes...")

fc <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE
)

count_matrix <- fc$counts

write.csv(
  count_matrix,
  file = file.path(dirs$results, "gene_counts.csv")
)

############################
# 8. Sample Metadata
############################
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c(rep("Control", length.out = ncol(count_matrix) / 2),
      rep("Treated", length.out = ncol(count_matrix) / 2))
  )
)

############################
# 9. Differential Expression (DESeq2)
############################
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ condition
)

# Filter low-expression genes
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]

write.csv(
  as.data.frame(res),
  file = file.path(dirs$results, "DESeq2_results.csv")
)

############################
# 10. Visualization
############################
vsd <- vst(dds, blind = FALSE)

# PCA
plotPCA(vsd, intgroup = "condition")

# Volcano
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  title = "Differential Expression: Treated vs Control"
)

# Heatmap (Top 20 genes)
top_genes <- head(rownames(res), 20)

pheatmap(
  assay(vsd)[top_genes, ],
  scale = "row",
  annotation_col = sample_info,
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
)

############################
# 11. Reproducibility Info
############################
writeLines(
  capture.output(sessionInfo()),
  con = file.path(dirs$results, "sessionInfo.txt")
)

############################################################
# End of Script
############################################################

