############################################################
# Bulk RNA-seq Pipeline | Step 03: Gene Quantification
# Author: Shafiya Sakina
# Dataset: GSE96870 — or airway (Bioconductor demo)
# Description: Gene-level read counting with featureCounts.
#              Falls back to airway count matrix for demo.
############################################################

suppressPackageStartupMessages({
  library(Rsubread)
  library(airway)
})

# ----- Paths -----
project_dir  <- getwd()
ref_dir      <- file.path(project_dir, "reference")
bam_dir      <- file.path(project_dir, "results", "bam")
counts_dir   <- file.path(project_dir, "results", "counts")

dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)

gtf_file  <- file.path(ref_dir, "hg38.gtf")
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# ----- featureCounts (if BAM files exist) -----
if (length(bam_files) > 0 && file.exists(gtf_file)) {

  message("Running featureCounts on ", length(bam_files), " BAM files...")

  fc <- featureCounts(
    files              = bam_files,
    annot.ext          = gtf_file,
    GTF.featureType    = "exon",
    GTF.attrType       = "gene_id",
    isGTFAnnotationFile = TRUE,
    isPairedEnd        = TRUE,
    nthreads           = 8,
    countMultiMappingReads = FALSE   # conservative: exclude multi-mappers
  )

  count_matrix <- fc$counts

  # Clean column names (remove .bam suffix)
  colnames(count_matrix) <- sub("\\.bam$", "", basename(colnames(count_matrix)))

  message("Count matrix: ", nrow(count_matrix), " genes x ",
          ncol(count_matrix), " samples")

  # Save raw counts
  write.csv(count_matrix,
            file = file.path(counts_dir, "raw_counts.csv"))

  # Save featureCounts stats
  write.csv(fc$stat,
            file = file.path(counts_dir, "featurecounts_stats.csv"))

# ----- Fallback: airway Bioconductor demo dataset -----
} else {
  message("BAM files or GTF not found. Loading airway demo dataset...")

  data(airway)

  count_matrix <- assay(airway, "counts")
  col_data     <- as.data.frame(colData(airway))

  message("Airway count matrix: ", nrow(count_matrix), " genes x ",
          ncol(count_matrix), " samples")
  message("Sample conditions: ",
          paste(col_data$dex, collapse = ", "))

  write.csv(count_matrix,
            file = file.path(counts_dir, "raw_counts.csv"))

  write.csv(col_data,
            file = file.path(counts_dir, "sample_metadata.csv"))
}

message("\nStep 03 complete. Count matrix saved to results/counts/")
