############################################################
# Bulk RNA-seq Pipeline | Step 01: Quality Control
# Author: Shafiya Sakina
# Dataset: GSE96870 (Human airway smooth muscle cells)
# Description: Run FastQC on raw FASTQ files and generate
#              a MultiQC summary report.
############################################################

suppressPackageStartupMessages({
  library(fastqcr)
})

# ----- Paths -----
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data")
qc_dir      <- file.path(project_dir, "results", "fastqc")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)

# ----- Download test dataset (airway — Bioconductor) -----
# For a quick demo without downloading raw FASTQs,
# we use the pre-built airway SummarizedExperiment.
# Skip to 02_alignment.R if you have your own FASTQs.

if (!requireNamespace("airway", quietly = TRUE))
  BiocManager::install("airway")

library(airway)
data(airway)
message("Airway dataset loaded: ", nrow(airway), " genes x ",
        ncol(airway), " samples")
message("Conditions: ", paste(unique(airway$dex), collapse = ", "))

# ----- FastQC on raw FASTQs (if available) -----
fastq_files <- list.files(data_dir,
                           pattern = "\\.fastq\\.gz$",
                           full.names = TRUE)

if (length(fastq_files) > 0) {
  message("Running FastQC on ", length(fastq_files), " files...")
  fastqc_install()
  fastqc(fq.dir = data_dir, outdir = qc_dir, threads = 4)

  # Aggregate QC results
  qc_report <- qc_aggregate(qc_dir)
  print(qc_report)

  # Flag any samples that fail basic QC
  failed <- qc_report[qc_report$status == "FAIL", ]
  if (nrow(failed) > 0) {
    warning("The following samples failed QC:\n",
            paste(failed$sample, collapse = "\n"))
  } else {
    message("All samples passed QC.")
  }
} else {
  message("No FASTQ files found in data/. ",
          "Using airway Bioconductor package for downstream steps.")
}

message("\nStep 01 complete. Check results/fastqc/ for QC reports.")
