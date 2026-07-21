############################################################
# Bulk RNA-seq Pipeline | Step 02: Alignment
# Author: Shafiya Sakina
# Dataset: GSE96870 (Human airway smooth muscle cells,
#          hg38 / GRCh38)
# Description: Build Rsubread genome index and align
#              paired-end FASTQ files to produce BAM files.
#              Skip this step if using airway count matrix.
############################################################

suppressPackageStartupMessages({
  library(Rsubread)
  library(R.utils)
})

# ----- Paths -----
project_dir <- getwd()
data_dir    <- file.path(project_dir, "data")
ref_dir     <- file.path(project_dir, "reference")
bam_dir     <- file.path(project_dir, "results", "bam")

dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bam_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Reference genome (human hg38, chr1 subset for demo) -----
# For full analysis replace with the complete GRCh38 genome.
# Full genome: ftp://ftp.ensembl.org/pub/release-109/fasta/
#              homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

genome_fa <- file.path(ref_dir, "hg38.fa")
gtf_file  <- file.path(ref_dir, "hg38.gtf")

if (!file.exists(genome_fa)) {
  stop(
    "Reference genome not found at: ", genome_fa, "\n",
    "Download from Ensembl: https://ftp.ensembl.org/pub/release-109/fasta/",
    "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz\n",
    "Then gunzip and place in reference/"
  )
}

if (!file.exists(gtf_file)) {
  stop(
    "GTF annotation not found at: ", gtf_file, "\n",
    "Download from Ensembl: https://ftp.ensembl.org/pub/release-109/gtf/",
    "homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz\n",
    "Then gunzip and place in reference/"
  )
}

# ----- Build genome index -----
index_base <- file.path(ref_dir, "hg38_index")

if (!file.exists(paste0(index_base, ".00.b.array"))) {
  message("Building Rsubread genome index (this takes ~10 min)...")
  buildindex(
    basename  = index_base,
    reference = genome_fa,
    memory    = 8000        # MB; reduce if RAM is limited
  )
  message("Index built: ", index_base)
} else {
  message("Index already exists, skipping build.")
}

# ----- Align paired-end reads -----
reads1 <- sort(list.files(data_dir, pattern = "_1\\.fastq\\.gz$",
                           full.names = TRUE))
reads2 <- sort(list.files(data_dir, pattern = "_2\\.fastq\\.gz$",
                           full.names = TRUE))

if (length(reads1) == 0) {
  message("No FASTQ files found in data/. ",
          "Skipping alignment — proceeding with airway count matrix.")
} else {
  stopifnot("Read1 and Read2 counts must match" =
              length(reads1) == length(reads2))

  message("Aligning ", length(reads1), " sample(s)...")

  align(
    index        = index_base,
    readfile1    = reads1,
    readfile2    = reads2,
    input_format = "gzFASTQ",
    output_format = "BAM",
    output_file  = file.path(bam_dir,
                              sub("_1\\.fastq\\.gz$", ".bam",
                                  basename(reads1))),
    nthreads     = 8,
    sortReadsByCoordinates = TRUE
  )

  # Mapping summary
  bam_files <- list.files(bam_dir, pattern = "\\.bam$",
                           full.names = TRUE)
  mapping_stats <- propmapped(bam_files)
  print(mapping_stats)

  write.csv(mapping_stats,
            file.path(project_dir, "results", "mapping_stats.csv"),
            row.names = FALSE)

  message("\nStep 02 complete. BAM files written to results/bam/")
}
