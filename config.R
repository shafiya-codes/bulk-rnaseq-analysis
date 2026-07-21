############################################################
# config.R — Master Configuration
# Author: Shafiya Sakina
#
# ✅ THIS IS THE ONLY FILE YOU NEED TO EDIT TO REUSE
#    THE PIPELINE ON A NEW DATASET.
#
# After editing, run:
#   source("run_pipeline.R")
############################################################


# ============================================================
# 1. PROJECT
# ============================================================
PROJECT_NAME  <- "airway_dex_study"   # Used in plot titles & output folder names
AUTHOR        <- "Shafiya Sakina"


# ============================================================
# 2. DATASET
# ============================================================

# --- Option A: Bioconductor package (quickest, no download) ---
# Supported: "airway" (human), "pasilla" (Drosophila), "fission" (yeast)
USE_BIOC_DATASET <- TRUE
BIOC_DATASET     <- "airway"          # change to "pasilla" or "fission"

# --- Option B: GEO accession (downloads count matrix via GEOquery) ---
# Set USE_BIOC_DATASET <- FALSE and fill these in:
GEO_ACCESSION    <- "GSE96870"
GEO_COUNT_FILE   <- ""                # filename inside GEO archive, or "" to auto-detect

# --- Option C: Your own local files ---
# Set USE_BIOC_DATASET <- FALSE and GEO_ACCESSION <- ""
LOCAL_COUNT_FILE    <- "data/counts.csv"
LOCAL_METADATA_FILE <- "data/metadata.csv"


# ============================================================
# 3. ORGANISM
# ============================================================
# Human  → "Hs"  | Mouse → "Mm" | Yeast → "Sc" | Fly → "Dm"

ORGANISM <- "Hs"

ORGANISM_CONFIG <- list(
  Hs = list(orgdb = "org.Hs.eg.db", from_type = "ENSEMBL",
            kegg_org = "hsa", latin = "Homo sapiens"),
  Mm = list(orgdb = "org.Mm.eg.db", from_type = "ENSEMBL",
            kegg_org = "mmu", latin = "Mus musculus"),
  Sc = list(orgdb = "org.Sc.sgd.db", from_type = "ORF",
            kegg_org = "sce", latin = "Saccharomyces cerevisiae"),
  Dm = list(orgdb = "org.Dm.eg.db", from_type = "ENSEMBL",
            kegg_org = "dme", latin = "Drosophila melanogaster")
)


# ============================================================
# 4. EXPERIMENTAL DESIGN
# ============================================================
CONDITION_COL    <- "dex"
CONDITION_LEVELS <- c("untrt", "trt")   # control FIRST
BATCH_COL        <- "cell"              # set "" to skip batch correction
REFERENCE_LEVEL  <- "untrt"


# ============================================================
# 5. DEG THRESHOLDS
# ============================================================
PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 1.0
LFC_SHRINKAGE <- "apeglm"   # "apeglm", "ashr", or "normal"
MIN_COUNT     <- 10
MIN_SAMPLES   <- 2


# ============================================================
# 6. VISUALIZATION
# ============================================================
TOP_N_HEATMAP    <- 50
TOP_N_VOLCANO    <- 20
PLOT_FORMAT      <- "png"
PLOT_DPI         <- 150
PLOT_WIDTH       <- 9
PLOT_HEIGHT      <- 7
CONDITION_COLORS <- c("#2166AC", "#D73027")   # control → treatment


# ============================================================
# 7. REFERENCE GENOME (only if aligning raw FASTQs)
# ============================================================
GENOME_FA  <- "reference/genome.fa"
GTF_FILE   <- "reference/annotation.gtf"
INDEX_BASE <- "reference/genome_index"
N_THREADS  <- 8


# ============================================================
# 8. OUTPUT PATHS — auto-generated, no need to edit
# ============================================================
RESULTS_DIR <- file.path("results", PROJECT_NAME)
COUNTS_DIR  <- file.path(RESULTS_DIR, "counts")
DEG_DIR     <- file.path(RESULTS_DIR, "DEGs")
PLOTS_DIR   <- file.path(RESULTS_DIR, "plots")
QC_DIR      <- file.path(RESULTS_DIR, "fastqc")

invisible(lapply(
  c(RESULTS_DIR, COUNTS_DIR, DEG_DIR, PLOTS_DIR, QC_DIR),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

ORG <- ORGANISM_CONFIG[[ORGANISM]]

message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
message("  Config loaded: ", PROJECT_NAME)
message("  Organism:      ", ORG$latin)
message("  Dataset:       ",
        if (USE_BIOC_DATASET) BIOC_DATASET
        else if (nchar(GEO_ACCESSION) > 0) GEO_ACCESSION
        else "local files")
message("  Contrast:      ", CONDITION_LEVELS[2], " vs ", CONDITION_LEVELS[1])
message("  padj < ", PADJ_CUTOFF, "  |LFC| >= ", LFC_CUTOFF)
message("  Output:        ", RESULTS_DIR)
message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
