############################################################
# Bulk RNA-seq Pipeline | Step 05: Batch Correction
# Author: Shafiya Sakina
# Dataset: GSE96870 / airway
# Description: Detect and correct for batch effects using
#              ComBat-seq (count-aware). Updates dds design
#              to include batch as a covariate.
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(ggplot2)
  library(limma)
})

# ----- Paths -----
project_dir <- getwd()
counts_dir  <- file.path(project_dir, "results", "counts")
plots_dir   <- file.path(project_dir, "results", "plots")

# ----- Load dds from previous step -----
dds_path <- file.path(counts_dir, "dds_raw.rds")
if (!file.exists(dds_path)) stop("Run 04_normalization.R first.")
dds <- readRDS(dds_path)

count_matrix <- counts(dds)
col_data     <- as.data.frame(colData(dds))

# ----- Check for batch variable -----
# The airway dataset has a 'cell' column (4 cell lines) that acts as batch.
# In your own data, replace 'cell' with your actual batch variable.

if ("cell" %in% colnames(col_data)) {
  batch     <- col_data$cell
  condition <- col_data$dex

  message("Batch variable detected: cell line (", 
          paste(unique(batch), collapse = ", "), ")")

  # ----- ComBat-seq: count-aware batch correction -----
  message("Running ComBat-seq...")
  corrected_counts <- ComBat_seq(
    counts    = count_matrix,
    batch     = batch,
    group     = condition
  )

  # Rebuild dds with corrected counts and batch in design
  dds_corrected <- DESeqDataSetFromMatrix(
    countData = corrected_counts,
    colData   = col_data,
    design    = ~ cell + dex     # cell = batch covariate
  )

  message("Batch correction applied. New design: ~ cell + dex")

  # ----- Visualize batch effect before/after -----
  # Before
  vst_before <- vst(dds, blind = TRUE)
  pca_before <- plotPCA(vst_before, intgroup = c("dex", "cell"),
                         returnData = TRUE)
  pct_before <- round(100 * attr(pca_before, "percentVar"))

  p_before <- ggplot(pca_before,
                      aes(PC1, PC2, color = dex, shape = cell)) +
    geom_point(size = 4) +
    labs(title = "PCA before batch correction",
         x = paste0("PC1: ", pct_before[1], "%"),
         y = paste0("PC2: ", pct_before[2], "%")) +
    theme_bw(base_size = 12)

  # After
  vst_after  <- vst(dds_corrected, blind = TRUE)
  pca_after  <- plotPCA(vst_after, intgroup = c("dex", "cell"),
                         returnData = TRUE)
  pct_after  <- round(100 * attr(pca_after, "percentVar"))

  p_after <- ggplot(pca_after,
                     aes(PC1, PC2, color = dex, shape = cell)) +
    geom_point(size = 4) +
    labs(title = "PCA after ComBat-seq batch correction",
         x = paste0("PC1: ", pct_after[1], "%"),
         y = paste0("PC2: ", pct_after[2], "%")) +
    theme_bw(base_size = 12)

  library(patchwork)
  combined <- p_before + p_after
  ggsave(file.path(plots_dir, "batch_correction_pca.png"),
         combined, width = 12, height = 5, dpi = 150)
  message("Batch correction PCA comparison saved.")

  saveRDS(dds_corrected,
          file = file.path(counts_dir, "dds_corrected.rds"))
  message("Corrected dds saved to results/counts/dds_corrected.rds")

} else {
  message("No batch variable found in metadata. ",
          "Skipping ComBat-seq — saving original dds as corrected.")
  saveRDS(dds,
          file = file.path(counts_dir, "dds_corrected.rds"))
}

message("\nStep 05 complete.")
