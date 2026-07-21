############################################################
# Bulk RNA-seq Pipeline | Step 06: Differential Expression
# Author: Shafiya Sakina
# Dataset: GSE96870 / airway (dexamethasone treated vs untreated)
# Description: Run DESeq2 Wald test, apply apeglm LFC shrinkage,
#              and export significant DEG tables.
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
})

# ----- Paths -----
project_dir <- getwd()
counts_dir  <- file.path(project_dir, "results", "counts")
deg_dir     <- file.path(project_dir, "results", "DEGs")

dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Parameters -----
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0      # |log2FoldChange| threshold

# ----- Load dds -----
dds_path <- file.path(counts_dir, "dds_corrected.rds")
if (!file.exists(dds_path)) {
  dds_path <- file.path(counts_dir, "dds_raw.rds")
  message("Corrected dds not found, using raw dds.")
}
if (!file.exists(dds_path)) stop("Run 04_normalization.R first.")

dds <- readRDS(dds_path)

# ----- Run DESeq2 -----
message("Running DESeq2 (Wald test)...")
dds <- DESeq(dds)

resultsNames(dds)

# ----- Extract results: treated vs untreated -----
res_raw <- results(
  dds,
  contrast  = c("dex", "trt", "untrt"),
  alpha     = PADJ_CUTOFF
)

# ----- LFC shrinkage (apeglm — reduces noise for low-count genes) -----
# Use the coefficient name from resultsNames(dds)
coef_name <- grep("dex_trt_vs_untrt", resultsNames(dds), value = TRUE)

res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  type = "apeglm"
)

message("DESeq2 summary (shrunken LFC):")
summary(res_shrunk, alpha = PADJ_CUTOFF)

# ----- Export results -----
# All genes
res_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj)

write.csv(res_df,
          file = file.path(deg_dir, "all_genes_DESeq2.csv"),
          row.names = FALSE)

# Significant DEGs only
sig_degs <- res_df %>%
  filter(!is.na(padj),
         padj < PADJ_CUTOFF,
         abs(log2FoldChange) >= LFC_CUTOFF)

write.csv(sig_degs,
          file = file.path(deg_dir, "significant_DEGs.csv"),
          row.names = FALSE)

message("\nResults summary:")
message("  Total genes tested:     ", nrow(res_df))
message("  Significant DEGs:       ", nrow(sig_degs))
message("  Upregulated (trt):      ",
        sum(sig_degs$log2FoldChange > 0))
message("  Downregulated (trt):    ",
        sum(sig_degs$log2FoldChange < 0))

# Save dds with DESeq2 results for visualization step
saveRDS(dds,          file.path(counts_dir, "dds_deseq2.rds"))
saveRDS(res_shrunk,   file.path(deg_dir,    "res_shrunk.rds"))

message("\nStep 06 complete. DEG tables saved to results/DEGs/")
