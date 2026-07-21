############################################################
# Bulk RNA-seq Pipeline | Step 07: Visualization
# Author: Shafiya Sakina
# Dataset: GSE96870 / airway
# Description: Generate publication-ready plots:
#   - Volcano plot (EnhancedVolcano)
#   - MA plot (apeglm-shrunken)
#   - Heatmap of top 50 DEGs (pheatmap)
#   - Dispersion estimates
#   - GO enrichment dot plot (clusterProfiler)
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

# ----- Paths -----
project_dir <- getwd()
counts_dir  <- file.path(project_dir, "results", "counts")
deg_dir     <- file.path(project_dir, "results", "DEGs")
plots_dir   <- file.path(project_dir, "results", "plots")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Load objects -----
dds       <- readRDS(file.path(counts_dir, "dds_deseq2.rds"))
res       <- readRDS(file.path(deg_dir,    "res_shrunk.rds"))
sig_degs  <- read.csv(file.path(deg_dir,   "significant_DEGs.csv"))

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0

# =========================================================
# 1. Volcano Plot
# =========================================================
message("Generating volcano plot...")

EnhancedVolcano(
  res,
  lab            = rownames(res),
  x              = "log2FoldChange",
  y              = "padj",
  pCutoff        = PADJ_CUTOFF,
  FCcutoff       = LFC_CUTOFF,
  title          = "Dexamethasone treated vs untreated",
  subtitle       = paste0("DESeq2 | apeglm LFC shrinkage | ",
                           "padj < ", PADJ_CUTOFF,
                           ", |LFC| ≥ ", LFC_CUTOFF),
  caption        = paste0("Total DEGs: ", nrow(sig_degs)),
  legendPosition = "bottom",
  col            = c("grey70", "grey70", "#2166AC", "#D73027"),
  colAlpha       = 0.7,
  pointSize      = 2.0,
  labSize        = 3.5,
  drawConnectors = TRUE,
  widthConnectors = 0.4
)

ggsave(file.path(plots_dir, "volcano_plot.png"),
       width = 9, height = 7, dpi = 150)
message("  Saved: volcano_plot.png")

# =========================================================
# 2. MA Plot
# =========================================================
message("Generating MA plot...")

png(file.path(plots_dir, "ma_plot.png"),
    width = 800, height = 600, res = 120)
plotMA(res,
       alpha = PADJ_CUTOFF,
       main  = "MA Plot — apeglm LFC shrinkage",
       ylim  = c(-5, 5),
       colSig = "#D73027",
       colNonSig = "grey70")
abline(h = c(-LFC_CUTOFF, LFC_CUTOFF), lty = 2, col = "steelblue")
dev.off()
message("  Saved: ma_plot.png")

# =========================================================
# 3. Dispersion Plot
# =========================================================
message("Generating dispersion plot...")

png(file.path(plots_dir, "dispersion_plot.png"),
    width = 700, height = 550, res = 120)
plotDispEsts(dds,
             main = "DESeq2 Dispersion Estimates")
dev.off()
message("  Saved: dispersion_plot.png")

# =========================================================
# 4. Heatmap — Top 50 DEGs
# =========================================================
message("Generating DEG heatmap...")

top50 <- sig_degs %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

# Use VST-normalized matrix
vst_mat <- assay(vst(dds, blind = FALSE))
top50   <- top50[top50 %in% rownames(vst_mat)]

if (length(top50) > 1) {
  heat_mat <- vst_mat[top50, ]

  # Z-score scale across samples
  heat_mat_scaled <- t(scale(t(heat_mat)))

  anno_col <- data.frame(
    condition = as.character(dds$dex),
    row.names = colnames(dds)
  )

  ann_colors <- list(
    condition = c(untrt = "#2166AC", trt = "#D73027")
  )

  png(file.path(plots_dir, "heatmap_top50.png"),
      width = 900, height = 1000, res = 120)
  pheatmap(
    heat_mat_scaled,
    annotation_col  = anno_col,
    annotation_colors = ann_colors,
    color           = colorRampPalette(
                        rev(brewer.pal(11, "RdBu")))(100),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    show_rownames   = TRUE,
    show_colnames   = TRUE,
    fontsize_row    = 7,
    main            = "Top 50 DEGs (VST, z-score scaled)"
  )
  dev.off()
  message("  Saved: heatmap_top50.png")
} else {
  message("  Not enough significant DEGs for heatmap.")
}

# =========================================================
# 5. GO Enrichment — Biological Process
# =========================================================
message("Running GO enrichment analysis...")

# Convert Ensembl IDs to Entrez (required by clusterProfiler)
gene_list <- sig_degs$gene_id

entrez_ids <- bitr(
  gene_list,
  fromType = "ENSEMBL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

if (nrow(entrez_ids) > 5) {
  go_results <- enrichGO(
    gene          = entrez_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",              # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  if (!is.null(go_results) && nrow(go_results) > 0) {
    # Dot plot
    go_dot <- dotplot(go_results,
                       showCategory = 15,
                       title = "GO Biological Process Enrichment",
                       font.size = 10)

    ggsave(file.path(plots_dir, "go_dotplot.png"),
           go_dot, width = 10, height = 7, dpi = 150)
    message("  Saved: go_dotplot.png")

    # Save GO table
    write.csv(as.data.frame(go_results),
              file.path(deg_dir, "GO_enrichment_BP.csv"),
              row.names = FALSE)
  } else {
    message("  No significant GO terms found.")
  }
} else {
  message("  Too few genes for GO enrichment.")
}

message("\n=================================================")
message("Step 07 complete. All plots saved to results/plots/")
message("  - pca_plot.png")
message("  - sample_distance_heatmap.png")
message("  - volcano_plot.png")
message("  - ma_plot.png")
message("  - dispersion_plot.png")
message("  - heatmap_top50.png")
message("  - batch_correction_pca.png")
message("  - go_dotplot.png")
message("=================================================\n")
