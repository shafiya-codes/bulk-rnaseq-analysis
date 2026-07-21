############################################################
# Bulk RNA-seq Pipeline | Step 04: Normalization
# Author: Shafiya Sakina
# Dataset: GSE96870 / airway (dex treated vs untreated)
# Description: Filter low-count genes, create DESeq2 object,
#              apply VST normalization, and run PCA.
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(airway)
})

# ----- Paths -----
project_dir <- getwd()
counts_dir  <- file.path(project_dir, "results", "counts")
plots_dir   <- file.path(project_dir, "results", "plots")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ----- Load count matrix + metadata -----
counts_file   <- file.path(counts_dir, "raw_counts.csv")
metadata_file <- file.path(counts_dir, "sample_metadata.csv")

if (file.exists(counts_file) && file.exists(metadata_file)) {
  count_matrix <- as.matrix(read.csv(counts_file, row.names = 1))
  col_data     <- read.csv(metadata_file, row.names = 1)
} else {
  message("Count files not found — loading airway directly...")
  data(airway)
  count_matrix <- assay(airway, "counts")
  col_data     <- as.data.frame(colData(airway))
}

# ----- Filter low-count genes -----
# Keep genes with >= 10 counts in at least 2 samples
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix <- count_matrix[keep, ]
message("Genes retained after filtering: ", nrow(count_matrix))

# ----- Build DESeq2 object -----
# Ensure 'condition' (dex) column is a factor
col_data$dex <- factor(col_data$dex,
                        levels = c("untrt", "trt"))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = col_data,
  design    = ~ dex
)

# ----- VST normalization -----
vst_data <- vst(dds, blind = TRUE)
message("VST normalization complete.")

# Save normalized matrix
vst_matrix <- assay(vst_data)
write.csv(vst_matrix,
          file = file.path(counts_dir, "vst_normalized.csv"))

# ----- PCA plot -----
pca_data <- plotPCA(vst_data, intgroup = "dex", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = dex, label = name)) +
  geom_point(size = 4, alpha = 0.85) +
  ggrepel::geom_text_repel(size = 3.2) +
  labs(
    title    = "PCA of VST-normalized RNA-seq samples",
    subtitle = "Colored by treatment condition (dex)",
    x        = paste0("PC1: ", pct_var[1], "% variance"),
    y        = paste0("PC2: ", pct_var[2], "% variance"),
    color    = "Treatment"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(plots_dir, "pca_plot.png"),
       pca_plot, width = 7, height = 5, dpi = 150)
message("PCA plot saved to results/plots/pca_plot.png")

# ----- Sample distance heatmap -----
library(pheatmap)
library(RColorBrewer)

sample_dists <- dist(t(vst_matrix))
dist_matrix  <- as.matrix(sample_dists)

anno_col <- data.frame(
  condition = col_data$dex,
  row.names = rownames(col_data)
)

png(file.path(plots_dir, "sample_distance_heatmap.png"),
    width = 700, height = 600, res = 120)
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = anno_col,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (VST)")
dev.off()
message("Sample distance heatmap saved.")

# Save dds for next step
saveRDS(dds, file = file.path(counts_dir, "dds_raw.rds"))
message("\nStep 04 complete. dds object saved to results/counts/dds_raw.rds")
