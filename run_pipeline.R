############################################################
# run_pipeline.R — Master Runner
# Author: Shafiya Sakina
#
# Usage:
#   Edit config.R, then run:
#     source("run_pipeline.R")
#
#   Or from terminal:
#     Rscript run_pipeline.R
#
#   Or run individual steps:
#     source("config.R")
#     source("scripts/04_normalization.R")
############################################################

start_time <- Sys.time()

cat("\n")
cat("╔══════════════════════════════════════════════╗\n")
cat("║     Bulk RNA-seq Pipeline  |  Shafiya Sakina ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

# Load config first
source("config.R")

# Define steps
steps <- list(
  list(id = "01", name = "Quality Control",          file = "scripts/01_qc.R"),
  list(id = "02", name = "Alignment",                file = "scripts/02_alignment.R"),
  list(id = "03", name = "Quantification",           file = "scripts/03_quantification.R"),
  list(id = "04", name = "Normalization",            file = "scripts/04_normalization.R"),
  list(id = "05", name = "Batch Correction",         file = "scripts/05_batch_correction.R"),
  list(id = "06", name = "Differential Expression",  file = "scripts/06_deseq2_analysis.R"),
  list(id = "07", name = "Visualization",            file = "scripts/07_visualization.R")
)

# Track results
step_log <- data.frame(
  step = character(), name = character(),
  status = character(), time_s = numeric(),
  stringsAsFactors = FALSE
)

# Run each step
for (s in steps) {
  cat(sprintf("\n[%s/07] %s\n", s$id, s$name))
  cat(strrep("-", 46), "\n")

  t0 <- Sys.time()

  tryCatch({
    source(s$file)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("Done in %.1fs\n", elapsed))
    step_log <- rbind(step_log, data.frame(
      step = s$id, name = s$name,
      status = "SUCCESS", time_s = elapsed
    ))
  }, error = function(e) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    step_log <<- rbind(step_log, data.frame(
      step = s$id, name = s$name,
      status = paste0("FAILED: ", conditionMessage(e)),
      time_s = elapsed
    ))
  })
}

# Summary
total_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)

cat("\n")
cat("╔══════════════════════════════════════════════╗\n")
cat("║               Pipeline Summary              ║\n")
cat("╚══════════════════════════════════════════════╝\n")

for (i in seq_len(nrow(step_log))) {
  icon <- if (grepl("SUCCESS", step_log$status[i])) "OK" else "FAIL"
  cat(sprintf("  [%s] [%s] %-28s  %5.1fs\n",
              icon, step_log$step[i],
              step_log$name[i], step_log$time_s[i]))
}

cat(sprintf("\n  Total time: %.1f min\n", total_time))
cat(sprintf("  Outputs:    %s\n\n", RESULTS_DIR))

write.csv(step_log,
          file.path(RESULTS_DIR, "pipeline_run_log.csv"),
          row.names = FALSE)
