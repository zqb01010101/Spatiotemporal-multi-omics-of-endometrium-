#!/usr/bin/env Rscript
# =============================================================================
# 18_export_organoid_training.R — Export stromal training data with organoid avg
#
# Read Y_sig and the new organoid avg. Keep labeled stromal candidates
# with all four scores for organoid training.
# Output: model_training/organoid_run/labeled_training_data_organoid.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "",
                                        commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
model_dir  <- normalizePath(file.path(script_dir, ".."))
app_dir    <- normalizePath(file.path(model_dir, "..", "peak_gene_label_app"))
root_dir   <- dirname(app_dir)
out_dir    <- file.path(model_dir, "organoid_run")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in c("load_inputs.R", "extract_pair_values.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}
cfg <- read_config(app_dir)

Y_sig <- readRDS(file.path(root_dir, "26_08_06_Y_sig.rds"))
avg   <- readRDS(file.path(root_dir,
  "26_09_03_avg_merged_list_change_PO9_11_13_to_organoid.rds"))
inputs <- list(Y_sig = Y_sig, avg = avg)
rm(Y_sig, avg); gc()

con <- init_db(cfg)
ann <- as.data.table(DBI::dbGetQuery(con, "
  SELECT c.candidate_id, c.context, c.peak, c.gene,
         a.score_pre_open, a.score_sync_open, a.score_sync_repress,
         a.score_unrelated, a.overall_confidence, c.batch_id
  FROM candidates c
  JOIN annotations a ON a.candidate_id = c.candidate_id
  WHERE c.context = 'stroma'
    AND a.score_pre_open IS NOT NULL
    AND a.score_sync_open IS NOT NULL
    AND a.score_sync_repress IS NOT NULL
    AND a.score_unrelated IS NOT NULL"))
DBI::dbDisconnect(con)
cat("stroma 已标注样本：", nrow(ann), "\n")

parts <- list()
CHUNK <- 100000L
for (k in seq_len(ceiling(nrow(ann) / CHUNK))) {
  idx <- seq((k - 1L) * CHUNK + 1L, min(k * CHUNK, nrow(ann)))
  vals <- extract_bulk_values(inputs, cfg, "stroma", ann$peak[idx], ann$gene[idx])
  vals[, key_ := NULL]
  parts[[k]] <- cbind(ann[idx, .(candidate_id, context, peak, gene,
                                  score_pre_open, score_sync_open,
                                  score_sync_repress, score_unrelated,
                                  overall_confidence, batch_id)], vals)
}
out <- data.table::rbindlist(parts)
data.table::fwrite(out, file.path(out_dir, "labeled_training_data_organoid.csv"))
cat("已写出：", file.path(out_dir, "labeled_training_data_organoid.csv"),
    "（", nrow(out), " 行）\n")
