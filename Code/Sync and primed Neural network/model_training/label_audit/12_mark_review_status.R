#!/usr/bin/env Rscript
# =============================================================================
# 12_mark_review_status.R — Mark batch candidates as pending review in the app
#
# Usage: Rscript label_audit/12_mark_review_status.R [batch CSV filename]
#       Default: model_training/data/review_batch_audit_v1.csv
# Call mark_review for each candidate to set review_status='review' and retain prior scores.
#       Use the app's pending-review filter to find candidates with that status.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "",
                                        commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
model_dir <- normalizePath(file.path(script_dir, ".."))
app_dir   <- normalizePath(file.path(model_dir, "..", "peak_gene_label_app"))

for (f in c("load_inputs.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}
cfg <- read_config(app_dir)

args <- commandArgs(trailingOnly = TRUE)
batch_name <- if (length(args)) args[1] else "review_batch_audit_v1.csv"
batch <- data.table::fread(file.path(model_dir, "data", batch_name))
cat("批次行数：", nrow(batch), "\n")

con <- init_db(cfg)
for (i in seq_len(nrow(batch))) {
  mark_review(con, batch$candidate_id[i], annotator = "audit_review_v1")
}
check <- DBI::dbGetQuery(con, "
  SELECT batch_id, COUNT(*) AS n_review_marked
  FROM candidates
  JOIN annotations USING (candidate_id)
  WHERE batch_id IS NOT NULL AND review_status = 'review'
  GROUP BY batch_id")
print(check[grepl("^audit_v1", check$batch_id), ])
DBI::dbDisconnect(con)
cat("全部", nrow(batch), "条候选已标记为待复核。\n")
