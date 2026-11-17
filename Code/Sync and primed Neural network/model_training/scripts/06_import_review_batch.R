#!/usr/bin/env Rscript
# =============================================================================
# 06_import_review_batch.R — Import a review batch into the labeling database
#
# For each candidate in data/review_batch_*.csv:
#   * Store the batch ID and selection score in candidates.batch_id and batch_score.
#   * Set queue_priority = 0 to move the candidate to the front of the queue.
#
# Usage: Rscript scripts/06_import_review_batch.R [batch CSV filename]
#       Default: data/review_batch_v1.csv. Example override: review_batch_v2.csv.
#
# 01_build_candidate_index.R rebuilds the candidates table.
# Rerun this script afterward to restore the batch metadata.
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
batch_name <- if (length(args)) args[1] else "review_batch_v1.csv"
batch <- data.table::fread(file.path(model_dir, "data", batch_name))
cat("复核批次行数：", nrow(batch), "\n")
stopifnot(!anyDuplicated(batch$candidate_id))

con <- init_db(cfg)
DBI::dbBegin(con)
for (i in seq_len(nrow(batch))) {
  DBI::dbExecute(
    con,
    "UPDATE candidates
     SET batch_id = ?, batch_score = ?, queue_priority = 0
     WHERE candidate_id = ?",
    params = list(batch$batch_id[i], as.numeric(batch$batch_score[i]),
                  batch$candidate_id[i])
  )
}
DBI::dbCommit(con)

check <- DBI::dbGetQuery(con, "
  SELECT batch_id, COUNT(*) AS n,
         SUM(CASE WHEN queue_priority = 0 THEN 1 ELSE 0 END) AS n_priority0
  FROM candidates
  WHERE batch_id IS NOT NULL
  GROUP BY batch_id")
print(check)
missing <- batch$candidate_id[!batch$candidate_id %in%
  DBI::dbGetQuery(con, "SELECT candidate_id FROM candidates WHERE batch_id IS NOT NULL")$candidate_id]
if (length(missing)) {
  cat("未匹配到候选：", length(missing), "\n")
  print(head(missing, 10))
} else {
  cat("全部", nrow(batch), "条候选已成功标记为复核批次。\n")
}
DBI::dbDisconnect(con)
