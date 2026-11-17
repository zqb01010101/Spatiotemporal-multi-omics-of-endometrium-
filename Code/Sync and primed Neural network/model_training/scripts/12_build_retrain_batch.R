#!/usr/bin/env Rscript
# =============================================================================
# 12_build_retrain_batch.R — Import retraining samples from 26_08_18_re_train_sample.rds
#
# Input columns include candidate_id, context, peak, gene, and score_pre_open.
# Database rules:
#   - Mark queued candidates for review (review_status='review').
#   - Add indexed candidates outside the queue as a new batch (queue_priority=0).
#   - Log and skip candidates absent from the candidates table.
#
# Output: model_training/data/retrain_batch_v1.csv
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

for (f in c("load_inputs.R", "database.R")) source(file.path(app_dir, "R", f))
cfg <- read_config(app_dir)

dat <- as.data.table(readRDS(file.path(root_dir, "26_08_18_re_train_sample.rds")))
required <- c("candidate_id", "context", "peak", "gene", "score_pre_open")
stopifnot(all(required %in% names(dat)))
cat("输入行数：", nrow(dat), "\n")
cat("context 分布：\n")
print(dat[, .N, by = context])
cat("score_pre_open 范围：", min(dat$score_pre_open), "-", max(dat$score_pre_open), "\n")

con <- init_db(cfg)
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL")$candidate_id
cand_tab <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))

batch <- dat[, .(candidate_id, context, peak, gene,
                 batch_score = as.numeric(score_pre_open))]
batch[, batch_id := "retrain_v1_pre_open"]
batch[, source := "26_08_18_re_train_sample"]
batch[, in_candidates := candidate_id %in% cand_tab$candidate_id]
batch[, in_queue := candidate_id %in% queue_ids]
batch[, status := data.table::fifelse(
  !in_candidates, "excluded",
  data.table::fifelse(in_queue, "review", "new"))]

batch_valid <- batch[status != "excluded"]
DBI::dbBegin(con)
for (i in seq_len(nrow(batch_valid))) {
  r <- batch_valid[i]
  DBI::dbExecute(con, "
    UPDATE candidates
    SET batch_id = ?, batch_score = ?, queue_priority = 0
    WHERE candidate_id = ?",
    params = list(r$batch_id, r$batch_score, r$candidate_id))
}
DBI::dbCommit(con)
for (i in seq_len(nrow(batch_valid))) {
  if (batch_valid$status[i] == "review") {
    mark_review(con, batch_valid$candidate_id[i], annotator = "retrain_v1")
  }
}

data.table::fwrite(batch, file.path(model_dir, "data", "retrain_batch_v1.csv"))
check <- DBI::dbGetQuery(con, "
  SELECT batch_id, COUNT(*) AS n_total,
         SUM(CASE WHEN a.review_status = 'review' THEN 1 ELSE 0 END) AS n_review
  FROM candidates c
  LEFT JOIN annotations a USING (candidate_id)
  WHERE c.batch_id = 'retrain_v1_pre_open'
  GROUP BY batch_id")
print(check)
DBI::dbDisconnect(con)

cat("\n状态统计：\n")
print(batch[, .N, by = status])
cat("\n已写出：", file.path(model_dir, "data", "retrain_batch_v1.csv"), "\n")
