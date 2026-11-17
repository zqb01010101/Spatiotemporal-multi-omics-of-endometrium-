#!/usr/bin/env Rscript
# =============================================================================
# 17_build_organoid_retrain_batch.R — Select training samples from 26_09_03_re_train_sample.rds
#
# Take the top 100 peak–gene pairs by S_pre_shape. Set the context to
#       stroma for organoid training.
# Mark queued candidates for review, enqueue indexed candidates, and skip unindexed candidates.
# Output: model_training/data/organoid_retrain_v1.csv
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

N_TAKE   <- 100L
BATCH_ID <- "organoid_retrain_v1"
ANNOTATOR <- "organoid_retrain_v1"

dat <- as.data.table(readRDS(file.path(root_dir, "26_09_03_re_train_sample.rds")))
stopifnot(all(c("peak", "gene", "S_pre_shape") %in% names(dat)))
cat("输入行数：", nrow(dat), "\n")

dat <- dat[order(-S_pre_shape)][1:N_TAKE]
dat[, candidate_id := paste("stroma", peak, gene, sep = "|")]
dat[, batch_score := as.numeric(S_pre_shape)]
batch <- dat[, .(candidate_id, context = "stroma", peak, gene, batch_score)]
batch[, batch_id := BATCH_ID]
batch[, source := "26_09_03_re_train_sample_S_pre_shape_top100"]

con <- init_db(cfg)
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL")$candidate_id
cand_tab <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))

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
    mark_review(con, batch_valid$candidate_id[i], annotator = ANNOTATOR)
  }
}

data.table::fwrite(batch, file.path(model_dir, "data", "organoid_retrain_v1.csv"))
check <- DBI::dbGetQuery(con, "
  SELECT batch_id, COUNT(*) AS n_total,
         SUM(CASE WHEN a.review_status = 'review' THEN 1 ELSE 0 END) AS n_review
  FROM candidates c
  LEFT JOIN annotations a USING (candidate_id)
  WHERE c.batch_id = ?", params = list(BATCH_ID))
print(check)
DBI::dbDisconnect(con)

cat("\n状态统计：\n")
print(batch[, .N, by = status])
cat("S_pre_shape 范围：", min(batch$batch_score), "-", max(batch$batch_score), "\n")
cat("\n已写出：", file.path(model_dir, "data", "organoid_retrain_v1.csv"), "\n")
