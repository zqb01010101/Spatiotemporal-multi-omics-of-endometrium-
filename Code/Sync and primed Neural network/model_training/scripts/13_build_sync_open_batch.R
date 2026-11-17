#!/usr/bin/env Rscript
# =============================================================================
# 13_build_sync_open_batch.R — Add near-threshold synchronous opening training samples
#
# Select the 100 lowest scores among candidates with score_sync_open > 7.0.
# Mark queued candidates for review, enqueue indexed candidates, and skip unindexed candidates.
# Output: model_training/data/sync_open_batch_v7.csv
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

for (f in c("load_inputs.R", "database.R")) source(file.path(app_dir, "R", f))
cfg <- read_config(app_dir)

SCORE_COL  <- "score_sync_open"
THRESHOLD  <- 7.0
N_TAKE     <- 100L
BATCH_ID   <- "sync_open_v7"
ANNOTATOR  <- "sync_open_v7"

con <- init_db(cfg)
scores <- data.table::fread(file.path(model_dir, "outputs",
                                      "all_candidate_scores_with_context.csv"))
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL")$candidate_id
cand_tab <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))

sub <- scores[get(SCORE_COL) > THRESHOLD]
data.table::setorderv(sub, SCORE_COL)
batch <- head(sub, N_TAKE)
batch[, batch_id := BATCH_ID]
batch[, batch_score := as.numeric(get(SCORE_COL))]
batch[, source := paste0(SCORE_COL, "_gt_", THRESHOLD, "_lowest")]
batch <- batch[, .(candidate_id, context, peak, gene, batch_id, batch_score, source)]

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

data.table::fwrite(batch, file.path(model_dir, "data", "sync_open_batch_v7.csv"))
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
cat("score_sync_open 范围：", min(batch$batch_score), "-", max(batch$batch_score), "\n")
cat("\n已写出：", file.path(model_dir, "data", "sync_open_batch_v7.csv"), "\n")
