#!/usr/bin/env Rscript
# =============================================================================
# 10_build_sync_open_batch.R — Add synchronous opening training examples
#
# Source 1: top 100 peak–gene pairs by S_final in 26_06_20_Sync_peak_gene_stroma.rds
# Source 2: all_candidate_scores_with_context.csv from the current model.
#         Select 100 candidates with score_sync_open in [7, 8.5], prioritizing unlabeled candidates.
#
# Database rules:
#   - Mark labeled or queued candidates for review (review_status='review').
#   - Add indexed candidates that are outside the queue as a new batch (queue_priority=0).
#   - Log and skip candidates absent from the candidates table.
#
# Output: model_training/data/sync_open_batch_v6.csv with all selections and statuses
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

N_TOP_S     <- 100L
N_MODEL     <- 100L
MODEL_RANGE <- c(7.0, 8.5)

con <- init_db(cfg)

# ---------- Source 1: top 100 stromal pairs by S_final ----------
sync <- as.data.table(readRDS(file.path(root_dir, "26_06_20_Sync_peak_gene_stroma.rds")))
topS <- sync[order(-S_final)][1:N_TOP_S]
topS[, candidate_id := paste("stroma", peak, gene, sep = "|")]
topS[, batch_id := "sync_open_v6_topS"]
topS[, batch_score := as.numeric(S_final)]
topS[, source := "sync_S_final_top100"]
topS <- topS[, .(candidate_id, context = "stroma", peak, gene, batch_id, batch_score, source)]

# ---------- Source 2: model score_sync_open in [7, 8.5], prioritizing unlabeled candidates ----------
model_scores <- data.table::fread(file.path(model_dir, "outputs",
                                            "all_candidate_scores_with_context.csv"))
range_set <- model_scores[score_sync_open >= MODEL_RANGE[1] &
                          score_sync_open <= MODEL_RANGE[2]]
range_set <- range_set[!candidate_id %in% topS$candidate_id]

ann_ids <- DBI::dbGetQuery(con, "
  SELECT candidate_id FROM annotations
  WHERE score_pre_open IS NOT NULL OR score_sync_open IS NOT NULL
     OR score_sync_repress IS NOT NULL OR score_unrelated IS NOT NULL")$candidate_id
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL")$candidate_id

range_set[, annotated := candidate_id %in% ann_ids]
range_set[, in_queue := candidate_id %in% queue_ids]
data.table::setorder(range_set, annotated, -score_sync_open)
topM <- head(range_set, N_MODEL)
topM[, batch_id := "sync_open_v6_model"]
topM[, batch_score := as.numeric(score_sync_open)]
topM[, source := "model_sync_open_7_8.5"]
topM <- topM[, .(candidate_id, context, peak, gene, batch_id, batch_score, source)]

# ---------- Merge selections ----------
batch <- data.table::rbindlist(list(topS, topM), use.names = TRUE, fill = TRUE)
batch <- unique(batch, by = "candidate_id")

cand_tab <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))
batch[, in_candidates := candidate_id %in% cand_tab$candidate_id]
batch[, annotated := candidate_id %in% ann_ids]
batch[, in_queue := candidate_id %in% queue_ids]
batch[, status := data.table::fifelse(
  annotated | in_queue, "review",
  data.table::fifelse(in_candidates, "new", "excluded"))]

# ---------- Update the database ----------
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
    mark_review(con, batch_valid$candidate_id[i], annotator = "sync_open_v6")
  }
}

# ---------- Write and check the manifest ----------
data.table::fwrite(batch, file.path(model_dir, "data", "sync_open_batch_v6.csv"))
check <- DBI::dbGetQuery(con, "
  SELECT batch_id,
         COUNT(*) AS n_total,
         SUM(CASE WHEN a.review_status = 'review' THEN 1 ELSE 0 END) AS n_review
  FROM candidates c
  LEFT JOIN annotations a USING (candidate_id)
  WHERE c.batch_id IN ('sync_open_v6_topS', 'sync_open_v6_model')
  GROUP BY batch_id")
print(check)
DBI::dbDisconnect(con)

cat("\n选择统计：\n")
print(batch[, .N, by = .(source, status)])
cat("\n已写出：", file.path(model_dir, "data", "sync_open_batch_v6.csv"), "\n")
