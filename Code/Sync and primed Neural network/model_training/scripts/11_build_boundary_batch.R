#!/usr/bin/env Rscript
# =============================================================================
# 11_build_boundary_batch.R — Add training labels near decision thresholds
#
# Selection rules:
#   1) Select the 100 lowest scores among candidates with score_sync_open > 7.5
#   2) Select the 50 lowest scores among candidates with score_pre_open > 7.5
# Database rules:
#   - Mark queued candidates for review (review_status='review').
#   - Add indexed candidates outside the queue as a new batch (queue_priority=0).
#   - Log and skip candidates absent from the candidates table.
#
# Output: model_training/data/boundary_batch_v6.csv
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

SPECS <- list(
  list(col = "score_sync_open", n = 100L, batch = "boundary_v6_sync_open"),
  list(col = "score_pre_open",  n = 50L,  batch = "boundary_v6_pre_open")
)
THRESHOLD <- 7.5

con <- init_db(cfg)
scores <- data.table::fread(file.path(model_dir, "outputs",
                                      "all_candidate_scores_with_context.csv"))
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL")$candidate_id
cand_tab <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))

selected <- list()
seen <- character()
for (sp in SPECS) {
  sub <- scores[get(sp$col) > THRESHOLD]
  sub <- sub[!candidate_id %in% seen]
  data.table::setorderv(sub, sp$col)
  take <- head(sub, sp$n)
  take[, batch_id := sp$batch]
  take[, batch_score := as.numeric(get(sp$col))]
  take[, source := paste0(sp$col, "_gt_", THRESHOLD, "_lowest")]
  take <- take[, .(candidate_id, context, peak, gene, batch_id, batch_score, source)]
  selected[[sp$col]] <- take
  seen <- c(seen, take$candidate_id)
}

batch <- data.table::rbindlist(selected, use.names = TRUE, fill = TRUE)
batch <- unique(batch, by = "candidate_id")
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
    mark_review(con, batch_valid$candidate_id[i], annotator = "boundary_v6")
  }
}

data.table::fwrite(batch, file.path(model_dir, "data", "boundary_batch_v6.csv"))
check <- DBI::dbGetQuery(con, "
  SELECT batch_id,
         COUNT(*) AS n_total,
         SUM(CASE WHEN a.review_status = 'review' THEN 1 ELSE 0 END) AS n_review
  FROM candidates c
  LEFT JOIN annotations a USING (candidate_id)
  WHERE c.batch_id IN ('boundary_v6_sync_open', 'boundary_v6_pre_open')
  GROUP BY batch_id")
print(check)
DBI::dbDisconnect(con)

cat("\n选择统计：\n")
print(batch[, .N, by = .(source, status)])
cat("\n已写出：", file.path(model_dir, "data", "boundary_batch_v6.csv"), "\n")
