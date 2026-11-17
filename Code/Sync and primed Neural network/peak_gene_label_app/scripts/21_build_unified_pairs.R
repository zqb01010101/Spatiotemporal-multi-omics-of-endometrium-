#!/usr/bin/env Rscript
# =============================================================================
# 21_build_unified_pairs.R — Build the unified pair table
#
# Each row describes one peak–gene pair to score:
#   Metadata: candidate_id, context, peak, gene, is_organoid, batch_id,
#           batch_score, review_status, four scores, and overall_confidence
#   Features: peak_PO3..peak_PO13 and gene_PO3..gene_PO13 (12 values)
#
# Set is_organoid by these rules:
#   - A batch_id starting with organoid_ indicates an organoid candidate.
#   - A pair present only in the organoid features also indicates an organoid candidate.
#   - Otherwise classify the pair as nonorganoid.
#
# Output: data_cache/pair_features_all.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "",
                                        commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
app_dir   <- normalizePath(dirname(script_dir))

for (f in c("load_inputs.R", "database.R")) source(file.path(app_dir, "R", f))
cfg <- read_config(app_dir)

con <- init_db(cfg)
queue <- as.data.table(DBI::dbGetQuery(con, "
  SELECT candidate_id, context, peak, gene, batch_id, batch_score, queue_priority
  FROM candidates WHERE queue_priority IS NOT NULL"))
ann <- as.data.table(DBI::dbGetQuery(con, "
  SELECT candidate_id, score_pre_open, score_sync_open, score_sync_repress,
         score_unrelated, overall_confidence, review_status
  FROM annotations"))
DBI::dbDisconnect(con)

nonorg <- data.table::fread(file.path(app_dir, "data_cache",
                                      "queue_features_nonorganoid.csv"))
org    <- data.table::fread(file.path(app_dir, "data_cache",
                                      "queue_features_organoid.csv"))

organoid_batches <- unique(queue$batch_id[grepl("^organoid_", queue$batch_id)])
queue[, is_organoid := data.table::fifelse(
  batch_id %in% organoid_batches |
    (candidate_id %in% org$candidate_id & !candidate_id %in% nonorg$candidate_id),
  1L, 0L)]

meta <- queue[, .(candidate_id, context, peak, gene, is_organoid,
                  batch_id, batch_score)]
meta[ann, on = "candidate_id", `:=`(
  score_pre_open = i.score_pre_open,
  score_sync_open = i.score_sync_open,
  score_sync_repress = i.score_sync_repress,
  score_unrelated = i.score_unrelated,
  overall_confidence = i.overall_confidence,
  review_status = i.review_status)]

org_rows <- meta[is_organoid == 1L]
non_rows <- meta[is_organoid == 0L]
org_feat <- org[org_rows, on = "candidate_id", nomatch = 0]
non_feat <- nonorg[non_rows, on = "candidate_id", nomatch = 0]

cols <- c("candidate_id", "context", "peak", "gene", "is_organoid",
          "batch_id", "batch_score", "review_status",
          "score_pre_open", "score_sync_open", "score_sync_repress",
          "score_unrelated", "overall_confidence",
          paste0("peak_", cfg$stages$labels), paste0("gene_", cfg$stages$labels))
unified <- data.table::rbindlist(list(non_feat[, ..cols], org_feat[, ..cols]),
                                 use.names = TRUE, fill = TRUE)
unified <- unique(unified, by = "candidate_id")
data.table::setorder(unified, candidate_id)

outfile <- file.path(app_dir, "data_cache", "pair_features_all.csv")
data.table::fwrite(unified, outfile)
cat("统一数据文件：", outfile, "（", nrow(unified), " 行）\n")
cat("is_organoid 分布：\n")
print(unified[, .N, by = is_organoid])
