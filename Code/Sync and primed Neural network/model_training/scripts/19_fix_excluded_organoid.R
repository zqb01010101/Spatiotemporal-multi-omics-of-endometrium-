#!/usr/bin/env Rscript
# =============================================================================
# 19_fix_excluded_organoid.R — Restore six pairs excluded by the earlier avg gene alignment
#
# TOMM6, KRTAP1-5, and TNFRSF6B have eligible pairs in res$mat and are present in organoid avg.
#       Their absence from the earlier avg excluded them from the original candidates table.
# Get TSS distances from dist$dist_mat and insert the pairs at the front of the queue as a new batch.
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

dat <- as.data.table(readRDS(file.path(root_dir, "26_09_03_re_train_sample.rds")))
top <- dat[order(-S_pre_shape)][1:100]
top[, sid := paste("stroma", peak, gene, sep = "|")]

con <- init_db(cfg)
have <- DBI::dbGetQuery(con, "SELECT candidate_id FROM candidates")$candidate_id
missing <- top[!sid %in% have]
cat("需补回：", nrow(missing), " 条\n")

dist_obj <- readRDS(file.path(root_dir, "26_08_06_PeakGene_TSSdist_1Mb.rds"))
dist_mat <- dist_obj$dist_mat
rm(dist_obj); gc()

now <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
rows <- lapply(seq_len(nrow(missing)), function(i) {
  p <- missing$peak[i]; g <- missing$gene[i]
  d <- if (p %in% rownames(dist_mat) && g %in% colnames(dist_mat)) {
    as.numeric(dist_mat[p, g])
  } else {
    NA_real_
  }
  data.table(
    candidate_id = paste("stroma", p, g, sep = "|"),
    context = "stroma", peak = p, gene = g, distance = d,
    mandatory = 0L, source_tri_pre = 0L, source_sync_open = 0L,
    source_sync_repress_heuristic = 0L, S_pre_shape = NA_real_,
    S_final = NA_real_, queue_priority = 0L, created_at = now,
    is_repeat = 0L, repeat_of = NA_character_,
    batch_id = "organoid_retrain_v1", batch_score = as.numeric(missing$S_pre_shape[i])
  )
})
rdt <- data.table::rbindlist(rows, use.names = TRUE)
DBI::dbBegin(con)
DBI::dbAppendTable(con, "candidates", rdt)
DBI::dbCommit(con)

check <- DBI::dbGetQuery(con, "
  SELECT batch_id, COUNT(*) AS n FROM candidates
  WHERE batch_id = 'organoid_retrain_v1' GROUP BY batch_id")
print(check)
DBI::dbDisconnect(con)
cat("已补回", nrow(rdt), "条到 candidates，并置顶为新批次。\n")
