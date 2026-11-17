#!/usr/bin/env Rscript
# =============================================================================
# 05_generate_review_batch.R — Generate another batch of review candidates
#
# Batch 1: synchronous repression (method v1)
#   Compare the peak profile with the inverted expression profile, 1 - gene.
#   Score candidates by Pearson cor(peak, 1 - gene). Require a positive correlation,
#   sufficient variation in both profiles, and no prior label. Keep the top 50 within 1 Mb.
#
# Batch 2: low confidence
#   Take the 50 lowest-confidence predictions in the labeling queue,
#   excluding labeled candidates and overlaps with batch 1.
#
# Outputs:
#   data/review_batch_v1.csv           Combined batch (100 entries with batch_id and batch_score)
#   data/review_sync_repress_top50.csv Top 50 synchronous repression candidates
#   data/review_low_confidence50.csv   50 lowest-confidence candidates
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "",
                                        commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
model_dir <- normalizePath(file.path(script_dir, ".."))
app_dir   <- normalizePath(file.path(model_dir, "..", "peak_gene_label_app"))

for (f in c("load_inputs.R", "extract_pair_values.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}

cfg <- read_config(app_dir)
stages <- cfg$stages$labels
MIN_RANGE   <- 0.05   # Minimum dynamic range for peak and gene profiles
MIN_CORR    <- 0.30   # Minimum synchronous repression score
N_SYNC      <- 50L
N_LOWCONF   <- 50L
CHUNK       <- 250000L

# Load avg and Y_sig only, skipping the large res and dist matrices to save about 4.4 GB
avg   <- readRDS(resolve_path(cfg, cfg$data_files$avg))
Y_sig <- readRDS(resolve_path(cfg, cfg$data_files$Y_sig))
inputs <- list(Y_sig = Y_sig, avg = avg)
rm(avg, Y_sig)
gc()
cat("已读取 avg 与 Y_sig（不加载 res/dist 大矩阵）\n")

con <- init_db(cfg)
cand <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))
annotated <- DBI::dbGetQuery(
  con, "SELECT candidate_id FROM annotations
        WHERE score_pre_open IS NOT NULL OR score_sync_open IS NOT NULL
           OR score_sync_repress IS NOT NULL OR score_unrelated IS NOT NULL")
DBI::dbDisconnect(con)
cat("全 1Mb 候选：", nrow(cand), "；已标注：", nrow(annotated), "\n")

# -----------------------------------------------------------------------------
# 1) Score synchronous repression with cor(peak, 1 - gene)
# -----------------------------------------------------------------------------
parts <- list()
for (ctx in c("stroma", "epithelial")) {
  sub <- cand[context == ctx]
  chunk_rows <- list()
  n_chunks <- ceiling(nrow(sub) / CHUNK)
  for (k in seq_len(n_chunks)) {
    idx <- seq((k - 1L) * CHUNK + 1L, min(k * CHUNK, nrow(sub)))
    vals <- extract_bulk_values(inputs, cfg, ctx, sub$peak[idx], sub$gene[idx])
    pcols <- paste0("peak_", stages)
    gcols <- paste0("gene_", stages)
    P <- as.matrix(vals[, pcols, with = FALSE])
    G <- as.matrix(vals[, gcols, with = FALSE])
    Q <- 1 - G

    pbar <- rowMeans(P)
    qbar <- rowMeans(Q)
    num  <- rowSums((P - pbar) * (Q - qbar))
    den  <- sqrt(rowSums((P - pbar)^2) * rowSums((Q - qbar)^2))
    corr <- ifelse(den > 1e-9, num / den, NA_real_)
    rng_p <- apply(P, 1L, max) - apply(P, 1L, min)
    rng_g <- apply(G, 1L, max) - apply(G, 1L, min)

    chunk_rows[[k]] <- data.table(
      candidate_id = sub$candidate_id[idx],
      context      = sub$context[idx],
      peak         = sub$peak[idx],
      gene         = sub$gene[idx],
      sync_repress_score = corr,
      range_peak   = rng_p,
      range_gene   = rng_g
    )
    if (k %% 2 == 0) cat(sprintf("[%s] chunk %d/%d 完成\n", ctx, k, n_chunks))
    rm(vals, P, G, Q); gc()
  }
  parts[[ctx]] <- data.table::rbindlist(chunk_rows)
}
scores <- data.table::rbindlist(parts)
rm(parts, chunk_rows); gc()

scores <- scores[!candidate_id %in% annotated$candidate_id]
scores <- scores[!is.na(sync_repress_score) &
                 sync_repress_score > MIN_CORR &
                 range_peak >= MIN_RANGE & range_gene >= MIN_RANGE]
data.table::setorder(scores, -sync_repress_score)
top_sr <- head(scores, N_SYNC)
cat("同步抑制 top", nrow(top_sr), "（分数范围 ",
    sprintf("%.3f ~ %.3f", min(top_sr$sync_repress_score), max(top_sr$sync_repress_score)), "）\n")

# -----------------------------------------------------------------------------
# 2) Select the lowest-confidence model predictions
# -----------------------------------------------------------------------------
qp <- data.table::fread(file.path(model_dir, "data", "queue_predictions.csv"))
qp <- qp[!candidate_id %in% annotated$candidate_id]
data.table::setorder(qp, confidence)
top_lc <- head(qp[, .(candidate_id, context, peak, gene, confidence, pred_class)], N_LOWCONF)
cat("低信心 top", nrow(top_lc), "（confidence 范围 ",
    sprintf("%.3f ~ %.3f", min(top_lc$confidence), max(top_lc$confidence)), "）\n")

# -----------------------------------------------------------------------------
# 3) Merge batches, giving synchronous repression priority and replenishing the low-confidence batch
# -----------------------------------------------------------------------------
sr_ids <- top_sr$candidate_id
top_lc <- top_lc[!candidate_id %in% sr_ids]
if (nrow(top_lc) > N_LOWCONF) top_lc <- head(top_lc, N_LOWCONF)

batch_sr <- data.table(
  candidate_id = top_sr$candidate_id,
  context      = top_sr$context,
  peak         = top_sr$peak,
  gene         = top_sr$gene,
  batch_id     = "review_v1_sync_repress",
  batch_score  = top_sr$sync_repress_score,
  note         = "sync_repress: cor(peak, 1-gene)"
)
batch_lc <- data.table(
  candidate_id = top_lc$candidate_id,
  context      = top_lc$context,
  peak         = top_lc$peak,
  gene         = top_lc$gene,
  batch_id     = "review_v1_low_confidence",
  batch_score  = top_lc$confidence,
  note         = "lowest model confidence (top1 prob)"
)
batch <- data.table::rbindlist(list(batch_sr, batch_lc))
stopifnot(!anyDuplicated(batch$candidate_id))
cat("合并复核批次：", nrow(batch), " 条（同步抑制 ", nrow(batch_sr),
    "，低信心 ", nrow(batch_lc), "）\n")

data.table::fwrite(batch, file.path(model_dir, "data", "review_batch_v1.csv"))
data.table::fwrite(top_sr, file.path(model_dir, "data", "review_sync_repress_top50.csv"))
data.table::fwrite(top_lc, file.path(model_dir, "data", "review_low_confidence50.csv"))
cat("已写出：\n  data/review_batch_v1.csv\n  data/review_sync_repress_top50.csv\n  data/review_low_confidence50.csv\n")
