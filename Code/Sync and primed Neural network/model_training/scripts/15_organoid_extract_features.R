#!/usr/bin/env Rscript
# =============================================================================
# 15_organoid_extract_features.R — Extract candidate features with organoid avg
#
# Replace the PO9, PO11, and PO13 expression columns in avg with organoid data.
# Optional first argument: context (stroma, epithelial, or all). Default: all.
# Output directory: model_training/organoid_run/
#   queue_features.csv              12 stage features for queued candidates
#   all_candidates_features.csv     12 stage features for all candidates
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
out_dir    <- file.path(model_dir, "organoid_run")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in c("load_inputs.R", "extract_pair_values.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}
cfg <- read_config(app_dir)

Y_sig <- readRDS(file.path(root_dir, "26_08_06_Y_sig.rds"))
avg   <- readRDS(file.path(root_dir,
  "26_09_03_avg_merged_list_change_PO9_11_13_to_organoid.rds"))
inputs <- list(Y_sig = Y_sig, avg = avg)
rm(Y_sig, avg); gc()

con <- init_db(cfg)
all_cand <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))
queue <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates
        WHERE queue_priority IS NOT NULL"))
DBI::dbDisconnect(con)

args <- commandArgs(trailingOnly = TRUE)
ctx_filter <- if (length(args)) tolower(args[1]) else "all"
if (!ctx_filter %in% c("all", "stroma", "epithelial")) stop("未知 context: ", ctx_filter)
if (ctx_filter != "all") {
  all_cand <- all_cand[context == ctx_filter]
  queue <- queue[context == ctx_filter]
}
cat("context 过滤：", ctx_filter, "\n")
cat("全部候选：", nrow(all_cand), "；队列候选：", nrow(queue), "\n")

CHUNK <- 250000L

write_features <- function(cand, outfile) {
  if (file.exists(outfile)) file.remove(outfile)
  for (ctx in c("stroma", "epithelial")) {
    sub <- cand[context == ctx]
    n_chunks <- ceiling(nrow(sub) / CHUNK)
    for (k in seq_len(n_chunks)) {
      idx <- seq((k - 1L) * CHUNK + 1L, min(k * CHUNK, nrow(sub)))
      vals <- extract_bulk_values(inputs, cfg, ctx, sub$peak[idx], sub$gene[idx])
      vals[, key_ := NULL]
      out <- cbind(sub[idx, .(candidate_id, context, peak, gene)], vals)
      data.table::fwrite(out, outfile, append = (k > 1L || ctx != "stroma"))
      if (k %% 2 == 0) cat(sprintf("[%s] chunk %d/%d 完成\n", ctx, k, n_chunks))
      rm(vals, out); gc()
    }
  }
  cat("已写出：", outfile, "\n")
}

write_features(all_cand, file.path(out_dir, "all_candidates_features.csv"))
write_features(queue,    file.path(out_dir, "queue_features.csv"))
