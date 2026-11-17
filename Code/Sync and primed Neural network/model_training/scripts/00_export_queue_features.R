#!/usr/bin/env Rscript
# =============================================================================
# 00_export_queue_features.R — Export 12 stage features for the 5,071 candidates in the labeling queue
#
# Score the full labeling queue after training to generate review PDFs
#       for the highest-scoring candidates in three classes and the lowest-confidence candidates.
# Output: model_training/data/queue_features.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
model_dir <- normalizePath(file.path(script_dir, ".."))
app_dir   <- normalizePath(file.path(model_dir, "..", "peak_gene_label_app"))

for (f in c("load_inputs.R", "validate_inputs.R", "build_candidates.R",
            "extract_pair_values.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}

cfg <- read_config(app_dir)
cat("加载输入矩阵（约 1 分钟）...\n")
inputs <- load_inputs(cfg)
validate_inputs(inputs, cfg)

con <- init_db(cfg)
q <- as.data.table(fetch_candidates(con, list(), include_unqueued = FALSE))
DBI::dbDisconnect(con)
cat("标注队列候选数：", nrow(q), "\n")

parts <- list()
for (ctx in c("stroma", "epithelial")) {
  sub <- q[context == ctx]
  if (!nrow(sub)) next
  vals <- extract_bulk_values(inputs, cfg, ctx, sub$peak, sub$gene)
  vals[, key_ := NULL]
  parts[[ctx]] <- cbind(sub[, .(candidate_id, context, peak, gene, distance)], vals)
}
feat <- data.table::rbindlist(parts)

outdir <- file.path(model_dir, "data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, "queue_features.csv")
data.table::fwrite(feat, out)
cat("已写出：", out, "（", nrow(feat), " 行，", ncol(feat), " 列）\n")
