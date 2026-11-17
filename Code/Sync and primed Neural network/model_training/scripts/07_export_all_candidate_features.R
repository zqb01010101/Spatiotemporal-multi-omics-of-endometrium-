#!/usr/bin/env Rscript
# =============================================================================
# 07_export_all_candidate_features.R — Export 12 stage features for all candidates
#
# Output: model_training/data/all_candidates_features.csv
#       （candidate_id, context, peak, gene + peak_PO3..PO13 + gene_PO3..PO13）
# Read avg and Y_sig only, and process in chunks to limit memory use.
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
CHUNK <- 250000L

avg   <- readRDS(resolve_path(cfg, cfg$data_files$avg))
Y_sig <- readRDS(resolve_path(cfg, cfg$data_files$Y_sig))
inputs <- list(Y_sig = Y_sig, avg = avg)
rm(avg, Y_sig); gc()

con <- init_db(cfg)
cand <- as.data.table(DBI::dbGetQuery(
  con, "SELECT candidate_id, context, peak, gene FROM candidates"))
DBI::dbDisconnect(con)
cat("全部候选：", nrow(cand), "\n")

outfile <- file.path(model_dir, "data", "all_candidates_features.csv")
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
