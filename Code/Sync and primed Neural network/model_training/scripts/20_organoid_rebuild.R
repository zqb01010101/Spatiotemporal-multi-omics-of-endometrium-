#!/usr/bin/env Rscript
# =============================================================================
# 20_organoid_rebuild.R — Rebuild the candidate index and features with organoid avg
#
# Keep stromal pairs with res$mat == 1 and genes present in organoid avg.
# Outputs in organoid_run/:
#   candidates_organoid.csv         Stromal organoid candidate index
#   all_candidates_features.csv     12 stage features for all candidates
#   queue_features.csv              Organoid features for queued candidates
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

res_obj <- readRDS(file.path(root_dir, "26_08_06_res_peak_gene_1Mb.rds"))
res_mat <- res_obj$mat
rm(res_obj); gc()
Y_sig <- readRDS(file.path(root_dir, "26_08_06_Y_sig.rds"))
new_avg <- readRDS(file.path(root_dir,
  "26_09_03_avg_merged_list_change_PO9_11_13_to_organoid.rds"))
inputs <- list(Y_sig = Y_sig, avg = new_avg)
rm(Y_sig, new_avg); gc()

genes_keep <- intersect(colnames(res_mat), rownames(inputs$avg))
cat("organoid 可用的 res 基因数：", length(genes_keep), "\n")

s <- Matrix::summary(res_mat)
s <- s[abs(s$x - 1) < 1e-12, , drop = FALSE]
rn <- rownames(res_mat); cn <- colnames(res_mat)
keep <- cn[s$j] %in% genes_keep
peaks <- rn[s$i[keep]]
genes <- cn[s$j[keep]]
cand <- unique(data.table(
  candidate_id = paste("stroma", peaks, genes, sep = "|"),
  context = "stroma", peak = peaks, gene = genes))
cat("organoid stroma 候选数：", nrow(cand), "\n")

con <- init_db(cfg)
queue_ids <- DBI::dbGetQuery(con,
  "SELECT candidate_id FROM candidates WHERE queue_priority IS NOT NULL AND context = 'stroma'")$candidate_id
DBI::dbDisconnect(con)
q <- cand[candidate_id %in% queue_ids]
cat("organoid 队列候选数：", nrow(q), "\n")

CHUNK <- 250000L
write_features <- function(cand, outfile) {
  if (file.exists(outfile)) file.remove(outfile)
  n_chunks <- ceiling(nrow(cand) / CHUNK)
  for (k in seq_len(n_chunks)) {
    idx <- seq((k - 1L) * CHUNK + 1L, min(k * CHUNK, nrow(cand)))
    vals <- extract_bulk_values(inputs, cfg, "stroma", cand$peak[idx], cand$gene[idx])
    vals[, key_ := NULL]
    out <- cbind(cand[idx, .(candidate_id, context, peak, gene)], vals)
    data.table::fwrite(out, outfile, append = (k > 1L))
    if (k %% 2 == 0) cat(sprintf("  chunk %d/%d 完成\n", k, n_chunks))
    rm(vals, out); gc()
  }
  cat("已写出：", outfile, "\n")
}

data.table::fwrite(cand[, .(candidate_id, context, peak, gene)],
                   file.path(out_dir, "candidates_organoid.csv"))
write_features(cand, file.path(out_dir, "all_candidates_features.csv"))
write_features(q,    file.path(out_dir, "queue_features.csv"))
