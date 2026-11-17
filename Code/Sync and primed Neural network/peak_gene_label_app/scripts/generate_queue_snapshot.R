#!/usr/bin/env Rscript
# =============================================================================
# generate_queue_snapshot.R — Export a feature snapshot of the app queue
#
# Usage: Rscript scripts/generate_queue_snapshot.R [organoid|nonorganoid]
#
# Read Y_sig and the selected avg once to extract 12 stage features per queued candidate.
# Write data_cache/queue_features_<version>.csv for the app and training scripts
# to avoid reloading the 2.2 GB matrices.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

script_dir <- normalizePath(dirname(sub("^--file=", "",
                                        commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])))
app_dir   <- normalizePath(dirname(script_dir))
root_dir  <- dirname(app_dir)

for (f in c("load_inputs.R", "extract_pair_values.R", "database.R")) {
  source(file.path(app_dir, "R", f))
}
cfg <- read_config(app_dir)

version <- tolower(commandArgs(trailingOnly = TRUE)[1] %||% "nonorganoid")
if (!version %in% c("organoid", "nonorganoid")) stop("未知版本：", version)

avg_file <- if (version == "organoid") {
  "26_09_03_avg_merged_list_change_PO9_11_13_to_organoid.rds"
} else {
  "26_08_06_avg_list_binary.rds"
}

Y_sig <- readRDS(file.path(root_dir, "26_08_06_Y_sig.rds"))
avg   <- readRDS(file.path(root_dir, avg_file))
inputs <- list(Y_sig = Y_sig, avg = avg)
rm(Y_sig, avg); gc()

con <- init_db(cfg)
queue <- as.data.table(DBI::dbGetQuery(con,
  "SELECT candidate_id, context, peak, gene FROM candidates
   WHERE queue_priority IS NOT NULL"))
DBI::dbDisconnect(con)
cat("队列候选数：", nrow(queue), "；版本：", version, "\n")

# Exclude candidates absent from the selected avg or Y_sig, such as organoid genes missing from the original avg
queue[, keep := peak %in% rownames(inputs$Y_sig) & gene %in% rownames(inputs$avg)]
cat("该版本可用的候选：", nrow(queue[keep == TRUE]), "；排除：", nrow(queue[keep == FALSE]), "\n")
queue <- queue[keep == TRUE][, keep := NULL]

parts <- list()
CHUNK <- 200000L
for (ctx in c("stroma", "epithelial")) {
  sub <- queue[context == ctx]
  n_chunks <- ceiling(nrow(sub) / CHUNK)
  for (k in seq_len(n_chunks)) {
    idx <- seq((k - 1L) * CHUNK + 1L, min(k * CHUNK, nrow(sub)))
    vals <- extract_bulk_values(inputs, cfg, ctx, sub$peak[idx], sub$gene[idx])
    vals[, key_ := NULL]
    parts[[length(parts) + 1L]] <- cbind(
      sub[idx, .(candidate_id, context, peak, gene)], vals)
  }
}
out <- data.table::rbindlist(parts)

out_dir <- file.path(app_dir, "data_cache")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(out_dir, paste0("queue_features_", version, ".csv"))
data.table::fwrite(out, outfile)
cat("已写出：", outfile, "（", nrow(out), " 行）\n")
