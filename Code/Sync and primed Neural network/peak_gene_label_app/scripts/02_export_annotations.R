#!/usr/bin/env Rscript
# =============================================================================
# 02_export_annotations.R — Export annotations and training data
#
# Outputs in peak_gene_label_app/labels/:
#   labeled_training_data.rds / .csv   Complete training samples by context, peak, and gene
#   annotations_export.rds / .csv     All annotations with derived hard-label fields
#
# Training data include:
#   Six raw peak values, six raw gene values, changes between stages, scaled values in [0, 1],
#   four scores, confidence, hard_label, mandatory and source flags, score evidence, and notes.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

detect_app_dir <- function() {
  args0 <- commandArgs(FALSE)
  f <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  if (length(f) && nzchar(f[1])) {
    normalizePath(dirname(dirname(f[1])))
  } else {
    normalizePath(getwd())
  }
}
app_dir <- detect_app_dir()
for (f in c("load_inputs.R", "validate_inputs.R", "build_candidates.R",
            "extract_pair_values.R", "database.R", "label_rules.R",
            "queue.R", "make_plots.R")) {
  source(file.path(app_dir, "R", f))
}

cfg_app_dir <- Sys.getenv("PEAK_GENE_CONFIG_APP_DIR", unset = app_dir)
cfg <- read_config(cfg_app_dir)
log_file <- open_log(cfg, "02_export_annotations")
log_msg("=== 02_export_annotations 开始 ===", file = log_file)

stages <- cfg$stages$labels
delta_pairs <- paste0(stages[-length(stages)], "_", stages[-1])

# ---- 1) Read annotations ----
con <- init_db(cfg)
ann <- fetch_all_annotations(con)
DBI::dbDisconnect(con)
log_msg("标注/复核/跳过记录数：", nrow(ann), file = log_file)

if (!nrow(ann)) {
  log_msg("数据库中没有标注记录，跳过训练数据导出。", file = log_file, level = "WARN")
  cat("没有标注记录。请先在 Shiny 中完成标注。\n")
  quit(save = "no", status = 0)
}

# ---- 2) Read profile matrices for on-demand lookup ----
inputs <- load_inputs(cfg)
validate_inputs(inputs, cfg)

# ---- 3) Build the training data ----
scored <- ann[!is.na(ann$score_pre_open) &
              !is.na(ann$score_sync_open) &
              !is.na(ann$score_sync_repress) &
              !is.na(ann$score_unrelated), ]
log_msg("四维评分完整的样本数：", nrow(scored), file = log_file)

if (nrow(scored)) {
  rows_list <- vector("list", nrow(scored))
  for (i in seq_len(nrow(scored))) {
    r <- scored[i, ]
    v <- extract_pair_values(inputs, cfg, r$context, r$peak, r$gene)
    p <- as.numeric(v$peak)
    g <- as.numeric(v$gene)
    dp <- diff(p)
    dg <- diff(g)
    pn <- minmax01(p)
    gn <- minmax01(g)
    hl <- derive_hard_label(
      c(r$score_pre_open, r$score_sync_open, r$score_sync_repress, r$score_unrelated),
      r$overall_confidence, cfg
    )
    row <- c(
      candidate_id = r$candidate_id,
      context = r$context,
      peak = r$peak,
      gene = r$gene,
      distance = r$distance,
      setNames(p, paste0("peak_", stages)),
      setNames(g, paste0("gene_", stages)),
      setNames(dp, paste0("peak_delta_", delta_pairs)),
      setNames(dg, paste0("gene_delta_", delta_pairs)),
      setNames(pn, paste0("peak_norm_", stages)),
      setNames(gn, paste0("gene_norm_", stages)),
      score_pre_open = r$score_pre_open,
      score_sync_open = r$score_sync_open,
      score_sync_repress = r$score_sync_repress,
      score_unrelated = r$score_unrelated,
      overall_confidence = r$overall_confidence,
      hard_label = hl$hard_label,
      top_score = hl$top_score,
      second_score = hl$second_score,
      score_margin = hl$score_margin,
      mandatory = r$mandatory,
      source_tri_pre = r$source_tri_pre,
      source_sync_open = r$source_sync_open,
      source_sync_repress_heuristic = r$source_sync_repress_heuristic,
      S_pre_shape = r$S_pre_shape,
      S_final = r$S_final,
      notes = r$notes %||% "",
      updated_at = r$updated_at
    )
    rows_list[[i]] <- as.data.frame(t(row), stringsAsFactors = FALSE)
  }
  train <- data.table::rbindlist(rows_list, fill = TRUE)
  for (num_col in c("distance", paste0("peak_", stages), paste0("gene_", stages),
                    paste0("peak_delta_", delta_pairs), paste0("gene_delta_", delta_pairs),
                    paste0("peak_norm_", stages), paste0("gene_norm_", stages),
                    "score_pre_open", "score_sync_open", "score_sync_repress",
                    "score_unrelated", "overall_confidence", "top_score",
                    "second_score", "score_margin")) {
    if (num_col %in% names(train)) {
      data.table::set(train, j = num_col, value = suppressWarnings(as.numeric(train[[num_col]])))
    }
  }
  log_msg("训练样本行数：", nrow(train), file = log_file)
} else {
  train <- data.table()
}

# ---- 4) Export all annotations with derived hard labels ----
hl_all <- lapply(seq_len(nrow(ann)), function(i) {
  r <- ann[i, ]
  derive_hard_label(
    c(r$score_pre_open, r$score_sync_open, r$score_sync_repress, r$score_unrelated),
    r$overall_confidence, cfg
  )
})
ann_export <- data.table::as.data.table(ann)
ann_export[, hard_label := vapply(hl_all, `[[`, character(1), "hard_label")]
ann_export[, top_score := vapply(hl_all, `[[`, numeric(1), "top_score")]
ann_export[, second_score := vapply(hl_all, `[[`, numeric(1), "second_score")]
ann_export[, score_margin := vapply(hl_all, `[[`, numeric(1), "score_margin")]

# ---- 5) Write files ----
labels_dir <- resolve_path(cfg, cfg$paths$labels)
dir.create(labels_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(train, file.path(labels_dir, "labeled_training_data.rds"))
saveRDS(ann_export, file.path(labels_dir, "annotations_export.rds"))
data.table::fwrite(train, file.path(labels_dir, "labeled_training_data.csv"))
data.table::fwrite(ann_export, file.path(labels_dir, "annotations_export.csv"))

log_msg("已写出：", file = log_file)
log_msg("  labels/labeled_training_data.rds / .csv", file = log_file)
log_msg("  labels/annotations_export.rds / .csv", file = log_file)
log_msg("=== 02_export_annotations 完成 ===", file = log_file)

cat("\n==== 标注导出完成 ====\n")
cat("训练样本数 :", nrow(train), "\n")
cat("标注记录数 :", nrow(ann_export), "\n")
cat("输出目录   :", labels_dir, "\n")
