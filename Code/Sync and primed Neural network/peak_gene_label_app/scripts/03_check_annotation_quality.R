#!/usr/bin/env Rscript
# =============================================================================
# 03_check_annotation_quality.R — Check annotation quality
#
# Usage:
#   Rscript scripts/03_check_annotation_quality.R
#   Rscript scripts/03_check_annotation_quality.R --insert-repeats 20
#
# --insert-repeats N: Insert N repeat candidates for agreement checks.
#   Link repeats to the original with candidates.is_repeat and repeat_of for review in the app.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
n_repeats <- NULL
if (length(args)) {
  idx <- grep("^--insert-repeats$", args)
  if (length(idx)) {
    n_repeats <- as.integer(args[idx + 1L])
  }
}

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
log_file <- open_log(cfg, "03_annotation_quality")
log_msg("=== 03_check_annotation_quality 开始 ===", file = log_file)

stages <- cfg$stages$labels
con <- init_db(cfg)

rep <- function(...) {
  msg <- paste0(..., collapse = " ")
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
  invisible(msg)
}

# -----------------------------------------------------------------------------
# 0) Optional: insert repeat candidates
# -----------------------------------------------------------------------------
if (!is.null(n_repeats)) {
  rep(sprintf("插入 %d 个一致性重复样本 ...", n_repeats))
  q <- DBI::dbGetQuery(con, "
    SELECT c.candidate_id, c.context, c.peak, c.gene, c.distance,
           c.S_pre_shape, c.S_final
    FROM candidates c
    JOIN annotations a ON a.candidate_id = c.candidate_id
    WHERE a.score_pre_open IS NOT NULL
      AND c.is_repeat = 0
    ORDER BY RANDOM()
    LIMIT ?", params = list(n_repeats))
  if (!nrow(q)) {
    rep("没有已标注样本可供插入重复样本，跳过。", level = "WARN")
  } else {
    now <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
    repeat_rows <- lapply(seq_len(nrow(q)), function(i) {
      r <- q[i, ]
      data.frame(
        candidate_id = paste0(r$candidate_id, "#REPEAT", i),
        context = r$context,
        peak = r$peak,
        gene = r$gene,
        distance = r$distance,
        mandatory = 1L,
        source_tri_pre = 0L,
        source_sync_open = 0L,
        source_sync_repress_heuristic = 0L,
        S_pre_shape = r$S_pre_shape,
        S_final = r$S_final,
        queue_priority = 1L,
        created_at = now,
        is_repeat = 1L,
        repeat_of = r$candidate_id,
        stringsAsFactors = FALSE
      )
    })
    rdt <- data.table::rbindlist(repeat_rows)
    DBI::dbBegin(con)
    DBI::dbAppendTable(con, "candidates", rdt)
    DBI::dbCommit(con)
    rep(sprintf("已插入 %d 个重复样本，示例：%s -> %s",
                nrow(rdt), rdt$candidate_id[1], rdt$repeat_of[1]))
  }
}

# -----------------------------------------------------------------------------
# 1) Load annotations
# -----------------------------------------------------------------------------
ann <- fetch_all_annotations(con)
scored <- ann[!is.na(ann$score_pre_open) &
              !is.na(ann$score_sync_open) &
              !is.na(ann$score_sync_repress) &
              !is.na(ann$score_unrelated), ]

rep("==== 1) 基础数量 ====")
rep("标注/复核/跳过记录总数 :", nrow(ann))
rep("四维评分完整样本数     :", nrow(scored))
rep("skipped 数             :", sum(!is.na(ann$review_status) & ann$review_status == "skipped"))
rep("review 数              :", sum(!is.na(ann$review_status) & ann$review_status == "review"))

if (!nrow(scored)) {
  rep("没有评分完整的标注，无法做质量统计。")
  DBI::dbDisconnect(con)
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# 2) Invalid or missing scores
# -----------------------------------------------------------------------------
score_cols <- c("score_pre_open", "score_sync_open", "score_sync_repress", "score_unrelated")
illegal <- sum(apply(scored[, score_cols], 1L, function(x) any(x < 0 | x > 10 | is.na(x))))
rep("==== 2) 缺失值 / 非法评分 ====")
rep("非法或缺失评分行数 :", illegal)
rep("信心缺失行数       :", sum(is.na(scored$overall_confidence)))

# -----------------------------------------------------------------------------
# 3) Class distribution
# -----------------------------------------------------------------------------
rep("==== 3) 四维评分分布（0-10） ====")
for (sc in score_cols) {
  v <- scored[[sc]]
  rep(sprintf("  %-20s mean=%.2f sd=%.2f min=%.0f max=%.0f",
              sc, mean(v, na.rm = TRUE), sd(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}

hl <- vapply(seq_len(nrow(scored)), function(i) {
  derive_hard_label(
    c(scored$score_pre_open[i], scored$score_sync_open[i],
      scored$score_sync_repress[i], scored$score_unrelated[i]),
    scored$overall_confidence[i], cfg
  )$hard_label
}, character(1))
scored$hard_label <- hl

rep("==== 4) hard_label 数量 ====")
hl_tab <- sort(table(scored$hard_label, useNA = "ifany"), decreasing = TRUE)
for (nm in names(hl_tab)) rep(sprintf("  %-12s %d", nm, hl_tab[[nm]]))
uncertain_ratio <- mean(scored$hard_label == "uncertain")
rep(sprintf("uncertain 比例 : %.3f", uncertain_ratio))

rep("==== 5) stroma / epithelial 标签数量 ====")
ct_tab <- table(scored$context)
for (nm in names(ct_tab)) rep(sprintf("  %-12s %d", nm, ct_tab[[nm]]))

rep("==== 6) mandatory 与背景样本数量 ====")
rep("  mandatory          :", sum(scored$mandatory == 1L))
rep("  sync_repress_heuristic:", sum(scored$source_sync_repress_heuristic == 1L))
rep("  背景（非强制非启发式）:", sum(scored$mandatory == 0L & scored$source_sync_repress_heuristic == 0L))
rep("  一致性重复样本      :", sum(scored$is_repeat == 1L))

# -----------------------------------------------------------------------------
# 7) Mean peak and gene profiles by label
# -----------------------------------------------------------------------------
inputs <- load_inputs(cfg)
validate_inputs(inputs, cfg)
rep("==== 7) 每个 hard_label 的平均 peak/gene 曲线 ====")
for (lab in names(hl_tab)) {
  sub <- scored[scored$hard_label == lab, ]
  if (!nrow(sub)) next
  pk_mat <- matrix(NA_real_, nrow(sub), length(stages))
  gn_mat <- matrix(NA_real_, nrow(sub), length(stages))
  for (i in seq_len(nrow(sub))) {
    v <- extract_pair_values(inputs, cfg, sub$context[i], sub$peak[i], sub$gene[i])
    pk_mat[i, ] <- v$peak
    gn_mat[i, ] <- v$gene
  }
  rep(sprintf("  [%s] n=%d  peak均值: %s", lab, nrow(sub),
              paste(round(colMeans(pk_mat, na.rm = TRUE), 4), collapse = " ")))
  rep(sprintf("  [%s] n=%d  gene均值: %s", lab, nrow(sub),
              paste(round(colMeans(gn_mat, na.rm = TRUE), 4), collapse = " ")))
}

# -----------------------------------------------------------------------------
# 8) Human score distributions by candidate source
# -----------------------------------------------------------------------------
rep("==== 8) 不同来源候选的人工评分分布 ====")
scored$source_group <- ifelse(scored$mandatory == 1L, "mandatory",
                       ifelse(scored$source_sync_repress_heuristic == 1L,
                              "sync_repress_heuristic",
                              ifelse(scored$is_repeat == 1L, "repeat", "random_background")))
for (sg in unique(scored$source_group)) {
  sub <- scored[scored$source_group == sg, ]
  rep(sprintf("  [%s] n=%d  预开放=%.2f 同步开放=%.2f 同步抑制=%.2f 无关=%.2f",
              sg, nrow(sub),
              mean(sub$score_pre_open, na.rm = TRUE),
              mean(sub$score_sync_open, na.rm = TRUE),
              mean(sub$score_sync_repress, na.rm = TRUE),
              mean(sub$score_unrelated, na.rm = TRUE)))
}

# -----------------------------------------------------------------------------
# 9) Annotation order over time
# -----------------------------------------------------------------------------
rep("==== 9) 标注时间顺序 ====")
ord <- scored[order(scored$updated_at), ]
rep("  最早更新时间 :", ord$updated_at[1])
rep("  最晚更新时间 :", ord$updated_at[nrow(ord)])
rep("  时间戳唯一值 :", length(unique(scored$updated_at)))

# -----------------------------------------------------------------------------
# 10) Duplicate candidates and class imbalance
# -----------------------------------------------------------------------------
rep("==== 10) 重复与失衡检查 ====")
dup <- sum(duplicated(scored$candidate_id))
rep("  重复 candidate_id 行数 :", dup)
lab_counts <- as.numeric(hl_tab)
if (length(lab_counts) >= 2) {
  imbalance <- max(lab_counts) / min(lab_counts)
  rep(sprintf("  标签数量最大/最小比 : %.2f", imbalance))
  if (imbalance > 10) rep("  !! 标签类别严重失衡（>10 倍），建议补充少数类样本", level = "WARN")
}

# -----------------------------------------------------------------------------
# 11) Repeat annotation agreement where labels are available
# -----------------------------------------------------------------------------
rep("==== 11) 一致性重复样本检查 ====")
reps <- ann[!is.na(ann$repeat_of) & !is.na(ann$score_pre_open), ]
if (!nrow(reps)) {
  rep("  暂无已标注的重复样本。可用 --insert-repeats N 插入后再在 App 中标注。")
} else {
  agree <- 0L
  for (i in seq_len(nrow(reps))) {
    orig <- scored[scored$candidate_id == reps$repeat_of[i], ]
    if (!nrow(orig)) next
    h1 <- derive_hard_label(
      c(reps$score_pre_open[i], reps$score_sync_open[i],
        reps$score_sync_repress[i], reps$score_unrelated[i]),
      reps$overall_confidence[i], cfg)$hard_label
    h2 <- derive_hard_label(
      c(orig$score_pre_open[1], orig$score_sync_open[1],
        orig$score_sync_repress[1], orig$score_unrelated[1]),
      orig$overall_confidence[1], cfg)$hard_label
    if (identical(h1, h2)) agree <- agree + 1L
  }
  rep(sprintf("  已标注重复样本 %d 条，hard_label 一致 %d 条", nrow(reps), agree))
}

DBI::dbDisconnect(con)
rep("=== 03_check_annotation_quality 完成 ===")
cat("\n质量检查报告已写入：", log_file, "\n")
