#!/usr/bin/env Rscript
# =============================================================================
# 01_build_candidate_index.R — Step 2: Build the candidate index and write it to SQLite
#
# Outputs:
#   * peak_gene_label_app/labels/annotations.sqlite (candidates table)
#   * logs/01_build_candidate_index_*.log (alignment, exclusion, and sampling logs)
#   * data_cache/candidate_summary.rds (structured summary)
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

cfg <- read_config(app_dir)
log_file <- open_log(cfg, "01_build_candidate_index")
log_msg("=== 01_build_candidate_index 开始 ===", file = log_file)

# ---- 1) Load read-only inputs once ----
inputs <- load_inputs(cfg, log_file)
validate_inputs(inputs, cfg, log_file)

# ---- 2) Build the candidate index ----
log_msg("开始构建候选索引 ...", file = log_file)
index_res <- build_candidate_index(inputs, cfg, log_file)
cand_dt <- index_res$cand
excl <- index_res$exclusions
if (is.null(excl) || !nrow(excl)) {
  excl <- data.table(reason = character(), n_excluded = integer())
}

summary_tab <- summarize_candidates(cand_dt)
log_msg(sprintf("候选关系总数：%d（stroma %d，epithelial %d）",
                summary_tab$total, summary_tab$total_stroma, summary_tab$total_epithelial),
        file = log_file)
log_msg(sprintf("mandatory 总数：%d（stroma %d，epithelial %d）",
                summary_tab$mandatory_total,
                summary_tab$mandatory_stroma,
                summary_tab$mandatory_epithelial),
        file = log_file)
log_msg(sprintf("同步抑制启发式候选：%d（stroma %d，epithelial %d）",
                summary_tab$heuristic_stroma + summary_tab$heuristic_epithelial,
                summary_tab$heuristic_stroma,
                summary_tab$heuristic_epithelial),
        file = log_file)
log_msg(sprintf("背景样本：%d（stroma %d，epithelial %d）",
                summary_tab$background_stroma + summary_tab$background_epithelial,
                summary_tab$background_stroma,
                summary_tab$background_epithelial),
        file = log_file)
log_msg(sprintf("标注队列总数（含优先级）：%d", summary_tab$queue_total), file = log_file)

if (nrow(excl)) {
  log_msg("==== mandatory 排除汇总 ====", file = log_file)
  for (i in seq_len(nrow(excl))) {
    if (is.na(excl$reason[i]) || excl$n_excluded[i] == 0L) next
    log_msg(sprintf("  [%s/%s] %s：%d", excl$file[i], excl$source[i],
                    excl$reason[i], excl$n_excluded[i]),
            file = log_file, level = "WARN")
  }
}

# ---- 3) Write to SQLite ----
log_msg("写入 SQLite ...", file = log_file)
con <- init_db(cfg)
n_inserted <- replace_candidates(con, cand_dt)
log_msg(sprintf("candidates 表写入行数：%d", n_inserted), file = log_file)

# ---- 4) Check the database ----
db_check <- DBI::dbGetQuery(con, "
  SELECT context,
         COUNT(*) AS n_total,
         SUM(mandatory) AS n_mandatory,
         SUM(CASE WHEN queue_priority IS NOT NULL THEN 1 ELSE 0 END) AS n_queued
  FROM candidates
  GROUP BY context")
log_msg("==== 数据库核对 ====", file = log_file)
print(db_check)

db_mandatory <- DBI::dbGetQuery(con, "
  SELECT COUNT(*) AS n FROM candidates WHERE mandatory = 1")
db_dupes <- DBI::dbGetQuery(con, "
  SELECT COUNT(*) - COUNT(DISTINCT candidate_id) AS n_dup FROM candidates")
log_msg(sprintf("数据库 mandatory 行数：%s；candidate_id 重复数：%s",
                db_mandatory$n, db_dupes$n_dup), file = log_file)
if (db_dupes$n_dup != 0L) {
  stop("candidates 表存在重复 candidate_id！")
}

DBI::dbDisconnect(con)

# ---- 5) Save the summary ----
out <- c(summary_tab, list(
  excluded = excl,
  database = db_path(cfg),
  built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
))
dir.create(file.path(app_dir, "data_cache"), showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path(app_dir, "data_cache", "candidate_summary.rds"))

log_msg("=== 01_build_candidate_index 完成 ===", file = log_file)

# Print the summary to the console
cat("\n==== 候选索引构建完成 ====\n")
cat("候选关系总数      :", summary_tab$total, "\n")
cat("mandatory stroma  :", summary_tab$mandatory_stroma, "\n")
cat("mandatory epithelial:", summary_tab$mandatory_epithelial, "\n")
cat("启发式同步抑制    :", summary_tab$heuristic_stroma + summary_tab$heuristic_epithelial, "\n")
cat("背景样本          :", summary_tab$background_stroma + summary_tab$background_epithelial, "\n")
cat("标注队列总数      :", summary_tab$queue_total, "\n")
cat("数据库            :", db_path(cfg), "\n")
cat("日志              :", log_file, "\n")
if (nrow(excl)) {
  cat("\n被排除关系（按原因）：\n")
  print(excl[n_excluded > 0])
}
