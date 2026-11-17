#!/usr/bin/env Rscript
# =============================================================================
# 00_inspect_inputs.R — Step 1: Inspect eight RDS inputs
#
# Features:
#   1. Read each RDS and report class, dimensions, names, and example row and column names.
#   2. Check for the dist_mat and mat list elements.
#   3. Check that avg and Y_sig contain all six stages for stroma and epithelium.
#   4. Check peak and gene names across Y_sig, res, dist, and avg.
#   5. Detect peak and gene columns in the four statistics tables.
#   6. Check how well peak and gene names in the four tables overlap the matrix names.
#   7. Generate config.yml. Leave __FILL_ME__ placeholders and stop if column detection is ambiguous.
#   8. Write reports to logs/ and data_cache/inspect_summary.rds.
#
# Early exit rule:
#   If peak name overlap for any two inputs is below 1%,
#   write the report and exit with status 2 so the input files can be corrected.
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(yaml)
})

options(warn = 1)

CRITICAL_OVERLAP_RATIO <- 0.01   # Minimum overlap fraction

# -----------------------------------------------------------------------------
# The script lives in <app>/scripts/. Inputs live one level above the app directory.
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  args0 <- commandArgs(FALSE)
  f <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  if (length(f) && nzchar(f[1])) normalizePath(dirname(f[1])) else normalizePath(getwd())
}
script_dir <- get_script_dir()
app_dir    <- normalizePath(dirname(script_dir))
root_dir   <- normalizePath(dirname(app_dir))

# -----------------------------------------------------------------------------
# Logs and reports
# -----------------------------------------------------------------------------
log_dir <- file.path(app_dir, "logs")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
report_file <- file.path(log_dir, sprintf("00_inspect_report_%s.txt", stamp))
log_fh <- file(report_file, open = "wt", encoding = "UTF-8")

report <- function(..., level = "INFO") {
  msg <- paste0("[", level, "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  cat(msg, "\n", file = log_fh)
  flush(log_fh)
}

report_block <- function(title, lines) {
  report(sprintf("==== %s ====", title))
  for (ln in lines) {
    cat(ln, "\n")
    cat(ln, "\n", file = log_fh)
  }
  flush(log_fh)
}

report_r_object <- function(obj, txt) {
  # Write a text summary of each object to the report
  out <- capture.output(print(obj))
  report_block(txt, out)
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# -----------------------------------------------------------------------------
# Define the read-only RDS inputs
# -----------------------------------------------------------------------------
input_files <- c(
  avg          = "26_08_06_avg_list_binary.rds",
  Y_sig        = "26_08_06_Y_sig.rds",
  dist         = "26_08_06_PeakGene_TSSdist_1Mb.rds",
  res          = "26_08_06_res_peak_gene_1Mb.rds",
  tri_stroma   = "26_08_06_tri_Model_result_stroma_1Mb.rds",
  tri_epi      = "26_08_06_tri_Model_result_epi_1Mb.rds",
  sync_stroma  = "26_06_20_Sync_peak_gene_stroma.rds",
  sync_epi     = "26_06_20_Sync_peak_gene_epi.rds"
)

# Stage columns and labels in required order
sec_stroma_cols <- c(
  "sec_PO3_stroma", "sec_PO5_stroma", "sec_PO7_stroma",
  "sec_PO9_stroma", "sec_PO11_stroma", "sec_PO13_stroma"
)
sec_epi_cols <- c(
  "sec_PO3_epithelial", "sec_PO5_epithelial", "sec_PO7_epithelial",
  "sec_PO9_epithelial", "sec_PO11_epithelial", "sec_PO13_epithelial"
)
stage_labels <- c("PO3", "PO5", "PO7", "PO9", "PO11", "PO13")

# Detect candidate column names
PEAK_CANDIDATES  <- c("peak", "Peak", "peak_name", "peak_id", "PeakID", "region")
GENE_CANDIDATES  <- c("Gene", "gene", "gene_name", "gene_id", "symbol")

report("脚本目录: ", script_dir)
report("应用目录: ", app_dir)
report("项目根目录: ", root_dir)
report("报告文件: ", report_file)

# -----------------------------------------------------------------------------
# Read all inputs
# -----------------------------------------------------------------------------
report("================ 开始读取输入文件 ================")
objs <- list()
for (nm in names(input_files)) {
  p <- file.path(root_dir, input_files[[nm]])
  if (!file.exists(p)) {
    report(sprintf("缺失文件: %s (%s)", nm, p), level = "ERROR")
    stop("缺少输入文件: ", p)
  }
  sz_mb <- round(file.info(p)$size / 1024^2, 1)
  report(sprintf("[读取] %-12s  (%s, %.1f MB)", nm, input_files[[nm]], sz_mb))
  objs[[nm]] <- readRDS(p)
  report(sprintf("[完成] %-12s  class=%s", nm, paste(class(objs[[nm]]), collapse = ",")))
  gc()
}

# -----------------------------------------------------------------------------
# Describe each object's structure
# -----------------------------------------------------------------------------
report("================ 对象结构概览 ================")

describe_obj <- function(obj, nm) {
  cls <- paste(class(obj), collapse = ",")
  if (!is.null(dim(obj))) {
    d <- paste(dim(obj), collapse = " x ")
  } else {
    d <- as.character(length(obj))
  }
  report(sprintf("%-12s class=%s  dim/length=%s", nm, cls, d))

  if (is.list(obj) && !is.data.frame(obj)) {
    report(sprintf("  names: %s", paste(names(obj), collapse = ", ")))
  }
  if (!is.null(dim(obj))) {
    rn <- rownames(obj)
    cn <- if (length(dim(obj)) >= 2) colnames(obj) else NULL
    if (!is.null(rn)) {
      report(sprintf("  rownames(%d): %s", length(rn), paste(utils::head(rn, 5), collapse = " | ")))
      dup <- sum(duplicated(rn))
      if (dup > 0) report(sprintf("  !! rownames 存在 %d 个重复", dup), level = "WARN")
    } else {
      report("  rownames: NULL")
    }
    if (!is.null(cn)) {
      report(sprintf("  colnames(%d): %s", length(cn), paste(utils::head(cn, 8), collapse = " | ")))
      dup <- sum(duplicated(cn))
      if (dup > 0) report(sprintf("  !! colnames 存在 %d 个重复", dup), level = "WARN")
    }
  }
  if (inherits(obj, "dgCMatrix")) {
    report(sprintf("  nnz(nonzero) = %s", format(length(obj@x), big.mark = ",")))
  }
}

for (nm in names(input_files)) describe_obj(objs[[nm]], nm)

report("================ 1) list 对象内部检查 ================")
for (nm in c("dist", "res")) {
  obj <- objs[[nm]]
  report(sprintf("%s 是否为 list: %s", nm, is.list(obj) && !is.data.frame(obj)))
  if (is.list(obj) && !is.data.frame(obj)) {
    report(sprintf("  $dist_mat 存在: %s", !is.null(obj$dist_mat)))
    report(sprintf("  $mat      存在: %s", !is.null(obj$mat)))
    if (!is.null(obj$dist_mat)) {
      describe_obj(obj$dist_mat, paste0(nm, "$dist_mat"))
    }
    if (!is.null(obj$mat)) {
      describe_obj(obj$mat, paste0(nm, "$mat"))
    }
  }
}

report("================ 2) 时期列检查 ================")
avg_cols <- colnames(objs$avg)
Y_cols   <- colnames(objs$Y_sig)
check_sec_cols <- function(cols, where) {
  miss <- setdiff(cols, where)
  if (length(miss) == 0) {
    report(sprintf("  %s 的六个时期列全部存在于 %s", paste(cols, collapse = ","), "对象"))
    TRUE
  } else {
    report(sprintf("  缺失时期列（%s）: %s", "对象", paste(miss, collapse = ", ")), level = "ERROR")
    FALSE
  }
}
ok_stroma_avg <- check_sec_cols(sec_stroma_cols, avg_cols)
ok_epi_avg    <- check_sec_cols(sec_epi_cols, avg_cols)
ok_stroma_Y   <- check_sec_cols(sec_stroma_cols, Y_cols)
ok_epi_Y      <- check_sec_cols(sec_epi_cols, Y_cols)
if (!all(ok_stroma_avg, ok_epi_avg, ok_stroma_Y, ok_epi_Y)) {
  stop("时期列缺失，无法继续。")
}
report(sprintf("  时期列顺序确认: stroma=%s", paste(sec_stroma_cols, collapse = " > ")))
report(sprintf("                   epithelial=%s", paste(sec_epi_cols, collapse = " > ")))

# -----------------------------------------------------------------------------
# 3) Check named matrix objects before detecting data-frame columns
# -----------------------------------------------------------------------------
report("================ 3) peak / gene 名称一致性检查 ================")

pk_Y    <- rownames(objs$Y_sig)
pk_res  <- rownames(objs$res$mat)
pk_dist <- rownames(objs$dist$dist_mat)
gn_avg  <- rownames(objs$avg)
gn_res  <- colnames(objs$res$mat)
gn_dist <- colnames(objs$dist$dist_mat)

peak_sets <- list(
  "Y_sig 行(peak)"       = pk_Y,
  "res$mat 行(peak)"     = pk_res,
  "dist$dist_mat 行(peak)" = pk_dist
)
gene_sets <- list(
  "avg 行(gene)"         = gn_avg,
  "res$mat 列(gene)"     = gn_res,
  "dist$dist_mat 列(gene)" = gn_dist
)

overlap_stats <- function(A, B) {
  nA <- length(A); nB <- length(B)
  inter <- length(intersect(A, B))
  ratio <- if (min(nA, nB) > 0) inter / min(nA, nB) else 0
  c(nA = nA, nB = nB, intersect = inter, ratio = round(ratio, 4))
}

report("  --- 同域名称重叠表（peak 与 peak、gene 与 gene；ratio = intersect / min(nA,nB)） ---")
matrix_overlap <- list()
stop_flag <- FALSE

show_mismatch_examples <- function(A, B, labelA, labelB) {
  report(sprintf("  !! %s 与 %s 几乎不重叠，示例（前 5 个互差）:", labelA, labelB), level = "ERROR")
  report(sprintf("    仅在 %s: %s", labelA, paste(utils::head(setdiff(A, B), 5), collapse = " | ")))
  report(sprintf("    仅在 %s: %s", labelB, paste(utils::head(setdiff(B, A), 5), collapse = " | ")))
}

check_within_domain <- function(sets, domain, ref_list) {
  for (i in seq_along(sets)) {
    for (j in seq_len(i - 1L)) {
      key <- paste(names(sets)[i], "<->", names(sets)[j])
      st <- overlap_stats(sets[[i]], sets[[j]])
      ref_list[[key]] <- st
      flag <- if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) "  <<< CRITICAL" else ""
      report(sprintf("    %-48s nA=%-9s nB=%-9s inter=%-9s ratio=%.4f%s",
                     key, st[["nA"]], st[["nB"]], st[["intersect"]], st[["ratio"]], flag))
      if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) {
        stop_flag <<- TRUE
        show_mismatch_examples(sets[[i]], sets[[j]], names(sets)[i], names(sets)[j])
      }
    }
  }
  ref_list
}

matrix_overlap <- check_within_domain(peak_sets, "peak", matrix_overlap)
matrix_overlap <- check_within_domain(gene_sets, "gene", matrix_overlap)

# Pairs whose names should match exactly
expect_same <- list(
  c("res$mat 行(peak)", "dist$dist_mat 行(peak)"),
  c("res$mat 列(gene)", "dist$dist_mat 列(gene)")
)
for (pr in expect_same) {
  A <- if (pr[1] %in% names(peak_sets)) peak_sets[[pr[1]]] else gene_sets[[pr[1]]]
  B <- if (pr[2] %in% names(peak_sets)) peak_sets[[pr[2]]] else gene_sets[[pr[2]]]
  if (!identical(A, B)) {
    report(sprintf("  !! %s 与 %s 不完全一致（应完全相同）: 相同=%d, 互差=%d/%d",
                   pr[1], pr[2], length(intersect(A, B)), length(setdiff(A, B)), length(setdiff(B, A))),
           level = "WARN")
  }
}

if (stop_flag) {
  report("================ 关键停点 ================", level = "ERROR")
  report("存在 peak/gene 名称几乎不重叠的文件对（重叠率 < 1%）。", level = "ERROR")
  report("按用户要求：停止全部后续工作，等待用户修复数据后重新运行本脚本。", level = "ERROR")
  close(log_fh)
  quit(save = "no", status = 2)
}

# -----------------------------------------------------------------------------
# 4) Detect peak and gene columns in the four data frames
# -----------------------------------------------------------------------------
report("================ 4) 统计结果 data.frame 的 peak/gene 列自动探测 ================")

detect_one <- function(nms, candidates, file_label, role_label) {
  hits <- nms[tolower(nms) %in% tolower(candidates)]
  hits_unique <- unique(hits)
  if (length(hits_unique) == 1) {
    report(sprintf("  [%s] %s 列 -> %s", file_label, role_label, hits_unique))
    return(hits_unique)
  }
  report(sprintf("  [%s] %s 列无法唯一判定（候选命中 %d 个）", file_label, role_label, length(hits_unique)), level = "ERROR")
  report(sprintf("  [%s] 全部列名: %s", file_label, paste(nms, collapse = ", ")), level = "ERROR")
  report_block(sprintf("[%s] 前 10 行", file_label), capture.output(print(utils::head(objs[[file_label]], 10))))
  return(NA_character_)
}

df_names <- c("tri_stroma", "tri_epi", "sync_stroma", "sync_epi")
detected <- list()
for (dfn in df_names) {
  nms <- colnames(objs[[dfn]])
  detected[[dfn]] <- list(
    peak = detect_one(nms, PEAK_CANDIDATES, dfn, "peak"),
    gene = detect_one(nms, GENE_CANDIDATES, dfn, "gene"),
    S_pre_shape = if ("S_pre_shape" %in% nms) "S_pre_shape" else NA_character_,
    S_final     = if ("S_final" %in% nms) "S_final" else NA_character_
  )
  report(sprintf("  [%s] S_pre_shape 列: %s ; S_final 列: %s",
                 dfn, detected[[dfn]]$S_pre_shape %||% "未找到",
                 detected[[dfn]]$S_final %||% "未找到"))
}

ambiguous <- sapply(detected, function(x) any(is.na(c(x$peak, x$gene))))
if (any(ambiguous)) {
  report("================ 列名探测不明确 ================", level = "ERROR")
  report("存在多个候选或无法判断的 peak/gene 列。已在 config.yml 中留下 __FILL_ME__ 字段。", level = "ERROR")
  report("请人工填写 config.yml 后重新运行本脚本；候选索引构建已停止。", level = "ERROR")
  # Write config.yml with placeholders, then exit
} else {
  report("全部 4 个 data.frame 的 peak/gene 列均唯一确定。")
}

# -----------------------------------------------------------------------------
# 5) Check peak and gene name overlap between data frames and matrices
# -----------------------------------------------------------------------------
report("================ 5) data.frame 与矩阵名称重叠检查 ================")
df_overlap <- list()
for (dfn in df_names) {
  pk <- detected[[dfn]]$peak
  gn <- detected[[dfn]]$gene
  if (is.na(pk) || is.na(gn)) next
  df_peaks <- as.character(objs[[dfn]][[pk]])
  df_genes <- as.character(objs[[dfn]][[gn]])
  report(sprintf("  [%s] 行数=%d, peak 列唯一值=%d, gene 列唯一值=%d",
                 dfn, nrow(objs[[dfn]]), length(unique(df_peaks)), length(unique(df_genes))))
  for (ref_name in c("res$mat 行(peak)", "Y_sig 行(peak)")) {
    ref <- peak_sets[[ref_name]]
    st <- overlap_stats(df_peaks, ref)
    key <- sprintf("%s[peak] <-> %s", dfn, ref_name)
    df_overlap[[key]] <- st
    flag <- if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) "  <<< CRITICAL" else ""
    report(sprintf("    %-52s nA=%-8s nB=%-8s inter=%-8s ratio=%.4f%s",
                   key, st[["nA"]], st[["nB"]], st[["intersect"]], st[["ratio"]], flag))
    if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) {
      stop_flag <- TRUE
      show_mismatch_examples(df_peaks, ref, sprintf("%s[peak]", dfn), ref_name)
    }
  }
  for (ref_name in c("res$mat 列(gene)", "avg 行(gene)")) {
    ref <- gene_sets[[ref_name]]
    st <- overlap_stats(df_genes, ref)
    key <- sprintf("%s[gene] <-> %s", dfn, ref_name)
    df_overlap[[key]] <- st
    flag <- if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) "  <<< CRITICAL" else ""
    report(sprintf("    %-52s nA=%-8s nB=%-8s inter=%-8s ratio=%.4f%s",
                   key, st[["nA"]], st[["nB"]], st[["intersect"]], st[["ratio"]], flag))
    if (st[["ratio"]] < CRITICAL_OVERLAP_RATIO && st[["nA"]] > 0 && st[["nB"]] > 0) {
      stop_flag <- TRUE
      show_mismatch_examples(df_genes, ref, sprintf("%s[gene]", dfn), ref_name)
    }
  }
}

if (stop_flag) {
  report("================ 关键停点 ================", level = "ERROR")
  report("data.frame 与矩阵之间的 peak/gene 名称几乎不重叠（重叠率 < 1%）。", level = "ERROR")
  report("按用户要求：停止全部后续工作，等待用户修复数据后重新运行本脚本。", level = "ERROR")
  close(log_fh)
  quit(save = "no", status = 2)
}

# -----------------------------------------------------------------------------
# 6) Write config.yml
# -----------------------------------------------------------------------------
report("================ 6) 生成 config.yml ================")

fill_or_placeholder <- function(x) if (is.na(x)) "__FILL_ME__" else x

cfg <- list(
  paths = list(
    root_dir  = "..",
    app_dir   = ".",
    database  = "labels/annotations.sqlite",
    data_cache = "data_cache",
    logs       = "logs",
    labels     = "labels"
  ),
  data_files = list(
    avg         = file.path("..", input_files[["avg"]]),
    Y_sig       = file.path("..", input_files[["Y_sig"]]),
    dist        = file.path("..", input_files[["dist"]]),
    res         = file.path("..", input_files[["res"]]),
    tri_stroma  = file.path("..", input_files[["tri_stroma"]]),
    tri_epi     = file.path("..", input_files[["tri_epi"]]),
    sync_stroma = file.path("..", input_files[["sync_stroma"]]),
    sync_epi    = file.path("..", input_files[["sync_epi"]])
  ),
  stages = list(
    labels = stage_labels,
    sec_stroma_cols = sec_stroma_cols,
    sec_epi_cols    = sec_epi_cols
  ),
  columns = list(
    peak = list(
      tri_stroma  = fill_or_placeholder(detected$tri_stroma$peak),
      tri_epi     = fill_or_placeholder(detected$tri_epi$peak),
      sync_stroma = fill_or_placeholder(detected$sync_stroma$peak),
      sync_epi    = fill_or_placeholder(detected$sync_epi$peak)
    ),
    gene = list(
      tri_stroma  = fill_or_placeholder(detected$tri_stroma$gene),
      tri_epi     = fill_or_placeholder(detected$tri_epi$gene),
      sync_stroma = fill_or_placeholder(detected$sync_stroma$gene),
      sync_epi    = fill_or_placeholder(detected$sync_epi$gene)
    ),
    S_pre_shape = list(
      tri_stroma = fill_or_placeholder(detected$tri_stroma$S_pre_shape),
      tri_epi    = fill_or_placeholder(detected$tri_epi$S_pre_shape)
    ),
    S_final = list(
      sync_stroma = fill_or_placeholder(detected$sync_stroma$S_final),
      sync_epi    = fill_or_placeholder(detected$sync_epi$S_final)
    )
  ),
  mandatory = list(
    tri_S_pre_shape_threshold = 0.2,
    sync_S_final_threshold    = 0.7
  ),
  sampling = list(
    random_seed = 20260806,
    n_background_per_context = 1000,
    sync_repress = list(
      sample_size            = 20000,
      max_per_context        = 500,
      dynamic_range_peak     = 0.05,
      dynamic_range_gene     = 0.05,
      min_drop_peak          = 0.05,
      min_drop_gene          = 0.05,
      max_window_gap         = 1
    )
  ),
  queue_priorities = list(
    mandatory             = 1,
    sync_repress_heuristic = 2,
    random_background      = 3
  ),
  label_rules = list(
    hard_label_min_top        = 7,
    hard_label_min_margin     = 2,
    hard_label_min_confidence = 6
  ),
  app = list(
    blind_default = TRUE,
    annotator_default = "annotator1",
    app_version = "1.0.0"
  ),
  qc = list(
    n_consistency_repeats = 20
  )
)

config_path <- file.path(app_dir, "config.yml")
write_yaml(cfg, config_path)
report("config.yml 已写入: ", config_path)

# -----------------------------------------------------------------------------
# 7) Save the inspection summary for downstream scripts
# -----------------------------------------------------------------------------
summary_out <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
  root_dir     = root_dir,
  app_dir      = app_dir,
  report_file  = report_file,
  object_meta  = lapply(objs, function(o) {
    list(
      class = class(o),
      dim   = if (!is.null(dim(o))) dim(o) else length(o),
      names = if (is.list(o) && !is.data.frame(o)) names(o) else NULL
    )
  }),
  sec_stroma_cols = sec_stroma_cols,
  sec_epi_cols    = sec_epi_cols,
  sec_cols_ok     = list(
    avg_stroma = ok_stroma_avg, avg_epi = ok_epi_avg,
    Y_stroma   = ok_stroma_Y,   Y_epi    = ok_epi_Y
  ),
  name_sets = list(
    peak_Y_sig = pk_Y, peak_res = pk_res, peak_dist = pk_dist,
    gene_avg = gn_avg, gene_res = gn_res, gene_dist = gn_dist
  ),
  matrix_overlap = matrix_overlap,
  df_overlap    = df_overlap,
  detected_columns = detected,
  ambiguous = ambiguous
)
saveRDS(summary_out, file.path(app_dir, "data_cache", "inspect_summary.rds"))
report("检查摘要已写入: ", file.path(app_dir, "data_cache", "inspect_summary.rds"))

report("================ 阶段一检查完成 ================")
report(sprintf("报告文件: %s", report_file))
close(log_fh)

if (any(ambiguous)) {
  quit(save = "no", status = 3)
}
