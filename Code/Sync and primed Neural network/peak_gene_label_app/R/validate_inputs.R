# =============================================================================
# validate_inputs.R — Validate input structures and configuration
# =============================================================================

validate_inputs <- function(inputs, cfg, log_file = NULL) {
  problems <- character()
  add_problem <- function(...) {
    problems <<- c(problems, paste0(..., collapse = " "))
  }

  # ---- Types and dimensions ----
  if (!inherits(inputs$avg, "dgCMatrix")) {
    add_problem("avg 不是 dgCMatrix，实际 class=", paste(class(inputs$avg), collapse = ","))
  }
  if (!inherits(inputs$Y_sig, "matrix")) {
    add_problem("Y_sig 不是普通 matrix，实际 class=", paste(class(inputs$Y_sig), collapse = ","))
  }
  if (!inherits(inputs$res_mat, "dgCMatrix")) {
    add_problem("res$mat 不是 dgCMatrix，实际 class=", paste(class(inputs$res_mat), collapse = ","))
  }
  if (!inherits(inputs$dist_mat, "dgCMatrix")) {
    add_problem("dist$dist_mat 不是 dgCMatrix，实际 class=", paste(class(inputs$dist_mat), collapse = ","))
  }
  if (!identical(dim(inputs$res_mat), dim(inputs$dist_mat))) {
    add_problem("res$mat 与 dist$dist_mat 维度不一致: ",
                paste(dim(inputs$res_mat), collapse = "x"), " vs ",
                paste(dim(inputs$dist_mat), collapse = "x"))
  }
  if (nrow(inputs$Y_sig) != nrow(inputs$res_mat)) {
    add_problem("Y_sig 行数 (", nrow(inputs$Y_sig), ") 与 res$mat 行数 (",
                nrow(inputs$res_mat), ") 不一致")
  }
  if (ncol(inputs$res_mat) != nrow(inputs$avg)) {
    add_problem("res$mat 列数 (", ncol(inputs$res_mat), ") 与 avg 行数 (",
                nrow(inputs$avg), ") 不一致")
  }
  if (length(inputs$common_peaks) == 0 || length(inputs$common_genes) == 0) {
    add_problem("名称对齐后公共 peak 或 gene 数量为 0，无法构建候选")
  }

  # ---- Stage columns ----
  avg_cols <- colnames(inputs$avg)
  y_cols   <- colnames(inputs$Y_sig)
  for (ctx in c("stroma", "epithelial")) {
    sec_cols <- sec_cols_for(cfg, ctx)
    for (c in sec_cols) {
      if (!c %in% avg_cols) add_problem("avg 缺少时期列: ", c)
      if (!c %in% y_cols)  add_problem("Y_sig 缺少时期列: ", c)
    }
  }

  # ---- Configured column names ----
  col_cfg <- cfg$columns
  dfs <- list(
    tri_stroma  = inputs$tri_stroma,
    tri_epi     = inputs$tri_epi,
    sync_stroma = inputs$sync_stroma,
    sync_epi    = inputs$sync_epi
  )
  for (nm in names(dfs)) {
    pk <- col_cfg$peak[[nm]]
    gn <- col_cfg$gene[[nm]]
    if (is.null(pk) || is.na(pk) || pk == "__FILL_ME__") {
      add_problem("config.yml columns.peak.", nm, " 未填写")
    } else if (!pk %in% colnames(dfs[[nm]])) {
      add_problem("config.yml columns.peak.", nm, " = ", pk, " 不在 data.frame 列名中")
    }
    if (is.null(gn) || is.na(gn) || gn == "__FILL_ME__") {
      add_problem("config.yml columns.gene.", nm, " 未填写")
    } else if (!gn %in% colnames(dfs[[nm]])) {
      add_problem("config.yml columns.gene.", nm, " = ", gn, " 不在 data.frame 列名中")
    }
    scol <- col_cfg$S_pre_shape[[nm]] %||% NULL
    fcol <- col_cfg$S_final[[nm]] %||% NULL
    if (is.null(scol) || is.na(scol) || scol == "__FILL_ME__") {
      if (nm %in% c("tri_stroma", "tri_epi")) {
        add_problem("config.yml columns.S_pre_shape.", nm, " 未填写")
      }
    }
    if (is.null(fcol) || is.na(fcol) || fcol == "__FILL_ME__") {
      if (nm %in% c("sync_stroma", "sync_epi")) {
        add_problem("config.yml columns.S_final.", nm, " 未填写")
      }
    }
  }

  # ---- Threshold validation ----
  num_ok <- function(x, label) {
    if (!is.numeric(x) || is.na(x) || x <= 0) {
      add_problem("config.yml 阈值非法（应 > 0）: ", label, " = ", x)
    }
  }
  num_ok(cfg$mandatory$tri_S_pre_shape_threshold, "mandatory.tri_S_pre_shape_threshold")
  num_ok(cfg$mandatory$sync_S_final_threshold, "mandatory.sync_S_final_threshold")
  for (k in c("dynamic_range_peak", "dynamic_range_gene", "min_drop_peak", "min_drop_gene")) {
    num_ok(cfg$sampling$sync_repress[[k]], paste0("sampling.sync_repress.", k))
  }

  if (length(problems)) {
    for (p in problems) log_msg(p, file = log_file, level = "ERROR")
    stop("输入/配置校验失败，共 ", length(problems), " 个问题。")
  }
  log_msg("输入与配置校验通过", file = log_file)
  invisible(TRUE)
}
