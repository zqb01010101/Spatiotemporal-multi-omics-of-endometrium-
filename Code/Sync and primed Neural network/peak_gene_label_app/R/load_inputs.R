# =============================================================================
# load_inputs.R — Read configuration, resolve paths, and align eight RDS inputs by name
#
# Input handling:
#   * Load input matrices once at startup and keep them read-only.
#   * Keep dgCMatrix objects sparse.
#   * Align rows and columns by name.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# Internal pair key with a nonprinting separator to avoid ambiguous IDs
make_pair_key <- function(peak, gene) paste(peak, gene, sep = "\u0001")

# -----------------------------------------------------------------------------
# Locate the app directory
#   - Prefer the PEAK_GENE_APP_DIR environment variable.
#   - Otherwise infer the app directory from Rscript --file as the parent of scripts/.
#   - Fall back to the working directory, as with shiny::runApp('.').
# -----------------------------------------------------------------------------
get_app_dir <- function() {
  env_dir <- Sys.getenv("PEAK_GENE_APP_DIR", unset = "")
  if (nzchar(env_dir)) {
    return(normalizePath(env_dir, mustWork = TRUE))
  }
  args0 <- commandArgs(FALSE)
  f <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  if (length(f) && nzchar(f[1])) {
    return(normalizePath(dirname(dirname(f[1])), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

read_config <- function(app_dir = get_app_dir()) {
  p <- file.path(app_dir, "config.yml")
  if (!file.exists(p)) {
    stop("找不到 config.yml：", p)
  }
  cfg <- yaml::read_yaml(p)
  cfg$app_dir <- app_dir
  cfg
}

# Resolve relative paths from app_dir
resolve_path <- function(cfg, p) {
  if (is.null(p) || length(p) == 0L || is.na(p)) {
    return(NA_character_)
  }
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) {
    return(p)
  }
  normalizePath(file.path(cfg$app_dir, p), mustWork = FALSE)
}

# Log to the console and optionally to a file
log_msg <- function(..., file = NULL, level = "INFO") {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                "] [", level, "] ", paste0(..., collapse = " "))
  cat(msg, "\n")
  if (!is.null(file)) {
    cat(msg, "\n", file = file, append = TRUE)
  }
  invisible(msg)
}

open_log <- function(cfg, prefix) {
  log_dir <- resolve_path(cfg, cfg$paths$logs)
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(log_dir, sprintf("%s_%s.log", prefix, format(Sys.time(), "%Y%m%d_%H%M%S")))
  file.create(f, showWarnings = FALSE)
  f
}

# Get the six stage columns for the context in config.yml order
sec_cols_for <- function(cfg, context) {
  if (context == "stroma") {
    return(cfg$stages$sec_stroma_cols)
  }
  if (context == "epithelial") {
    return(cfg$stages$sec_epi_cols)
  }
  stop("未知 context: ", context)
}

# Reuse the original object when names already align
subset_by_name <- function(m, rows = NULL, cols = NULL) {
  need_r <- !is.null(rows) && !identical(rownames(m), rows)
  need_c <- !is.null(cols) && !identical(colnames(m), cols)
  if (need_r && need_c) {
    return(m[rows, cols, drop = FALSE])
  }
  if (need_r) {
    return(m[rows, , drop = FALSE])
  }
  if (need_c) {
    return(m[, cols, drop = FALSE])
  }
  m
}

# -----------------------------------------------------------------------------
# Read eight RDS inputs and align matrices by name
# -----------------------------------------------------------------------------
load_inputs <- function(cfg, log_file = NULL) {
  df <- cfg$data_files
  read_one <- function(key) {
    p <- resolve_path(cfg, df[[key]])
    log_msg("读取 ", key, " <- ", p, file = log_file)
    obj <- readRDS(p)
    log_msg("完成 ", key, " class=", paste(class(obj), collapse = ","),
            " dim=", paste(dim(obj), collapse = "x"), file = log_file)
    obj
  }

  inputs <- list(
    avg          = read_one("avg"),
    Y_sig        = read_one("Y_sig"),
    dist         = read_one("dist"),
    res          = read_one("res"),
    tri_stroma   = read_one("tri_stroma"),
    tri_epi      = read_one("tri_epi"),
    sync_stroma  = read_one("sync_stroma"),
    sync_epi     = read_one("sync_epi")
  )

  align_inputs(inputs, cfg, log_file)
}

align_inputs <- function(inputs, cfg, log_file = NULL) {
  res_mat  <- inputs$res$mat
  dist_mat <- inputs$dist$dist_mat
  Y_sig    <- inputs$Y_sig
  avg      <- inputs$avg

  if (is.null(res_mat) || is.null(dist_mat)) {
    stop("res$mat / dist$dist_mat 缺失")
  }

  common_peaks <- Reduce(intersect, list(rownames(Y_sig), rownames(res_mat), rownames(dist_mat)))
  common_genes <- Reduce(intersect, list(rownames(avg), colnames(res_mat), colnames(dist_mat)))

  miss_peak <- setdiff(rownames(res_mat), common_peaks)
  miss_gene <- setdiff(colnames(res_mat), common_genes)
  log_msg("名称对齐：公共 peak = ", length(common_peaks),
          "，公共 gene = ", length(common_genes), file = log_file)
  if (length(miss_peak)) {
    log_msg("被排除的 peak（不在 Y_sig/dist）：", length(miss_peak),
            "，示例：", paste(utils::head(miss_peak, 5), collapse = " | "),
            file = log_file, level = "WARN")
  }
  if (length(miss_gene)) {
    log_msg("被排除的 gene（不在 avg/dist）：", length(miss_gene),
            "，示例：", paste(utils::head(miss_gene, 5), collapse = " | "),
            file = log_file, level = "WARN")
  }

  inputs$res_mat  <- subset_by_name(res_mat, common_peaks, common_genes)
  inputs$dist_mat <- subset_by_name(dist_mat, common_peaks, common_genes)
  inputs$Y_sig    <- subset_by_name(Y_sig, common_peaks)
  inputs$avg      <- subset_by_name(avg, common_genes)

  inputs$common_peaks <- common_peaks
  inputs$common_genes <- common_genes
  inputs$meta <- list(
    stage_labels    = cfg$stages$labels,
    sec_stroma_cols = cfg$stages$sec_stroma_cols,
    sec_epi_cols    = cfg$stages$sec_epi_cols
  )
  inputs
}
