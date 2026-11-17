# =============================================================================
# extract_pair_values.R — Retrieve six-stage profiles for an individual peak–gene pair
#
# Read six peak and six gene values on demand by peak ID, gene ID, and context.
# Keep the input matrices read-only and sparse.
# =============================================================================

# Return list(peak = six values, gene = six values, stages = labels) for a single pair
extract_pair_values <- function(inputs, cfg, context, peak, gene) {
  sec_cols <- sec_cols_for(cfg, context)
  stages   <- cfg$stages$labels

  if (!peak %in% rownames(inputs$Y_sig)) {
    stop("peak 不在 Y_sig 中: ", peak)
  }
  if (!gene %in% rownames(inputs$avg)) {
    stop("gene 不在 avg 中: ", gene)
  }

  peak_vals <- as.numeric(inputs$Y_sig[peak, sec_cols, drop = FALSE])
  gene_vals <- as.numeric(inputs$avg[gene, sec_cols, drop = FALSE])

  if (length(peak_vals) != length(sec_cols) || length(gene_vals) != length(sec_cols)) {
    stop("时期值提取失败: ", context, " | ", peak, " | ", gene)
  }
  names(peak_vals) <- stages
  names(gene_vals) <- stages
  list(peak = peak_vals, gene = gene_vals, stages = stages)
}

# Extract pairs in bulk for heuristic sampling and export, returning a data.table:
# key, peak_PO3..peak_PO13, gene_PO3..gene_PO13
extract_bulk_values <- function(inputs, cfg, context, peaks, genes) {
  sec_cols <- sec_cols_for(cfg, context)
  stages   <- cfg$stages$labels
  stopifnot(length(peaks) == length(genes))
  n <- length(peaks)
  if (n == 0L) {
    return(data.table::data.table())
  }

  pk_idx <- match(peaks, rownames(inputs$Y_sig))
  gn_idx <- match(genes, rownames(inputs$avg))
  if (anyNA(pk_idx) || anyNA(gn_idx)) {
    stop("批量提取遇到未知 peak/gene 名称")
  }

  pv <- inputs$Y_sig[pk_idx, sec_cols, drop = FALSE]
  gv <- inputs$avg[gn_idx, sec_cols, drop = FALSE]
  if (inherits(pv, "dgeMatrix") || inherits(pv, "dgCMatrix")) pv <- as.matrix(pv)
  if (inherits(gv, "dgeMatrix") || inherits(gv, "dgCMatrix")) gv <- as.matrix(gv)

  out <- data.table::data.table(key_ = make_pair_key(peaks, genes))
  for (k in seq_along(stages)) {
    out[[paste0("peak_", stages[k])]] <- pv[, k]
    out[[paste0("gene_", stages[k])]] <- gv[, k]
  }
  out
}
