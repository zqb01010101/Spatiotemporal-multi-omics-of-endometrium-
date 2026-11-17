# =============================================================================
# make_plots.R — Four plot types and a raw-value table
#
# 1. Raw peak profile on its own y-axis
# 2. Raw gene profile on its own y-axis
# 3. Overlay raw peak and gene profiles on a shared y-axis
# 4. Plot diff(peak) and diff(gene) for adjacent stages, with labels such as PO3 to PO5
# =============================================================================

minmax01 <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (r[2] == r[1]) {
    return(rep(0, length(x)))
  }
  (x - r[1]) / (r[2] - r[1])
}

interval_labels <- function(stages) {
  paste0(stages[-length(stages)], "\u2192", stages[-1])
}

.series_frame <- function(values, stages) {
  data.frame(
    stage = factor(stages, levels = stages),
    value = as.numeric(values),
    stringsAsFactors = FALSE
  )
}

plot_raw_peak <- function(values, stages) {
  d <- .series_frame(values, stages)
  ggplot2::ggplot(d, ggplot2::aes(x = stage, y = value)) +
    ggplot2::geom_line(group = 1, color = "#c0392b", linewidth = 0.9) +
    ggplot2::geom_point(color = "#c0392b", size = 2.5) +
    ggplot2::labs(x = "时期", y = "Peak accessibility",
                  title = "原始 Peak 曲线") +
    ggplot2::theme_minimal(base_size = 13)
}

plot_raw_gene <- function(values, stages) {
  d <- .series_frame(values, stages)
  ggplot2::ggplot(d, ggplot2::aes(x = stage, y = value)) +
    ggplot2::geom_line(group = 1, color = "#2471a3", linewidth = 0.9) +
    ggplot2::geom_point(color = "#2471a3", size = 2.5) +
    ggplot2::labs(x = "时期", y = "Gene expression proportion",
                  title = "原始 Gene 曲线") +
    ggplot2::theme_minimal(base_size = 13)
}

plot_overlay <- function(peak_values, gene_values, stages) {
  p <- as.numeric(peak_values)
  g <- as.numeric(gene_values)
  d <- rbind(
    data.frame(stage = factor(stages, levels = stages), value = p, series = "Peak"),
    data.frame(stage = factor(stages, levels = stages), value = g, series = "Gene")
  )
  ggplot2::ggplot(d, ggplot2::aes(x = stage, y = value, color = series, group = series)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_color_manual(values = c(Peak = "#c0392b", Gene = "#2471a3")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "时期", y = "原始值",
                  title = "原始值叠加曲线（Peak accessibility 与 Gene expression proportion，纵轴固定 0-1）",
                  color = "序列") +
    ggplot2::theme_minimal(base_size = 13)
}

plot_deltas <- function(peak_values, gene_values, stages) {
  dp <- diff(as.numeric(peak_values))
  dg <- diff(as.numeric(gene_values))
  labs <- interval_labels(stages)
  d <- rbind(
    data.frame(interval = factor(labs, levels = labs), delta = dp, series = "Peak"),
    data.frame(interval = factor(labs, levels = labs), delta = dg, series = "Gene")
  )
  ggplot2::ggplot(d, ggplot2::aes(x = interval, y = delta, fill = series)) +
    ggplot2::geom_col(position = "dodge", width = 0.65) +
    ggplot2::scale_fill_manual(values = c(Peak = "#c0392b", Gene = "#2471a3")) +
    ggplot2::facet_wrap(~series, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "相邻时期区间", y = "变化量 diff()",
                  title = "相邻时期变化图", fill = "序列") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}

# Raw values and changes between stages
make_value_table <- function(peak_values, gene_values, stages) {
  p <- as.numeric(peak_values)
  g <- as.numeric(gene_values)
  out <- data.frame(
    时期 = stages,
    peak_accessibility = round(p, 6),
    gene_expression = round(g, 6),
    stringsAsFactors = FALSE
  )
  ilabs <- interval_labels(stages)
  dp <- diff(p)
  dg <- diff(g)
  delta_df <- data.frame(
    时期 = ilabs,
    peak_delta = round(dp, 6),
    gene_delta = round(dg, 6),
    stringsAsFactors = FALSE
  )
  list(raw = out, deltas = delta_df)
}
