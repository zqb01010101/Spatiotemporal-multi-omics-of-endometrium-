test_that("叠加图使用原始 peak/gene 值（不进行 0-1 标准化）", {
  source_modules()
  stages <- c("PO3", "PO5", "PO7", "PO9", "PO11", "PO13")

  pk <- c(0.02, 0.9, 0.95, 0.5, 0.1, 0.3)
  gn <- c(0.6, 0.4, 0.2, 0.1, 0.05, 0.0)
  p <- plot_overlay(pk, gn, stages)
  d <- ggplot2::ggplot_build(p)$data[[1]]
  expect_setequal(round(d$y, 6), round(c(pk, gn), 6))
  expect_equal(p$coordinates$limits$y, c(0, 1))   # Fixed y-axis range of [0, 1]

  # Constant profiles retain their original values
  p2 <- plot_overlay(rep(5, 6), rep(3, 6), stages)
  d2 <- ggplot2::ggplot_build(p2)$data[[1]]
  expect_true(all(d2$y > 0))
  expect_setequal(unique(d2$y), c(5, 3))
})
