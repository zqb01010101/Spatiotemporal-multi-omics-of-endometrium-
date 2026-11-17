test_that("stroma/epithelial 时期列顺序严格且标签顺序正确", {
  fx <- make_test_fixture("04")
  source_modules()

  expect_identical(sec_cols_for(fx$cfg, "stroma"),
                   fx$cfg$stages$sec_stroma_cols)
  expect_identical(sec_cols_for(fx$cfg, "epithelial"),
                   fx$cfg$stages$sec_epi_cols)
  expect_identical(fx$cfg$stages$labels,
                   c("PO3", "PO5", "PO7", "PO9", "PO11", "PO13"))
  expect_length(sec_cols_for(fx$cfg, "stroma"), 6L)
  expect_length(sec_cols_for(fx$cfg, "epithelial"), 6L)
  expect_error(sec_cols_for(fx$cfg, "unknown"))
})
