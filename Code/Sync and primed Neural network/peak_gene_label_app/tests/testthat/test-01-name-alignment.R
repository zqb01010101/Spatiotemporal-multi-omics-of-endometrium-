test_that("输入名称对齐：按名称取公共 peak/gene，缺失 gene 被排除", {
  fx <- make_test_fixture("01")
  source_modules()

  inputs <- load_inputs(fx$cfg)
  out <- capture.output(ok <- validate_inputs(inputs, fx$cfg))
  expect_true(ok)

  expect_true(inherits(inputs$avg, "dgCMatrix"))
  expect_true(inherits(inputs$Y_sig, "matrix"))
  expect_true(inherits(inputs$res_mat, "dgCMatrix"))
  expect_true(inherits(inputs$dist_mat, "dgCMatrix"))

  expect_setequal(rownames(inputs$Y_sig), fx$info$peaks)
  expect_setequal(rownames(inputs$res_mat), fx$info$peaks)
  expect_setequal(rownames(inputs$dist_mat), fx$info$peaks)
  expect_setequal(colnames(inputs$res_mat), fx$info$genes)
  expect_setequal(colnames(inputs$dist_mat), fx$info$genes)
  expect_setequal(rownames(inputs$avg), fx$info$genes)

  expect_false("GENE_NOT_IN_AVG" %in% inputs$common_genes)
  expect_true(all(fx$cfg$stages$sec_stroma_cols %in% colnames(inputs$Y_sig)))
  expect_true(all(fx$cfg$stages$sec_epi_cols %in% colnames(inputs$Y_sig)))
  expect_true(all(fx$cfg$stages$sec_stroma_cols %in% colnames(inputs$avg)))
  expect_true(all(fx$cfg$stages$sec_epi_cols %in% colnames(inputs$avg)))
})
