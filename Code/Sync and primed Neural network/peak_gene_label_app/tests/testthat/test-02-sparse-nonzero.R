test_that("稀疏矩阵非零提取：只取 res 中值为 1 的对，距离正确，距离 0 保留", {
  fx <- make_test_fixture("02")
  source_modules()
  inputs <- load_inputs(fx$cfg)

  cand <- extract_res_candidates(inputs$res_mat, inputs$dist_mat)
  s <- Matrix::summary(inputs$res_mat)
  expected_n <- sum(s$x == 1)

  expect_equal(nrow(cand), expected_n)
  expect_equal(anyDuplicated(cand$key_), 0L)

  ri <- match(cand$peak, rownames(inputs$res_mat))
  ci <- match(cand$gene, colnames(inputs$res_mat))
  expect_true(all(inputs$res_mat[cbind(ri, ci)] == 1))

  dvals <- inputs$dist_mat[cbind(ri, ci)]
  expect_equal(cand$distance, as.numeric(dvals))
  expect_true(any(cand$distance == 0))
})
