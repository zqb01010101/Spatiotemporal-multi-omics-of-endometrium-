test_that("min-max 标准化：常量序列全为 0，正常序列归一化正确", {
  source_modules()
  expect_equal(minmax01(rep(5, 6)), rep(0, 6))
  expect_equal(minmax01(c(1, 2, 3)), c(0, 0.5, 1))
  expect_equal(minmax01(c(0, 0, 0)), rep(0, 3))
})
