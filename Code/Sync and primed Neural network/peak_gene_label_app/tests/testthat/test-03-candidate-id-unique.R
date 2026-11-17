test_that("candidate_id 唯一，context 拆分正确，优先级合法", {
  fx <- make_test_fixture("03")
  source_modules()
  inputs <- load_inputs(fx$cfg)
  res <- build_candidate_index(inputs, fx$cfg)
  cand <- res$cand

  expect_false(anyDuplicated(cand$candidate_id) > 0)
  n_ones <- sum(Matrix::summary(inputs$res_mat)$x == 1)
  expect_equal(nrow(cand), 2L * n_ones)
  expect_equal(nrow(cand[context == "stroma"]), n_ones)
  expect_equal(nrow(cand[context == "epithelial"]), n_ones)

  expect_true(all(cand$mandatory %in% c(0L, 1L)))
  expect_true(all(is.na(cand$queue_priority) | cand$queue_priority %in% c(1L, 2L, 3L)))
  # Mandatory rows must have priority 1
  expect_true(all(cand[mandatory == 1L]$queue_priority == 1L))
  # Heuristic candidates per context must stay within the configured cap
  expect_lte(nrow(cand[context == "stroma" & source_sync_repress_heuristic == 1L]),
             fx$cfg$sampling$sync_repress$max_per_context)
  expect_lte(nrow(cand[context == "epithelial" & source_sync_repress_heuristic == 1L]),
             fx$cfg$sampling$sync_repress$max_per_context)
  # Background candidates per context must stay within the configured count
  expect_lte(nrow(cand[context == "stroma" & queue_priority == 3L]),
             fx$cfg$sampling$n_background_per_context)
  expect_lte(nrow(cand[context == "epithelial" & queue_priority == 3L]),
             fx$cfg$sampling$n_background_per_context)
  # The exclusion log contains both expected reasons
  reasons <- res$exclusions$reason
  expect_true("gene 不在 avg" %in% reasons)
  expect_true("peak 不在 Y_sig" %in% reasons)
})
