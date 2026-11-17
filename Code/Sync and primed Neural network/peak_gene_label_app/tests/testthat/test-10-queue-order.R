test_that("队列顺序：mandatory 优先、组内稳定、可恢复", {
  fx <- make_test_fixture("10")
  source_modules()
  built <- build_fixture_db(fx)
  queued <- built$cand[!is.na(queue_priority)]

  q1 <- build_queue(as.data.table(queued))
  q2 <- build_queue(as.data.table(queued))
  expect_identical(q1$candidate_id, q2$candidate_id)   # Stable queue order
  expect_true(all(diff(q1$queue_priority) >= 0))        # Monotonic priority order
  expect_true(all(q1$queue_priority[seq_len(sum(q1$mandatory))] == 1L))

  # Mandatory candidates are shuffled within their priority group
  mand_ids <- q1$candidate_id[q1$mandatory == 1L]
  expect_false(identical(mand_ids, sort(mand_ids)))
})
