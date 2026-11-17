test_that("App 重启后恢复进度：标注持久化，队列定位到第一个未标注样本", {
  fx <- make_test_fixture("07")
  source_modules()
  built <- build_fixture_db(fx)
  con <- init_db(fx$cfg)

  # Save the first queued candidate so it appears in the annotation queue
  qd0 <- build_queue(as.data.table(fetch_candidates(con, list(), include_unqueued = FALSE)))
  cid <- qd0$candidate_id[1]
  save_annotation(con, cid,
                  c(pre_open = 8, sync_open = 1, sync_repress = 0, unrelated = 0),
                  overall_confidence = 8, notes = "持久化测试")
  DBI::dbDisconnect(con)

  # Reopen the database to simulate a restart
  con2 <- init_db(fx$cfg)
  ann <- annotation_row(con2, cid)
  expect_false(is.null(ann))
  expect_equal(ann$score_pre_open, 8)
  expect_equal(ann$notes, "持久化测试")

  qd <- build_queue(as.data.table(fetch_candidates(con2, list(), include_unqueued = FALSE)))
  pos <- resume_position(qd)
  expect_equal(pos, 2L)                     # Resume at position 2 after labeling the first candidate
  expect_true(qd$candidate_id[1] == cid)
  expect_equal(sum(!is.na(qd$score_pre_open)), 1L)
  DBI::dbDisconnect(con2)
})
