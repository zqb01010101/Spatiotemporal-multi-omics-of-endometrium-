test_that("SQLite 保存与更新：重复保存更新原记录，不产生重复行", {
  fx <- make_test_fixture("06")
  source_modules()
  built <- build_fixture_db(fx)
  con <- init_db(fx$cfg)

  cid <- built$cand$candidate_id[1]
  save_annotation(con, cid,
                  c(pre_open = 8, sync_open = 2, sync_repress = 1, unrelated = 0),
                  overall_confidence = 7, notes = "第一次")
  expect_equal(nrow(get_annotation(con, cid)), 1L)

  save_annotation(con, cid,
                  c(pre_open = 1, sync_open = 9, sync_repress = 0, unrelated = 0),
                  overall_confidence = 8, notes = "更新")
  n_rows <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations")$n
  expect_equal(n_rows, 1L)
  row2 <- get_annotation(con, cid)
  expect_equal(row2$score_sync_open, 9)
  expect_equal(row2$score_pre_open, 1)
  expect_equal(row2$notes, "更新")
  expect_equal(row2$overall_confidence, 8)
  expect_false(is.na(row2$updated_at))

  DBI::dbDisconnect(con)
})
