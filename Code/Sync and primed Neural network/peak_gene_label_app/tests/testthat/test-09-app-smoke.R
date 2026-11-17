test_that("Shiny App：渲染、盲标、保存落盘与证据显示", {
  fx <- make_test_fixture("09")
  source_modules()
  build_fixture_db(fx)
  Sys.setenv(PEAK_GENE_FIXTURE_DIR = fx$dir)
  source(file.path(app_dir, "app.R"))

  shiny::testServer(server, {
    session$setInputs(
      f_context = "all",
      score_pre_open = 8, score_sync_open = 1,
      score_sync_repress = 0, score_unrelated = 0,
      confidence = 8, notes = "smoke"
    )
    expect_false(is.null(output$status_boxes))
    expect_true(is.null(output$evidence_panel))   # Blinded mode hides evidence by default
    session$setInputs(save_only = 1)
    expect_false(is.null(output$evidence_panel))  # Evidence is visible after saving
  })

  # Verify the annotation was persisted
  con <- init_db(fx$cfg)
  n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations WHERE score_pre_open IS NOT NULL")$n
  expect_equal(n, 1L)
  DBI::dbDisconnect(con)
})
