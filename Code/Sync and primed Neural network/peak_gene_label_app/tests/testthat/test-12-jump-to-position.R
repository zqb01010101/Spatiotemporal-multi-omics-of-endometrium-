test_that("跳转功能：输入序号直达当前队列第 N 个，超界自动截断", {
  fx <- make_test_fixture("12")
  source_modules()
  build_fixture_db(fx)
  Sys.setenv(PEAK_GENE_FIXTURE_DIR = fx$dir)
  source(file.path(app_dir, "app.R"))

  parse_pos <- function(html) {
    m <- regmatches(html, regexpr("[0-9]+ / [0-9]+", html))
    expect_length(m, 1L)
    as.integer(strsplit(m, " / ", fixed = TRUE)[[1]])
  }

  shiny::testServer(server, {
    session$setInputs(f_context = "all")
    p0 <- parse_pos(as.character(output$status_boxes))
    expect_equal(p0[1], 1L)          # Resume at position 1 when all candidates are unlabeled

    session$setInputs(jump_to = 5)
    session$setInputs(jump_btn = 1)
    p1 <- parse_pos(as.character(output$status_boxes))
    expect_equal(p1[1], 5L)          # Jump to candidate 5
    expect_equal(input$jump_to, 5)   # Keep the input field synchronized with the queue position

    session$setInputs(jump_to = 99999)
    session$setInputs(jump_btn = 1)
    p2 <- parse_pos(as.character(output$status_boxes))
    expect_equal(p2[1], p2[2])       # Clamp beyond-end positions to the last candidate
    expect_equal(p2[2], p0[2])       # The queue size stays the same
  })
})
