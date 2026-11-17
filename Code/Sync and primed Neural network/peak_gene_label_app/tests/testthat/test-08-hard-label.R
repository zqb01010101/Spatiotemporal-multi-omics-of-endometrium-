test_that("hard_label 推导规则：最高分/分差/信心阈值", {
  source_modules()
  cfg <- list(label_rules = list(
    hard_label_min_top = 7,
    hard_label_min_margin = 2,
    hard_label_min_confidence = 6
  ))

  r1 <- derive_hard_label(c(2, 9, 1, 0), 8, cfg)
  expect_equal(r1$hard_label, "sync_open")
  expect_equal(r1$top_score, 9)
  expect_equal(r1$second_score, 2)
  expect_equal(r1$score_margin, 7)

  r2 <- derive_hard_label(c(8, 7, 0, 0), 8, cfg)   # Score margin below threshold
  expect_equal(r2$hard_label, "uncertain")

  r3 <- derive_hard_label(c(8, 2, 0, 0), 5, cfg)   # Confidence below threshold
  expect_equal(r3$hard_label, "uncertain")

  r4 <- derive_hard_label(c(0, 0, 0, 9), 7, cfg)   # Unrelated pair
  expect_equal(r4$hard_label, "unrelated")

  r5 <- derive_hard_label(c(NA, 8, 0, 0), 8, cfg)  # Missing scores
  expect_equal(r5$hard_label, NA_character_)
})
