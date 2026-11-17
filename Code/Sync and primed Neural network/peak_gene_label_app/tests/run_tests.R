#!/usr/bin/env Rscript
# =============================================================================
# Run all unit tests
#   Run from the peak_gene_label_app directory:
#     Rscript tests/run_tests.R
# =============================================================================

app_dir <- normalizePath(getwd())
stopifnot(file.exists(file.path(app_dir, "config.yml")))
Sys.setenv(PEAK_GENE_APP_DIR = app_dir)

library(testthat)
res <- tryCatch(
  testthat::test_dir(
    file.path(app_dir, "tests", "testthat"),
    reporter = "summary",
    stop_on_failure = TRUE,
    stop_on_warning = FALSE
  ),
  error = function(e) {
    cat("\n测试未全部通过：", conditionMessage(e), "\n")
    quit(save = "no", status = 1)
  }
)
cat("\n全部测试通过。\n")
