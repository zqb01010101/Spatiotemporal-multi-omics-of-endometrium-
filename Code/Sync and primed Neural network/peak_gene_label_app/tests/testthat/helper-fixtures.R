# Shared test helpers to locate the app, create fixtures, and load R modules

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(yaml)
})

app_dir <- Sys.getenv("PEAK_GENE_APP_DIR", unset = normalizePath(getwd()))
source(file.path(app_dir, "tests", "helpers", "make_fixtures.R"))

source_modules <- function() {
  for (f in c("load_inputs.R", "validate_inputs.R", "build_candidates.R",
              "extract_pair_values.R", "database.R", "label_rules.R",
              "queue.R", "make_plots.R")) {
    source(file.path(app_dir, "R", f))
  }
}

make_test_fixture <- function(tag, ...) {
  d <- file.path(tempdir(), paste0("pg_fixture_", tag))
  if (dir.exists(d)) unlink(d, recursive = TRUE)
  info <- make_fixtures(d, ...)
  cfg <- yaml::read_yaml(file.path(d, "config.yml"))
  cfg$app_dir <- d
  list(dir = d, cfg = cfg, info = info)
}

build_fixture_db <- function(fx) {
  inputs <- load_inputs(fx$cfg)
  res <- build_candidate_index(inputs, fx$cfg)
  con <- init_db(fx$cfg)
  replace_candidates(con, res$cand)
  DBI::dbDisconnect(con)
  list(inputs = inputs, cand = res$cand, exclusions = res$exclusions)
}
