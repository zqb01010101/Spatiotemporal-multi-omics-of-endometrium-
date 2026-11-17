# =============================================================================
# database.R — SQLite candidates and annotations tables
#
# Persistence rules:
#   * Write in transactions and persist on save.
#   * Update the existing record on repeated saves (UPSERT).
#   * Persist annotations so work can resume after an unexpected shutdown.
# =============================================================================

db_path <- function(cfg) {
  resolve_path(cfg, cfg$paths$database)
}

init_db <- function(cfg) {
  p <- db_path(cfg)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(con, "PRAGMA synchronous = NORMAL")
  DBI::dbExecute(con, "PRAGMA foreign_keys = ON")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS candidates (
      candidate_id                   TEXT PRIMARY KEY,
      context                        TEXT NOT NULL,
      peak                           TEXT NOT NULL,
      gene                           TEXT NOT NULL,
      distance                       REAL,
      mandatory                      INTEGER NOT NULL DEFAULT 0,
      source_tri_pre                 INTEGER NOT NULL DEFAULT 0,
      source_sync_open               INTEGER NOT NULL DEFAULT 0,
      source_sync_repress_heuristic  INTEGER NOT NULL DEFAULT 0,
      S_pre_shape                    REAL,
      S_final                        REAL,
      queue_priority                 INTEGER,
      created_at                     TEXT,
      is_repeat                      INTEGER NOT NULL DEFAULT 0,
      repeat_of                      TEXT,
      batch_id                       TEXT,
      batch_score                    REAL
    )")

  DBI::dbExecute(con, "
    CREATE TABLE IF NOT EXISTS annotations (
      candidate_id       TEXT PRIMARY KEY,
      score_pre_open     REAL,
      score_sync_open    REAL,
      score_sync_repress REAL,
      score_unrelated    REAL,
      overall_confidence REAL,
      review_status      TEXT,
      notes              TEXT,
      annotator          TEXT,
      created_at         TEXT,
      updated_at         TEXT,
      app_version        TEXT
    )")

  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_context ON candidates(context)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_priority ON candidates(queue_priority)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_peak ON candidates(peak)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_gene ON candidates(gene)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_mandatory ON candidates(mandatory)")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ann_status ON annotations(review_status)")

  # Migrate older databases by adding batch columns
  ensure_column(con, "candidates", "batch_id", "batch_id TEXT")
  ensure_column(con, "candidates", "batch_score", "batch_score REAL")
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cand_batch ON candidates(batch_id)")
  con
}

# Use ALTER TABLE only for missing columns
ensure_column <- function(con, table, column, ddl) {
  cols <- DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", table))$name
  if (!column %in% cols) {
    DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s", table, ddl))
  }
  invisible(TRUE)
}

# Rebuild candidates by clearing the table and inserting rows in transactional batches
replace_candidates <- function(con, cand_dt, chunk_size = 50000L) {
  DBI::dbExecute(con, "DELETE FROM candidates")
  total <- nrow(cand_dt)
  n_done <- 0L
  for (start in seq(1L, total, by = chunk_size)) {
    end <- min(start + chunk_size - 1L, total)
    DBI::dbBegin(con)
    DBI::dbAppendTable(con, "candidates", cand_dt[start:end])
    DBI::dbCommit(con)
    n_done <- n_done + (end - start + 1L)
  }
  n_done
}

get_annotation <- function(con, candidate_id) {
  DBI::dbGetQuery(
    con,
    "SELECT * FROM annotations WHERE candidate_id = ?",
    params = list(candidate_id)
  )
}

annotation_row <- function(con, candidate_id) {
  res <- get_annotation(con, candidate_id)
  if (nrow(res) == 0L) NULL else as.list(res[1, ])
}

# Save or update an annotation immediately (UPSERT)
save_annotation <- function(con, candidate_id, scores, overall_confidence,
                            review_status = NULL, notes = NULL,
                            annotator = NULL, app_version = NULL) {
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
  old <- annotation_row(con, candidate_id)
  created_at <- if (!is.null(old) && !is.na(old$created_at)) old$created_at else now
  if (is.null(notes) && !is.null(old)) notes <- old$notes
  if (is.null(annotator) && !is.null(old)) annotator <- old$annotator

  DBI::dbBegin(con)
  DBI::dbExecute(
    con,
    "INSERT INTO annotations
       (candidate_id, score_pre_open, score_sync_open, score_sync_repress,
        score_unrelated, overall_confidence, review_status, notes, annotator,
        created_at, updated_at, app_version)
     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
     ON CONFLICT(candidate_id) DO UPDATE SET
       score_pre_open = excluded.score_pre_open,
       score_sync_open = excluded.score_sync_open,
       score_sync_repress = excluded.score_sync_repress,
       score_unrelated = excluded.score_unrelated,
       overall_confidence = excluded.overall_confidence,
       review_status = excluded.review_status,
       notes = excluded.notes,
       annotator = excluded.annotator,
       updated_at = excluded.updated_at,
       app_version = excluded.app_version",
    params = list(
      candidate_id,
      as.numeric(scores[["pre_open"]]),
      as.numeric(scores[["sync_open"]]),
      as.numeric(scores[["sync_repress"]]),
      as.numeric(scores[["unrelated"]]),
      as.numeric(overall_confidence),
      review_status %||% NA_character_,
      notes %||% NA_character_,
      annotator %||% NA_character_,
      created_at,
      now,
      app_version %||% NA_character_
    )
  )
  DBI::dbCommit(con)
  invisible(TRUE)
}

# Mark as skipped with review_status='skipped' and NULL scores
mark_skipped <- function(con, candidate_id, annotator = NULL, app_version = NULL) {
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
  old <- annotation_row(con, candidate_id)
  created_at <- if (!is.null(old) && !is.na(old$created_at)) old$created_at else now
  DBI::dbBegin(con)
  DBI::dbExecute(
    con,
    "INSERT INTO annotations
       (candidate_id, score_pre_open, score_sync_open, score_sync_repress,
        score_unrelated, overall_confidence, review_status, notes, annotator,
        created_at, updated_at, app_version)
     VALUES (?, NULL, NULL, NULL, NULL, NULL, 'skipped', NULL, ?, ?, ?, ?)
     ON CONFLICT(candidate_id) DO UPDATE SET
       review_status = 'skipped',
       updated_at = excluded.updated_at,
       annotator = excluded.annotator",
    params = list(
      candidate_id,
      annotator %||% NA_character_,
      created_at,
      now,
      app_version %||% NA_character_
    )
  )
  DBI::dbCommit(con)
  invisible(TRUE)
}

# Mark for review, updating only the status when scores already exist
mark_review <- function(con, candidate_id, annotator = NULL) {
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
  DBI::dbBegin(con)
  DBI::dbExecute(
    con,
    "INSERT INTO annotations
       (candidate_id, review_status, annotator, created_at, updated_at)
     VALUES (?, 'review', ?, ?, ?)
     ON CONFLICT(candidate_id) DO UPDATE SET
       review_status = 'review', updated_at = excluded.updated_at,
       annotator = excluded.annotator",
    params = list(candidate_id, annotator %||% NA_character_, now, now)
  )
  DBI::dbCommit(con)
  invisible(TRUE)
}

# List candidates with annotation status. queue.R handles sorting by priority and stable hash.
fetch_candidates <- function(con, filters = list(), include_unqueued = FALSE) {
  where <- character()
  params <- list()
  add_where <- function(cond, value) {
    where <<- c(where, cond)
    if (!is.null(value)) {
      params <<- c(params, list(value))
    }
  }

  if (!is.null(filters$context) && length(filters$context) && !identical(filters$context, "all")) {
    add_where("c.context = ?", filters$context)
  }
  if (!include_unqueued) {
    add_where("c.queue_priority IS NOT NULL", NULL)
  }
  if (!is.null(filters$source) && length(filters$source) && !identical(filters$source, "all")) {
    src <- filters$source
    if (src == "mandatory") {
      add_where("c.mandatory = 1", NULL)
    } else if (src == "sync_repress_heuristic") {
      add_where("c.source_sync_repress_heuristic = 1", NULL)
    } else if (src == "random_background") {
      add_where("c.source_sync_repress_heuristic = 0 AND c.mandatory = 0", NULL)
    } else {
      stop("未知候选来源: ", src)
    }
  }
  if (!is.null(filters$only_unannotated) && isTRUE(filters$only_unannotated)) {
    add_where("a.candidate_id IS NULL", NULL)
  }
  if (!is.null(filters$only_skipped) && isTRUE(filters$only_skipped)) {
    add_where("a.review_status = 'skipped'", NULL)
  }
  if (!is.null(filters$only_review) && isTRUE(filters$only_review)) {
    add_where("a.review_status = 'review'", NULL)
  }
  if (!is.null(filters$only_batch) && isTRUE(filters$only_batch)) {
    add_where("c.batch_id IS NOT NULL", NULL)
  }
  if (!is.null(filters$gene_pattern) && nzchar(filters$gene_pattern)) {
    add_where("c.gene LIKE ?", paste0("%", filters$gene_pattern, "%"))
  }
  if (!is.null(filters$peak_pattern) && nzchar(filters$peak_pattern)) {
    add_where("c.peak LIKE ?", paste0("%", filters$peak_pattern, "%"))
  }
  if (!is.null(filters$dist_min) && !is.na(filters$dist_min)) {
    add_where("c.distance >= ?", as.numeric(filters$dist_min))
  }
  if (!is.null(filters$dist_max) && !is.na(filters$dist_max)) {
    add_where("c.distance <= ?", as.numeric(filters$dist_max))
  }

  sql <- paste0(
    "SELECT c.candidate_id, c.context, c.peak, c.gene, c.distance,
            c.mandatory, c.source_tri_pre, c.source_sync_open,
            c.source_sync_repress_heuristic, c.S_pre_shape, c.S_final,
            c.queue_priority, c.is_repeat, c.repeat_of,
            c.batch_id, c.batch_score,
            a.score_pre_open, a.score_sync_open, a.score_sync_repress,
            a.score_unrelated, a.overall_confidence, a.review_status,
            a.notes, a.updated_at
     FROM candidates c
     LEFT JOIN annotations a ON a.candidate_id = c.candidate_id",
    if (length(where)) paste0(" WHERE ", paste(where, collapse = " AND ")) else ""
  )
  if (length(params)) {
    DBI::dbGetQuery(con, sql, params = params)
  } else {
    DBI::dbGetQuery(con, sql)
  }
}

count_candidates <- function(con, filters = list(), include_unqueued = FALSE) {
  res <- fetch_candidates(con, filters = filters, include_unqueued = include_unqueued)
  nrow(res)
}

count_annotated <- function(con, filters = list()) {
  filters$only_unannotated <- NULL
  res <- fetch_candidates(con, filters = filters)
  sum(!is.na(res$score_pre_open) & !is.na(res$score_sync_open) &
        !is.na(res$score_sync_repress) & !is.na(res$score_unrelated))
}

fetch_all_annotations <- function(con) {
  DBI::dbGetQuery(
    con,
    "SELECT c.candidate_id, c.context, c.peak, c.gene, c.distance,
            c.mandatory, c.source_tri_pre, c.source_sync_open,
            c.source_sync_repress_heuristic, c.S_pre_shape, c.S_final,
            c.queue_priority, c.is_repeat, c.repeat_of,
            c.batch_id, c.batch_score,
            a.score_pre_open, a.score_sync_open, a.score_sync_repress,
            a.score_unrelated, a.overall_confidence, a.review_status,
            a.notes, a.annotator, a.created_at, a.updated_at, a.app_version
     FROM candidates c
     LEFT JOIN annotations a ON a.candidate_id = c.candidate_id
     WHERE a.score_pre_open IS NOT NULL
        OR a.score_sync_open IS NOT NULL
        OR a.score_sync_repress IS NOT NULL
        OR a.score_unrelated IS NOT NULL
        OR a.review_status IN ('skipped', 'review')
     ORDER BY a.updated_at")
}
