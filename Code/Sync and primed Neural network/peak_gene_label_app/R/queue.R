# =============================================================================
# queue.R — Labeling queue ordering and recovery
#
# Ordering rules from config.yml queue_priorities:
#   mandatory(1) > sync_repress_heuristic(2) > random_background(3)
# Within each priority, sort pseudorandomly by the MD5 hash of candidate_id.
# This keeps order stable across sessions and reduces runs of similar candidates.
# =============================================================================

stable_hash_order <- function(ids) {
  h <- vapply(ids, function(id) {
    digest::digest(id, algo = "md5", serialize = FALSE)
  }, character(1))
  order(h)
}

# Return rows of the candidates table in queue order using queue_priority and candidate_id
build_queue <- function(cand_dt) {
  dt <- data.table::as.data.table(cand_dt)
  if (!nrow(dt)) {
    return(dt)
  }
  dt <- dt[!is.na(queue_priority)]
  dt[, h := vapply(candidate_id, function(id) {
    digest::digest(id, algo = "md5", serialize = FALSE)
  }, character(1))]
  data.table::setorder(dt, queue_priority, h)
  dt[, h := NULL]
  dt[]
}

# Resume at the first candidate without a label, skip, or review status
resume_position <- function(queue_dt) {
  if (!nrow(queue_dt)) {
    return(0L)
  }
  done <- !is.na(queue_dt$score_pre_open) |
          !is.na(queue_dt$score_sync_open) |
          !is.na(queue_dt$score_sync_repress) |
          !is.na(queue_dt$score_unrelated) |
          (!is.na(queue_dt$review_status) & queue_dt$review_status %in% c("skipped", "review"))
  first_open <- which(!done)
  if (length(first_open)) first_open[1] else 0L
}

# Count all candidates and labeled candidates under the current filters
queue_progress <- function(queue_dt) {
  n_total <- nrow(queue_dt)
  n_done <- sum(!is.na(queue_dt$score_pre_open) |
                !is.na(queue_dt$score_sync_open) |
                !is.na(queue_dt$score_sync_repress) |
                !is.na(queue_dt$score_unrelated))
  n_review <- sum(!is.na(queue_dt$review_status) & queue_dt$review_status == "review")
  n_skipped <- sum(!is.na(queue_dt$review_status) & queue_dt$review_status == "skipped")
  list(
    total = n_total,
    done = n_done,
    review = n_review,
    skipped = n_skipped,
    remaining = n_total - n_done - n_skipped
  )
}
