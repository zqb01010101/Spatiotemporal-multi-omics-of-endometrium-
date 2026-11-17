# =============================================================================
# label_rules.R — Derive hard_label from the four annotation scores
#
# Thresholds come from config.yml:
#   Top score >= hard_label_min_top
#   and top score minus second-highest score >= hard_label_min_margin
#   and overall_confidence >= hard_label_min_confidence
# Assign the highest-scoring class when all thresholds pass. Otherwise use "uncertain".
# Also return top_score, second_score, and score_margin.
# =============================================================================

SCORE_NAMES <- c(
  pre_open     = "score_pre_open",
  sync_open    = "score_sync_open",
  sync_repress = "score_sync_repress",
  unrelated    = "score_unrelated"
)

derive_hard_label <- function(scores, overall_confidence, cfg) {
  s <- as.numeric(scores)
  if (length(s) != 4L || anyNA(s)) {
    return(list(
      hard_label = NA_character_,
      top_score = NA_real_,
      second_score = NA_real_,
      score_margin = NA_real_
    ))
  }
  ord <- order(s, decreasing = TRUE)
  top <- s[ord[1]]
  second <- s[ord[2]]
  margin <- top - second

  label <- "uncertain"
  if (!is.na(top) && !is.na(second) && !is.na(overall_confidence) &&
      top >= cfg$label_rules$hard_label_min_top &&
      margin >= cfg$label_rules$hard_label_min_margin &&
      overall_confidence >= cfg$label_rules$hard_label_min_confidence) {
    label <- names(SCORE_NAMES)[ord[1]]
  }

  list(
    hard_label   = label,
    top_score    = top,
    second_score = second,
    score_margin = margin
  )
}

# Validate numeric scores from 0 to 10. NA is allowed, and all NA means unlabeled.
validate_scores <- function(scores) {
  s <- as.numeric(scores)
  if (any(!is.na(s) & (s < 0 | s > 10))) {
    stop("评分必须在 0-10 之间")
  }
  invisible(TRUE)
}
