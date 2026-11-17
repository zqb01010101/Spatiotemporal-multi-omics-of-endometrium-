# =============================================================================
# build_candidates.R — Build the candidate pair index
#
# Candidates are sparse peak–gene entries with res_peak_gene_1Mb$mat == 1.
#       Treat stromal and epithelial contexts of the same pair as separate samples.
# Queue priority: mandatory (1), sync_repress_heuristic (2), then random_background (3).
# Read all thresholds from config.yml.
# =============================================================================

# -----------------------------------------------------------------------------
# 1) Select entries equal to 1 in res$mat and attach TSS distances via sparse element-wise multiplication
# -----------------------------------------------------------------------------
extract_res_candidates <- function(res_mat, dist_mat, log_file = NULL) {
  s <- Matrix::summary(res_mat)
  log_msg("res$mat 非零元素：", nrow(s), file = log_file)
  bad <- s$x[abs(s$x - 1) > 1e-12]
  if (length(bad)) {
    log_msg("res$mat 中存在 ", length(bad), " 个非 1 的非零值，已排除",
            file = log_file, level = "WARN")
    s <- s[abs(s$x - 1) < 1e-12, , drop = FALSE]
  }
  rn <- rownames(res_mat)
  cn <- colnames(res_mat)
  dt <- data.table::data.table(
    peak = rn[s$i],
    gene = cn[s$j],
    distance = NA_real_,
    key_ = make_pair_key(rn[s$i], cn[s$j])
  )
  data.table::setkey(dt, key_)

  # Multiply dist_mat by the sparse binary mask P to keep distances only for candidate pairs
  P <- Matrix::sparseMatrix(i = s$i, j = s$j, x = 1,
                            dims = dim(res_mat), dimnames = dimnames(res_mat))
  dc <- dist_mat * P
  sd <- Matrix::summary(dc)
  pos_key <- make_pair_key(rn[sd$i], cn[sd$j])
  dist_vec <- rep(0, nrow(dt))
  m <- match(dt$key_, pos_key)
  hit <- !is.na(m)
  dist_vec[hit] <- sd$x[m[hit]]
  dt[, distance := dist_vec]

  n_zero_dist <- sum(dt$distance == 0)
  log_msg("候选 peak–gene 对（res 中值为 1）：", nrow(dt),
          "，其中距离为 0 的 ", n_zero_dist, " 对", file = log_file)
  dt[]
}

# -----------------------------------------------------------------------------
# 2) Build the mandatory queue from tri (S_pre_shape > threshold) and sync (S_final > threshold), validating IDs
# -----------------------------------------------------------------------------
add_mandatory_sources <- function(cand_dt, inputs, cfg, log_file = NULL) {
  common_peaks <- rownames(inputs$Y_sig)
  common_genes <- rownames(inputs$avg)
  row_set <- cand_dt[, unique(row_key_)]
  excl_rows <- list()

  specs <- list(
    list(file = "tri_stroma",  context = "stroma",
         source = "tri_pre", s_col = cfg$columns$S_pre_shape$tri_stroma,
         threshold = cfg$mandatory$tri_S_pre_shape_threshold),
    list(file = "tri_epi",     context = "epithelial",
         source = "tri_pre", s_col = cfg$columns$S_pre_shape$tri_epi,
         threshold = cfg$mandatory$tri_S_pre_shape_threshold),
    list(file = "sync_stroma", context = "stroma",
         source = "sync_open", s_col = cfg$columns$S_final$sync_stroma,
         threshold = cfg$mandatory$sync_S_final_threshold),
    list(file = "sync_epi",    context = "epithelial",
         source = "sync_open", s_col = cfg$columns$S_final$sync_epi,
         threshold = cfg$mandatory$sync_S_final_threshold)
  )

  for (sp in specs) {
    df <- inputs[[sp$file]]
    pk_col <- cfg$columns$peak[[sp$file]]
    gn_col <- cfg$columns$gene[[sp$file]]
    d <- data.table::data.table(
      context = sp$context,
      peak = as.character(df[[pk_col]]),
      gene = as.character(df[[gn_col]]),
      S    = as.numeric(df[[sp$s_col]])
    )
    d <- d[!is.na(S) & S > sp$threshold]
    log_msg(sprintf("[%s] %s>%.3f 的行数：%d / %d",
                    sp$file, sp$s_col, sp$threshold, nrow(d), nrow(df)),
            file = log_file)

    d[, row_key_ := paste(context, peak, gene, sep = "|")]
    d[, pk_ok := peak %chin% common_peaks]
    d[, gn_ok := gene %chin% common_genes]
    d[, pair_ok := row_key_ %chin% row_set]
    d[, reason := fifelse(!pk_ok, "peak 不在 Y_sig",
                   fifelse(!gn_ok, "gene 不在 avg",
                   fifelse(!pair_ok, "res$mat 中对应值非 1/不存在", NA_character_)))]

    ok <- d[pk_ok & gn_ok & pair_ok]
    bad <- d[!pk_ok | !gn_ok | !pair_ok]
    log_msg(sprintf("[%s] mandatory 通过：%d；排除：%d",
                    sp$file, nrow(ok), nrow(bad)), file = log_file)
    if (nrow(bad)) {
      tab <- bad[, .N, by = reason][order(-N)]
      excl_rows[[length(excl_rows) + 1L]] <- data.table::data.table(
        file = sp$file,
        context = sp$context,
        source = sp$source,
        reason = tab$reason,
        n_excluded = tab$N
      )
      for (i in seq_len(nrow(tab))) {
        log_msg(sprintf("[%s] 排除原因「%s」：%d 行，示例：%s",
                        sp$file, tab$reason[i], tab$N[i],
                        paste(utils::head(bad[reason == tab$reason[i],
                                              paste(peak, gene, sep = "|")], 3),
                              collapse = " ; ")),
                file = log_file, level = "WARN")
      }
    }

    if (nrow(ok)) {
      # For duplicate pairs from the same source, keep the highest score
      ok <- ok[, .(row_key_ = row_key_, S = max(S)), by = row_key_]
      if (sp$source == "tri_pre") {
        cand_dt[ok, on = "row_key_", `:=`(
          mandatory = 1L,
          source_tri_pre = 1L,
          queue_priority = as.integer(cfg$queue_priorities$mandatory),
          S_pre_shape = data.table::fifelse(
            is.na(S_pre_shape) | is.na(i.S), data.table::fcoalesce(S_pre_shape, i.S),
            pmax(S_pre_shape, i.S))
        )]
      } else {
        cand_dt[ok, on = "row_key_", `:=`(
          mandatory = 1L,
          source_sync_open = 1L,
          queue_priority = as.integer(cfg$queue_priorities$mandatory),
          S_final = data.table::fifelse(
            is.na(S_final) | is.na(i.S), data.table::fcoalesce(S_final, i.S),
            pmax(S_final, i.S))
        )]
      }
    }
    if (!nrow(bad)) {
      excl_rows[[length(excl_rows) + 1L]] <- data.table::data.table(
        file = sp$file, context = sp$context, source = sp$source,
        reason = NA_character_, n_excluded = 0L
      )
    }
  }
  list(cand = cand_dt, exclusions = data.table::rbindlist(excl_rows, fill = TRUE))
}

# -----------------------------------------------------------------------------
# 3) Synchronous repression heuristic for sampling only
#    Both profiles must have sufficient dynamic range and at least one pronounced decrease.
#    The steepest decreases must occur in the same or adjacent intervals.
# -----------------------------------------------------------------------------
add_sync_repress_heuristic <- function(cand_dt, inputs, cfg, log_file = NULL) {
  sc <- cfg$sampling$sync_repress
  seed <- cfg$sampling$random_seed
  ctx_idx <- 0L
  for (ctx in c("stroma", "epithelial")) {
    ctx_idx <- ctx_idx + 1L
    sec_cols <- sec_cols_for(cfg, ctx)
    pool <- cand_dt[context == ctx & mandatory == 0L]
    set.seed(seed + ctx_idx * 1000L)
    n_sample <- min(sc$sample_size, nrow(pool))
    smp <- pool[sample(.N, n_sample)]
    if (!nrow(smp)) {
      log_msg(sprintf("[heuristic/%s] 无可抽样候选", ctx), file = log_file)
      next
    }

    vals <- extract_bulk_values(inputs, cfg, ctx, smp$peak, smp$gene)
    pcols <- paste0("peak_", cfg$stages$labels)
    gcols <- paste0("gene_", cfg$stages$labels)
    pmat <- as.matrix(vals[, pcols, with = FALSE])
    gmat <- as.matrix(vals[, gcols, with = FALSE])

    d <- data.table::data.table(row_key_ = smp$row_key_)
    d[, range_peak := apply(pmat, 1L, function(x) diff(range(x)))]
    d[, range_gene := apply(gmat, 1L, function(x) diff(range(x)))]
    dp <- t(apply(pmat, 1L, diff))
    dg <- t(apply(gmat, 1L, diff))
    d[, min_drop_peak := apply(dp, 1L, min)]
    d[, min_drop_gene := apply(dg, 1L, min)]
    d[, i_peak := apply(dp, 1L, which.min)]
    d[, i_gene := apply(dg, 1L, which.min)]
    d[, hit := range_peak >= sc$dynamic_range_peak &
                range_gene >= sc$dynamic_range_gene &
                min_drop_peak <= -sc$min_drop_peak &
                min_drop_gene <= -sc$min_drop_gene &
                abs(i_peak - i_gene) <= sc$max_window_gap]
    hits <- d[(hit)]
    log_msg(sprintf("[heuristic/%s] 抽样 %d，命中 %d", ctx, nrow(smp), nrow(hits)),
            file = log_file)

    if (nrow(hits)) {
      set.seed(seed + ctx_idx * 1000L + 1L)
      take_n <- min(sc$max_per_context, nrow(hits))
      take <- hits[sample(.N, take_n)]
      cand_dt[take, on = "row_key_", `:=`(
        source_sync_repress_heuristic = 1L,
        queue_priority = as.integer(cfg$queue_priorities$sync_repress_heuristic)
      )]
      log_msg(sprintf("[heuristic/%s] 纳入队列：%d", ctx, take_n), file = log_file)
    }
  }
  cand_dt[]
}

# -----------------------------------------------------------------------------
# 4) Draw a stratified random background sample from remaining pairs in each context
# -----------------------------------------------------------------------------
add_background <- function(cand_dt, cfg, log_file = NULL) {
  seed <- cfg$sampling$random_seed
  n_bg <- cfg$sampling$n_background_per_context
  ctx_idx <- 0L
  for (ctx in c("stroma", "epithelial")) {
    ctx_idx <- ctx_idx + 1L
    pool <- cand_dt[context == ctx & mandatory == 0L &
                      source_sync_repress_heuristic == 0L]
    set.seed(seed + ctx_idx * 2000L)
    take_n <- min(n_bg, nrow(pool))
    take <- pool[sample(.N, take_n)]
    cand_dt[take, on = "row_key_",
            queue_priority := as.integer(cfg$queue_priorities$random_background)]
    log_msg(sprintf("[background/%s] 纳入队列：%d（池 %d）", ctx, take_n, nrow(pool)),
            file = log_file)
  }
  cand_dt[]
}

# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
build_candidate_index <- function(inputs, cfg, log_file = NULL) {
  base <- extract_res_candidates(inputs$res_mat, inputs$dist_mat, log_file)
  base[, `:=`(
    mandatory = 0L,
    source_tri_pre = 0L,
    source_sync_open = 0L,
    source_sync_repress_heuristic = 0L,
    S_pre_shape = NA_real_,
    S_final = NA_real_,
    queue_priority = NA_integer_
  )]

  build_ctx <- function(ctx) {
    d <- data.table::copy(base)
    d[, context := ctx]
    d[]
  }
  all <- data.table::rbindlist(list(
    build_ctx("stroma"),
    build_ctx("epithelial")
  ))
  data.table::alloc.col(all, 40L)
  all[, row_key_ := paste(context, peak, gene, sep = "|")]
  data.table::setkey(all, row_key_)

  mand <- add_mandatory_sources(all, inputs, cfg, log_file)
  all <- mand$cand
  exclusions <- mand$exclusions
  all <- add_sync_repress_heuristic(all, inputs, cfg, log_file)
  all <- add_background(all, cfg, log_file)

  all[, `:=`(
    candidate_id = row_key_,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
    is_repeat = 0L,
    repeat_of = NA_character_
  )]
  data.table::setcolorder(all, c(
    "candidate_id", "context", "peak", "gene", "distance",
    "mandatory", "source_tri_pre", "source_sync_open",
    "source_sync_repress_heuristic", "S_pre_shape", "S_final",
    "queue_priority", "created_at", "is_repeat", "repeat_of"
  ))
  data.table::set(all, j = c("key_", "row_key_"), value = NULL)
  list(cand = all, exclusions = exclusions)
}

# Summary statistics for scripts and tests
summarize_candidates <- function(cand_dt) {
  n_total <- nrow(cand_dt)
  list(
    total = n_total,
    total_stroma = nrow(cand_dt[context == "stroma"]),
    total_epithelial = nrow(cand_dt[context == "epithelial"]),
    mandatory_stroma = nrow(cand_dt[context == "stroma" & mandatory == 1L]),
    mandatory_epithelial = nrow(cand_dt[context == "epithelial" & mandatory == 1L]),
    mandatory_total = nrow(cand_dt[mandatory == 1L]),
    heuristic_stroma = nrow(cand_dt[context == "stroma" & source_sync_repress_heuristic == 1L]),
    heuristic_epithelial = nrow(cand_dt[context == "epithelial" & source_sync_repress_heuristic == 1L]),
    background_stroma = nrow(cand_dt[context == "stroma" & queue_priority == 3L]),
    background_epithelial = nrow(cand_dt[context == "epithelial" & queue_priority == 3L]),
    queue_total = nrow(cand_dt[!is.na(queue_priority)]),
    source_tri_pre = nrow(cand_dt[source_tri_pre == 1L]),
    source_sync_open = nrow(cand_dt[source_sync_open == 1L])
  )
}
