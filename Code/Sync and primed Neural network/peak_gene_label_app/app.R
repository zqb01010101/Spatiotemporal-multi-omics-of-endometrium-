# =============================================================================
# app.R — Shiny app for annotating peak–gene temporal relationships
#
# Launch: R -e "shiny::runApp('.', launch.browser = TRUE)"
#
# Features:
#   * Four 0–10 scores, confidence, and notes. Saves persist immediately in SQLite.
#   * Peak and gene profiles, raw-value overlays, changes between adjacent stages, and a value table.
#   * Blinded mode hides S_pre_shape and S_final before annotation, with an option to reveal them.
#   * Filters for label status, mandatory candidates, pending review, skipped candidates, context, source, search, and distance.
#   * Jump directly to candidate N in the filtered queue.
#   * Shortcuts: Ctrl+Enter to save and advance, arrow keys to navigate, S to skip, and R to mark for review.
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(DT)
})

detect_app_dir <- function() {
  env_dir <- Sys.getenv("PEAK_GENE_APP_DIR", unset = "")
  if (nzchar(env_dir)) return(normalizePath(env_dir, mustWork = TRUE))
  args0 <- commandArgs(FALSE)
  f <- sub("^--file=", "", args0[grepl("^--file=", args0)])
  if (length(f) && nzchar(f[1])) {
    return(normalizePath(dirname(dirname(f[1])), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

app_dir <- detect_app_dir()
for (f in c("load_inputs.R", "validate_inputs.R", "build_candidates.R",
            "extract_pair_values.R", "make_plots.R", "database.R",
            "queue.R", "label_rules.R")) {
  source(file.path(app_dir, "R", f))
}

cfg <- read_config(app_dir)
app_version <- cfg$app$app_version %||% "1.0.0"

# -----------------------------------------------------------------------------
# Load read-only inputs once at app startup
# -----------------------------------------------------------------------------
fixture_dir <- Sys.getenv("PEAK_GENE_FIXTURE_DIR", unset = "")
use_snapshot <- !nzchar(fixture_dir)
if (nzchar(fixture_dir)) {
  # Test mode uses small fixtures
  cfg <- yaml::read_yaml(file.path(fixture_dir, "config.yml"))
  cfg$app_dir <- normalizePath(fixture_dir)
  cat("测试模式：使用夹具数据 ", cfg$app_dir, "\n")
  inputs <- load_inputs(cfg)
  validate_inputs(inputs, cfg)
} else {
  # Production mode reads the unified pair table instead of the large matrices
  unified_file <- file.path(app_dir, "data_cache", "pair_features_all.csv")
  if (!file.exists(unified_file)) {
    stop("缺少统一配对数据文件：", unified_file,
         "，请先运行 scripts/21_build_unified_pairs.R")
  }
  pair_features <- data.table::fread(unified_file)
  data.table::setkey(pair_features, candidate_id)
  cat("已加载统一配对数据：", unified_file, "（", nrow(pair_features), " 条）\n")
}
cat("输入加载完成。\n")

queue_source_label <- function(r) {
  if (!is.null(r$batch_id) && length(r$batch_id) && !is.na(r$batch_id)) {
    return(paste0("复核批次:", r$batch_id))
  }
  if (isTRUE(r$is_repeat == 1)) return("一致性重复样本")
  if (isTRUE(r$mandatory == 1)) return("mandatory")
  if (isTRUE(r$source_sync_repress_heuristic == 1)) return("sync_repress_heuristic")
  if (!is.na(r$queue_priority) && r$queue_priority == 3) return("random_background")
  "其他"
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- page_sidebar(
  title = div(
    style = "display:flex; align-items:center; gap:12px;",
    "Peak–Gene 时序关系人工标注平台",
    span(style = "font-size:13px; color:#888;", "盲标默认开启 · 保存即落盘")
  ),
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 300,
    h5("筛选"),
    checkboxInput("f_unannotated", "只看未标注", value = FALSE),
    checkboxInput("f_mandatory", "只看 mandatory", value = FALSE),
    checkboxInput("f_review", "只看待复核", value = FALSE),
    checkboxInput("f_skipped", "只看已跳过", value = FALSE),
    checkboxInput("f_batch", "只看新复核批次", value = FALSE),
    selectInput("f_context", "context",
                choices = c("全部" = "all", "stroma" = "stroma", "epithelial" = "epithelial")),
    selectInput("f_source", "候选来源",
                choices = c("全部" = "all",
                            "mandatory" = "mandatory",
                            "同步抑制启发式" = "sync_repress_heuristic",
                            "随机背景" = "random_background")),
    textInput("f_gene", "gene 搜索", placeholder = "如 HES4"),
    textInput("f_peak", "peak 搜索", placeholder = "如 chr1:802245-802745"),
    numericInput("f_dist_min", "距离下限 (bp)", value = NA, step = 1000),
    numericInput("f_dist_max", "距离上限 (bp)", value = NA, step = 1000),
    actionButton("f_reset", "重置筛选", width = "100%"),
    hr(),
    textInput("annotator", "标注者", value = cfg$app$annotator_default %||% "annotator1"),
    p(style = "font-size:12px; color:#999;",
      "快捷键：Ctrl+Enter 保存并下一条 · ← 上一条 · → 下一条 · S 跳过 · R 复核")
  ),

  tags$head(
    tags$script(HTML("
      $(document).on('keydown', function(e) {
        var tag = (e.target.tagName || '').toLowerCase();
        var typing = tag === 'input' || tag === 'textarea' || tag === 'select' ||
                     (e.target.isContentEditable === true);
        if (e.ctrlKey && e.key === 'Enter') {
          e.preventDefault();
          Shiny.setInputValue('key_event', 'save_next', {priority: 'event'});
          return;
        }
        if (typing) return;
        if (e.key === 'ArrowLeft') {
          e.preventDefault();
          Shiny.setInputValue('key_event', 'prev', {priority: 'event'});
        } else if (e.key === 'ArrowRight') {
          e.preventDefault();
          Shiny.setInputValue('key_event', 'next', {priority: 'event'});
        } else if (e.key === 's' || e.key === 'S') {
          e.preventDefault();
          Shiny.setInputValue('key_event', 'skip', {priority: 'event'});
        } else if (e.key === 'r' || e.key === 'R') {
          e.preventDefault();
          Shiny.setInputValue('key_event', 'review', {priority: 'event'});
        }
      });
    "))
  ),

  card(
    card_header("当前关系", class = "py-1"),
    uiOutput("status_boxes"),
    uiOutput("evidence_panel")
  ),

  navset_card_underline(
    nav_panel("图形",
        layout_columns(
          card(plotlyOutput("plot_peak", height = 400)),
          card(plotlyOutput("plot_gene", height = 400)),
          card(plotlyOutput("plot_overlay", height = 400)),
          card(plotlyOutput("plot_deltas", height = 460))
        )),
    nav_panel("数值表",
        card(
          h5("原始数值"),
          DTOutput("raw_table"),
          h5("相邻时期变化"),
          DTOutput("delta_table")
        ))
  ),

  card(
    card_header("人工标签（四个关系评分不强制总和为 10）"),
    layout_columns(
      sliderInput("score_pre_open", "预开放", 0, 10, 0, step = 1),
      sliderInput("score_sync_open", "同步开放", 0, 10, 0, step = 1),
      sliderInput("score_sync_repress", "同步抑制", 0, 10, 0, step = 1),
      sliderInput("score_unrelated", "无关系 / 无法判断", 0, 10, 0, step = 1),
      sliderInput("confidence", "总体判断信心", 0, 10, 0, step = 1)
    ),
    textAreaInput("notes", "备注", width = "100%", rows = 3),
    actionButton("save_next", "保存并下一条", width = "100%", class = "btn-primary"),
    layout_columns(col_widths = c(4, 4, 4),
      actionButton("save_only", "仅保存"),
      actionButton("prev", "上一条"),
      actionButton("next_btn", "下一条")),
    layout_columns(col_widths = c(4, 4, 4),
      actionButton("skip", "暂时跳过"),
      actionButton("review", "标记稍后复核"),
      actionButton("clear_scores", "清空当前评分")),
    layout_columns(col_widths = c(8, 4),
      numericInput("jump_to", "跳到当前队列第 N 个关系",
                   value = 1, min = 1, step = 1),
      div(style = "display:flex; align-items:flex-end;",
          actionButton("jump_btn", "跳转", width = "100%"))),
    actionButton("toggle_evidence", "显示统计证据", width = "100%")
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {
  con <- init_db(cfg)
  onSessionEnded(function() {
    try(DBI::dbDisconnect(con), silent = TRUE)
  })

  queue_state <- reactiveVal(list(ids = character(), pos = 0L, total = 0L, dt = data.frame()))
  current_id <- reactiveVal(NULL)
  show_evidence <- reactiveVal(isFALSE(cfg$app$blind_default %||% TRUE))

  filters <- reactive({
    num_input <- function(x) {
      v <- suppressWarnings(as.numeric(x))
      if (length(v) == 0L || is.na(v)) NA_real_ else v
    }
    list(
      context = input$f_context %||% "all",
      source = input$f_source %||% "all",
      only_unannotated = isTRUE(input$f_unannotated),
      only_review = isTRUE(input$f_review),
      only_skipped = isTRUE(input$f_skipped),
      only_batch = isTRUE(input$f_batch),
      gene_pattern = trimws(input$f_gene %||% ""),
      peak_pattern = trimws(input$f_peak %||% ""),
      dist_min = num_input(input$f_dist_min),
      dist_max = num_input(input$f_dist_max)
    )
  })

  refresh_queue <- reactive({
    flt <- filters()
    q <- fetch_candidates(con, flt, include_unqueued = FALSE)
    qd <- as.data.table(q)
    if (nrow(qd)) build_queue(qd) else qd
  })

  observe({
    qd <- refresh_queue()
    if (!nrow(qd)) {
      queue_state(list(ids = character(), pos = 0L, total = 0L, dt = qd))
      current_id(NULL)
      return()
    }
    cur <- current_id()
    if (!is.null(cur) && cur %in% qd$candidate_id) {
      pos <- which(qd$candidate_id == cur)[1]
    } else {
      pos <- resume_position(qd)
      if (pos == 0L) pos <- 1L
      current_id(qd$candidate_id[pos])
    }
    queue_state(list(ids = qd$candidate_id, pos = pos, total = nrow(qd), dt = qd))
  })

  go_to <- function(delta) {
    st <- queue_state()
    if (!st$total) return()
    n <- length(st$ids)
    newpos <- ((st$pos - 1L + delta) %% n) + 1L
    current_id(st$ids[newpos])
  }

  observeEvent(input$prev, go_to(-1L))
  observeEvent(input$next_btn, go_to(1L))

  # Jump to candidate N in the filtered queue
  observeEvent(input$jump_btn, {
    st <- queue_state()
    if (!st$total) {
      showNotification("当前队列为空，无法跳转", type = "warning")
      return()
    }
    n <- suppressWarnings(as.integer(input$jump_to))
    if (length(n) == 0L || is.na(n) || n < 1L) n <- 1L
    n <- min(n, st$total)
    current_id(st$ids[n])
    showNotification(sprintf("已跳转到当前队列第 %d / %d 个关系", n, st$total),
                     type = "message")
  })

  # Update the jump field when the queue position changes
  observe({
    st <- queue_state()
    if (st$total > 0L) {
      updateNumericInput(session, "jump_to", value = st$pos)
    }
  })

  observeEvent(input$f_reset, {
    updateCheckboxInput(session, "f_unannotated", value = FALSE)
    updateCheckboxInput(session, "f_mandatory", value = FALSE)
    updateCheckboxInput(session, "f_review", value = FALSE)
    updateCheckboxInput(session, "f_skipped", value = FALSE)
    updateCheckboxInput(session, "f_batch", value = FALSE)
    updateSelectInput(session, "f_context", selected = "all")
    updateSelectInput(session, "f_source", selected = "all")
    updateTextInput(session, "f_gene", value = "")
    updateTextInput(session, "f_peak", value = "")
    updateNumericInput(session, "f_dist_min", value = NA)
    updateNumericInput(session, "f_dist_max", value = NA)
  })

  current_row <- reactive({
    id <- current_id()
    st <- queue_state()
    if (is.null(id) || !nrow(st$dt)) return(NULL)
    st$dt[candidate_id == id][1]
  })

  current_values <- reactive({
    r <- current_row()
    if (is.null(r)) return(NULL)
    if (use_snapshot) {
      row <- pair_features[candidate_id == r$candidate_id]
      if (!nrow(row)) return(NULL)
      stages <- cfg$stages$labels
      return(list(
        peak = as.numeric(row[, paste0("peak_", stages), with = FALSE]),
        gene = as.numeric(row[, paste0("gene_", stages), with = FALSE]),
        stages = stages
      ))
    }
    extract_pair_values(inputs, cfg, r$context, r$peak, r$gene)
  })

  # Restore saved annotations on navigation. Unlabeled candidates start blinded.
  observe({
    r <- current_row()
    if (is.null(r)) return()
    ann <- annotation_row(con, r$candidate_id)
    if (!is.null(ann) && !is.na(ann$score_pre_open)) {
      updateSliderInput(session, "score_pre_open", value = ann$score_pre_open)
      updateSliderInput(session, "score_sync_open", value = ann$score_sync_open)
      updateSliderInput(session, "score_sync_repress", value = ann$score_sync_repress)
      updateSliderInput(session, "score_unrelated", value = ann$score_unrelated)
      updateSliderInput(session, "confidence", value = ann$overall_confidence)
      updateTextAreaInput(session, "notes", value = ann$notes %||% "")
      show_evidence(TRUE)   # Reveal statistical evidence after saving
    } else {
      updateSliderInput(session, "score_pre_open", value = 0)
      updateSliderInput(session, "score_sync_open", value = 0)
      updateSliderInput(session, "score_sync_repress", value = 0)
      updateSliderInput(session, "score_unrelated", value = 0)
      updateSliderInput(session, "confidence", value = 0)
      updateTextAreaInput(session, "notes", value = "")
      show_evidence(isFALSE(cfg$app$blind_default %||% TRUE))
    }
  })

  do_save <- function(advance = FALSE) {
    r <- current_row()
    if (is.null(r)) {
      showNotification("没有可保存的关系", type = "warning")
      return(invisible(FALSE))
    }
    scores <- c(
      pre_open = as.numeric(input$score_pre_open),
      sync_open = as.numeric(input$score_sync_open),
      sync_repress = as.numeric(input$score_sync_repress),
      unrelated = as.numeric(input$score_unrelated)
    )
    if (length(scores) != 4L || anyNA(scores)) {
      showNotification("评分输入未就绪，请稍候重试", type = "error")
      return(invisible(FALSE))
    }
    validate_scores(scores)
    save_annotation(
      con, r$candidate_id, scores,
      overall_confidence = as.numeric(input$confidence),
      review_status = NULL,
      notes = input$notes,
      annotator = input$annotator %||% cfg$app$annotator_default,
      app_version = app_version
    )
    show_evidence(TRUE)
    showNotification("已保存：四条评分 + 信心已写入数据库", type = "message")
    if (advance) go_to(1L)
    invisible(TRUE)
  }

  observeEvent(input$save_next, do_save(TRUE))
  observeEvent(input$save_only, do_save(FALSE))

  observeEvent(input$skip, {
    r <- current_row()
    if (!is.null(r)) {
      mark_skipped(con, r$candidate_id,
                   annotator = input$annotator %||% cfg$app$annotator_default,
                   app_version = app_version)
      showNotification("已暂时跳过（记录在数据库）", type = "message")
    }
    go_to(1L)
  })

  observeEvent(input$review, {
    r <- current_row()
    if (!is.null(r)) {
      mark_review(con, r$candidate_id,
                  annotator = input$annotator %||% cfg$app$annotator_default)
      showNotification("已标记稍后复核", type = "message")
    }
  })

  observeEvent(input$clear_scores, {
    updateSliderInput(session, "score_pre_open", value = 0)
    updateSliderInput(session, "score_sync_open", value = 0)
    updateSliderInput(session, "score_sync_repress", value = 0)
    updateSliderInput(session, "score_unrelated", value = 0)
    updateSliderInput(session, "confidence", value = 0)
    updateTextAreaInput(session, "notes", value = "")
  })

  observeEvent(input$toggle_evidence, {
    show_evidence(!show_evidence())
  })

  observeEvent(input$key_event, {
    switch(input$key_event,
           save_next = do_save(TRUE),
           prev = go_to(-1L),
           "next" = go_to(1L),
           skip = {
             r <- current_row()
             if (!is.null(r)) {
               mark_skipped(con, r$candidate_id,
                            annotator = input$annotator %||% cfg$app$annotator_default,
                            app_version = app_version)
             }
             go_to(1L)
           },
           review = {
             r <- current_row()
             if (!is.null(r)) {
               mark_review(con, r$candidate_id,
                           annotator = input$annotator %||% cfg$app$annotator_default)
             }
           },
           NULL)
  })

  # ---------------- Render outputs ----------------
  output$status_boxes <- renderUI({
    tryCatch({
      r <- current_row()
      st <- queue_state()
      if (is.null(r)) {
        return(div(class = "text-muted", "队列为空：请调整筛选条件。"))
      }
      ov <- overall_progress()
      badge <- function(label, value, strong = FALSE) {
        tags$span(
          class = "badge border text-dark px-2 py-1 me-1 mb-1",
          style = "font-size:13px; white-space:nowrap; background:#f8f9fa;",
          tags$span(class = "text-secondary", paste0(label, ": ")),
          if (strong) tags$strong(value) else value
        )
      }
      div(
        class = "d-flex flex-wrap align-items-center",
        badge("context", r$context),
        badge("peak", r$peak),
        badge("gene", r$gene),
        badge("distance",
              if (!is.na(r$distance)) paste0(format(r$distance, big.mark = ","), " bp") else "NA"),
        badge("mandatory", if (isTRUE(r$mandatory == 1)) "是" else "否",
              strong = isTRUE(r$mandatory == 1)),
        badge("来源", queue_source_label(r)),
        if (!is.na(r$review_status) && r$review_status == "review") {
          badge("状态", "待复核", strong = TRUE)
        } else if (!is.na(r$review_status) && r$review_status == "skipped") {
          badge("状态", "已跳过")
        },
        if (!is.na(r$batch_id)) badge("批次", r$batch_id),
        if (use_snapshot) {
          pf <- pair_features[candidate_id == r$candidate_id]
          if (nrow(pf)) {
            badge("数据", if (pf$is_organoid[1] == 1) "organoid" else "nonorganoid")
          }
        },
        tags$span(
          class = "ms-auto text-secondary",
          style = "font-size:13px; white-space:nowrap;",
          sprintf("当前筛选 %d / %d · 总体已标注 %d / %d",
                  st$pos, st$total, ov$done, ov$total)
        )
      )
    }, error = function(e) {
      div(class = "alert alert-danger",
          "界面渲染错误：", conditionMessage(e))
    })
  })

  output$evidence_panel <- renderUI({
    r <- current_row()
    if (is.null(r)) return(NULL)
    if (!show_evidence()) {
      return(NULL)
    }
    card(
      class = "border-info",
      card_header("统计证据（盲标模式默认隐藏）"),
      layout_columns(
        value_box(title = "S_pre_shape", value = if (!is.na(r$S_pre_shape)) round(r$S_pre_shape, 4) else "—"),
        value_box(title = "S_final", value = if (!is.na(r$S_final)) round(r$S_final, 4) else "—"),
        value_box(title = "来源标记",
                  value = paste0(
                    ifelse(isTRUE(r$source_tri_pre == 1), "tri_pre ", ""),
                    ifelse(isTRUE(r$source_sync_open == 1), "sync_open ", ""),
                    ifelse(isTRUE(r$source_sync_repress_heuristic == 1), "heuristic", "")
                  ))
      ),
      if (!is.na(r$repeat_of)) p("一致性重复样本，对应原始候选：", r$repeat_of)
    )
  })

  output$plot_peak <- renderPlotly({
    v <- current_values()
    if (is.null(v)) return(NULL)
    ggplotly(plot_raw_peak(v$peak, v$stages))
  })

  output$plot_gene <- renderPlotly({
    v <- current_values()
    if (is.null(v)) return(NULL)
    ggplotly(plot_raw_gene(v$gene, v$stages))
  })

  output$plot_overlay <- renderPlotly({
    v <- current_values()
    if (is.null(v)) return(NULL)
    ggplotly(plot_overlay(v$peak, v$gene, v$stages))
  })

  output$plot_deltas <- renderPlotly({
    v <- current_values()
    if (is.null(v)) return(NULL)
    ggplotly(plot_deltas(v$peak, v$gene, v$stages))
  })

  output$raw_table <- renderDT({
    v <- current_values()
    if (is.null(v)) return(NULL)
    t <- make_value_table(v$peak, v$gene, v$stages)
    datatable(t$raw, options = list(dom = "t", pageLength = 10),
              rownames = FALSE) |>
      formatRound(2:3, digits = 6)
  })

  output$delta_table <- renderDT({
    v <- current_values()
    if (is.null(v)) return(NULL)
    t <- make_value_table(v$peak, v$gene, v$stages)
    datatable(t$deltas, options = list(dom = "t", pageLength = 10),
              rownames = FALSE) |>
      formatRound(2:3, digits = 6)
  })

  overall_progress <- reactive({
    total_q <- count_candidates(con, list(), include_unqueued = FALSE)
    done <- count_annotated(con, list())
    list(total = total_q, done = done)
  })
}

shinyApp(ui = ui, server = server)
