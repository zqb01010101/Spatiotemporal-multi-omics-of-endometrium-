
library(ArchR)
library(Seurat)
library(BiocGenerics)
library(dplyr)
library(pheatmap)
library(SummarizedExperiment)
library(SeuratObject)
library(GenomicRanges)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(BSgenome)
library(ggseqlogo)
library(clusterProfiler)
library(ggseqlogo)
library(FNN)
library(Matrix)
library(future)
library(future.apply)
library(tidyr)


#### Helper functions ####
## Plot the distribution of a dataset
plot_distribution <- function(x,
                              type = c("hist", "violin"),
                              bins = NULL,
                              add_density = TRUE,
                              rug = TRUE,
                              na_rm = TRUE,
                              title = NULL) {
  stopifnot(is.numeric(x))
  
  # Preprocessing
  if (na_rm) x <- x[is.finite(x)]
  if (length(x) < 2) stop("向量长度不足（或全部为 NA/非有限数）。")
  
  # Select the number of bins automatically
  if (is.null(bins)) {
    bins <- tryCatch(nclass.FD(x), error = function(e) NA)
    if (is.na(bins) || bins < 5) bins <- max(5, ceiling(log2(length(x)) + 1))
  }
  
  # Hide the density curve for constant data
  if (sd(x) == 0) add_density <- FALSE
  
  type <- match.arg(type)
  library(ggplot2)
  df <- data.frame(value = x)
  
  if (type == "violin") {
    p <- ggplot(df, aes(x = "", y = value)) +
      geom_violin(trim = FALSE) +
      geom_boxplot(width = 0.12, outlier.size = 0.8) +
      xlab(NULL) + ylab("value") +
      ggtitle(ifelse(is.null(title), "Violin plot", title)) +
      theme_minimal(base_size = 12)
    if (rug) p <- p + geom_rug(aes(x = NULL, y = value), sides = "l", alpha = 0.4)
  } else { # hist
    p <- ggplot(df, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)), bins = bins, alpha = 0.75) +
      ggtitle(ifelse(is.null(title), "Histogram", title)) +
      xlab("value") + ylab("density") +
      theme_minimal(base_size = 12)
    if (add_density) p <- p + geom_density(adjust = 1.0, linewidth = 0.8)
    if (rug) p <- p + geom_rug(alpha = 0.4)
  }
  
  return(p)
}


get_tss_peaks_for_genes <- function(
    proj,
    genes,
    upstream = 1000,
    downstream = 1000
) {
  # 1. Get gene annotations
  ga <- getGeneAnnotation(proj)
  ga <- ga$genes
  m <- mcols(ga)
  
  # 2. Find the gene symbol column, preferring ga$symbol
  candidate_cols <- c("symbol", "gene_name", "name")
  symbol_col <- candidate_cols[candidate_cols %in% colnames(m)][1]
  if (is.na(symbol_col)) {
    stop("在 geneAnnotation 中找不到 symbol / gene_name / name 列，请先查看 mcols(getGeneAnnotation(proj))。")
  }
  gene_symbol <- m[[symbol_col]]
  
  # 3. Match gene symbols exactly, so NCAM excludes NCAM-AS1
  sel <- gene_symbol %in% genes
  if (!any(sel)) {
    stop("指定的基因都没有在 geneAnnotation 中找到，请检查基因名与物种/基因注释是否一致。")
  }
  
  ga_sub <- ga[sel]
  ga_sub$symbol <- gene_symbol[sel]   # Add an explicit symbol column for later use
  
  # 4. Define TSS windows
  # Use upstream = 0 and downstream = 1 for a 1 bp TSS window
  tss_windows <- promoters(ga_sub, upstream = upstream, downstream = downstream)
  
  # 5. Get peaks
  peaks <- getPeakSet(proj)
  
  # 6. Find peaks overlapping the TSS windows
  hits <- findOverlaps(peaks, tss_windows)
  if (length(hits) == 0L) {
    warning("没有任何 peak 落在指定的 TSS 窗口内，请尝试加大 upstream/downstream。")
    return(NULL)
  }
  
  peaks_hits <- peaks[queryHits(hits)]
  tss_hits   <- tss_windows[subjectHits(hits)]
  
  # 7. Build a data frame
  res <- data.frame(
    gene       = tss_hits$symbol,
    peak_chr   = as.character(seqnames(peaks_hits)),
    peak_start = start(peaks_hits),
    peak_end   = end(peaks_hits),
    peak_id    = paste0(seqnames(peaks_hits), ":", start(peaks_hits), "-", end(peaks_hits)),
    tss_chr    = as.character(seqnames(tss_hits)),
    tss_start  = start(tss_hits),
    tss_end    = end(tss_hits),
    strand     = as.character(strand(tss_hits)),
    stringsAsFactors = FALSE
  )
  
  # Optional: de-duplicate gene–peak_id pairs across overlapping transcript TSS windows
  res <- unique(res)
  
  return(res)
}



#### Correlate peak accessibility with gene expression
# Y_z: peaks by samples matrix with peak IDs as row names
#  <- : genes by samples matrix with gene symbols as row names
# peak_gene_mat: binary peak-by-gene candidate matrix, e.g. res_peak_gene_1Mb$mat
# use_cols: columns for correlation, e.g. stromal sec_* columns
# method: "pearson" or "spearman"
make_peak_gene_cor_matrix_subset <- function(Y_z,
                                             avg_list_binary,
                                             peak_gene_mat,
                                             use_cols,
                                             method = "pearson") {
  ## 1. Align peaks and genes
  peak_ids <- intersect(rownames(peak_gene_mat), rownames(Y_z))
  gene_ids <- intersect(colnames(peak_gene_mat), rownames(avg_list_binary))
  
  if (length(peak_ids) == 0) {
    stop("Y_z 与 peak_gene_mat 在行名（peak）上没有交集。")
  }
  if (length(gene_ids) == 0) {
    stop("avg_list_binary 与 peak_gene_mat 在基因名上没有交集。")
  }
  
  peak_gene_mat_sub <- peak_gene_mat[peak_ids, gene_ids, drop = FALSE]
  
  ## 2. Keep specified use_cols shared by both matrices
  common_cols <- intersect(use_cols, intersect(colnames(Y_z), colnames(avg_list_binary)))
  if (length(common_cols) < 2) {
    stop("用于计算的列太少（<2），请检查 use_cols 是否在矩阵列名中。")
  }
  
  Y_mat <- as.matrix(Y_z[peak_ids, common_cols, drop = FALSE])
  G_mat <- as.matrix(avg_list_binary[gene_ids, common_cols, drop = FALSE])
  
  ## 3. Correlate only 1 Mb candidate pairs (nonzero entries in peak_gene_mat_sub)
  m_t <- as(peak_gene_mat_sub, "dgTMatrix")
  n_pairs <- length(m_t@x)
  message("在列子集 {", paste(common_cols, collapse = ", "),
          "} 上计算候选 peak-gene 对数量：", n_pairs)
  
  cors <- numeric(n_pairs)
  for (k in seq_len(n_pairs)) {
    i <- m_t@i[k] + 1L  # row index (peak)
    j <- m_t@j[k] + 1L  # col index (gene)
    
    v_peak <- Y_mat[i, ]
    v_gene <- G_mat[j, ]
    
    cors[k] <- suppressWarnings(
      cor(as.numeric(v_peak),
          as.numeric(v_gene),
          use    = "pairwise.complete.obs",
          method = method)
    )
  }
  
  ## 4. Build a sparse correlation matrix with the structure of peak_gene_mat_sub
  PeakGene_cor_subset <- sparseMatrix(
    i = m_t@i + 1L,
    j = m_t@j + 1L,
    x = cors,
    dims     = dim(peak_gene_mat_sub),
    dimnames = dimnames(peak_gene_mat_sub)
  )
  
  message("相关性矩阵维度：",
          nrow(PeakGene_cor_subset), " peaks × ",
          ncol(PeakGene_cor_subset), " genes")
  
  return(PeakGene_cor_subset)
}







### Rescale residual r to downweight small values and upweight large ones, with f(0.2) = 0.1 and f(0.8) = 0.95
make_S_f_2anchors <- function(x1=0.2, y1=0.1, x2=0.8, y2=0.95) {
  stopifnot(0 < x1, x1 < x2, x2 < 1,
            0 < y1, y1 < y2, y2 < 1)
  
  # Normalize the logistic curve so f(0) = 0 and f(1) = 1
  f_norm <- function(x, a, k) {
    L <- function(z) 1 / (1 + exp(-k * (z - a)))
    (L(x) - L(0)) / (L(1) - L(0))
  }
  
  # Use a = plogis(t1) and k = exp(t2) to keep 0 < a < 1 and k > 0
  obj <- function(theta) {
    a <- plogis(theta[1])
    k <- exp(theta[2])
    (f_norm(x1, a, k) - y1)^2 + (f_norm(x2, a, k) - y2)^2
  }
  
  opt <- optim(par = c(qlogis(0.5), log(10)), fn = obj, method = "Nelder-Mead")
  a <- plogis(opt$par[1])
  k <- exp(opt$par[2])
  
  f <- function(x) {
    x <- pmin(pmax(x, 0), 1)  # Optional: clamp to [0, 1]
    f_norm(x, a, k)
  }
  attr(f, "params") <- c(a = a, k = k)
  f
}





## =========================================================
## 1) Extract candidate peak–gene pairs from the binary matrix
##    - Supports base matrices and sparse matrices such as Matrix::dgCMatrix
## =========================================================
get_candidate_pairs <- function(pair_mat) {
  if (!is.null(attr(pair_mat, "class")) && any(grepl("^dg|^ds|^dt|^Matrix", class(pair_mat)))) {
    # sparse Matrix
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("检测到稀疏矩阵，但未安装 Matrix 包。请先安装：install.packages('Matrix')")
    }
    s <- Matrix::summary(pair_mat)  # columns: i, j, x  (i/j are 1-based)
    s <- s[s$x != 0, , drop = FALSE]
    data.frame(
      peak = rownames(pair_mat)[s$i],
      gene = colnames(pair_mat)[s$j],
      stringsAsFactors = FALSE
    )
  } else {
    # base matrix
    idx <- which(pair_mat != 0, arr.ind = TRUE)
    data.frame(
      peak = rownames(pair_mat)[idx[, 1]],
      gene = colnames(pair_mat)[idx[, 2]],
      stringsAsFactors = FALSE
    )
  }
}

## =========================================================
## 2) Closed-form OLS fit and statistics for one (X, Y) pair
##    Model y = k*x + a
## =========================================================
fit_ols_1d <- function(x, y, min_points = 4) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  n <- length(x)
  
  # Require at least four observations by default
  if (n < min_points) {
    return(list(n = n, k = NA_real_, a = NA_real_, r2 = NA_real_, adj_r2 = NA_real_,
                rmse = NA_real_, cor = NA_real_, p_k = NA_real_, p_a = NA_real_, p_model = NA_real_))
  }
  
  xbar <- mean(x); ybar <- mean(y)
  xc <- x - xbar
  yc <- y - ybar
  Sxx <- sum(xc^2)
  Sxy <- sum(xc * yc)
  
  # The slope is undefined when x is constant (Sxx = 0)
  if (!is.finite(Sxx) || Sxx <= 0) {
    return(list(n = n, k = NA_real_, a = NA_real_, r2 = NA_real_, adj_r2 = NA_real_,
                rmse = NA_real_, cor = suppressWarnings(stats::cor(x, y)),
                p_k = NA_real_, p_a = NA_real_, p_model = NA_real_))
  }
  
  k <- Sxy / Sxx
  a <- ybar - k * xbar
  
  yhat <- k * x + a
  resid <- y - yhat
  SSE <- sum(resid^2)
  SST <- sum((y - ybar)^2)
  
  # R-squared, including the constant-y case
  if (SST <= 0) {
    r2 <- if (SSE <= 0) 1 else 0
  } else {
    r2 <- 1 - SSE / SST
  }
  
  # Adjusted R-squared for one predictor (p = 1)
  # adjR2 = 1 - (1-R2)*(n-1)/(n-p-1) = 1 - (1-R2)*(n-1)/(n-2)
  if (n > 2 && is.finite(r2)) {
    adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - 2)
  } else {
    adj_r2 <- NA_real_
  }
  
  rmse <- sqrt(SSE / n)
  cor_xy <- suppressWarnings(stats::cor(x, y))
  
  # Slope and intercept t-tests and overall F-test, with n > 2
  if (n > 2) {
    df <- n - 2
    MSE <- SSE / df
    
    se_k <- sqrt(MSE / Sxx)
    t_k  <- k / se_k
    p_k  <- 2 * stats::pt(-abs(t_k), df = df)
    
    se_a <- sqrt(MSE * (1 / n + (xbar^2) / Sxx))
    t_a  <- a / se_a
    p_a  <- 2 * stats::pt(-abs(t_a), df = df)
    
    # Overall F = (SSR/1)/(SSE/df), expressed through R-squared
    if (is.finite(r2) && r2 < 1) {
      F <- (r2 / (1 - r2)) * df
      p_model <- stats::pf(F, df1 = 1, df2 = df, lower.tail = FALSE)
    } else if (isTRUE(all.equal(r2, 1))) {
      p_model <- 0
    } else {
      p_model <- NA_real_
    }
  } else {
    p_k <- p_a <- p_model <- NA_real_
  }
  
  list(n = n, k = k, a = a, r2 = r2, adj_r2 = adj_r2,
       rmse = rmse, cor = cor_xy, p_k = p_k, p_a = p_a, p_model = p_model)
}

## =========================================================
## 3) Fit all candidate peak–gene pairs
## =========================================================
fit_peak_gene_pairs_ols <- function(
    Y_sig,
    avg_list_binary,
    time_order,
    pair_mat,
    min_points = 4,
    verbose = TRUE
) {
  # Align columns in time order
  cols <- intersect(time_order, intersect(colnames(Y_sig), colnames(avg_list_binary)))
  if (length(cols) < 3) stop("共同列太少：请检查 time_order 与两个矩阵的列名是否一致。")
  
  # Candidate pairs
  pairs <- get_candidate_pairs(pair_mat)
  if (nrow(pairs) == 0) {
    return(data.frame())
  }
  
  # Keep peaks and genes present in their matrices to avoid NA indices
  pairs <- pairs[pairs$peak %in% rownames(Y_sig) & pairs$gene %in% rownames(avg_list_binary), , drop = FALSE]
  if (verbose) message("Candidate pairs after name filtering: ", nrow(pairs))
  
  # Preallocate output
  m <- nrow(pairs)
  out <- data.frame(
    peak = pairs$peak,
    gene = pairs$gene,
    n = integer(m),
    k = numeric(m),
    a = numeric(m),
    r2 = numeric(m),
    adj_r2 = numeric(m),
    rmse = numeric(m),
    cor = numeric(m),
    p_k = numeric(m),
    p_a = numeric(m),
    p_model = numeric(m),
    stringsAsFactors = FALSE
  )
  out[,] <- NA  # Initialize all values to NA
  
  out$peak <- pairs$peak
  out$gene <-  pairs$gene
  
  # Fit each pair across six time points
  for (i in seq_len(m)) {
    pk <- pairs$peak[i]
    gn <- pairs$gene[i]
    x <- Y_sig[pk, cols]
    y <- avg_list_binary[gn, cols]
    
    fit <- fit_ols_1d(x, y, min_points = min_points)
    
    out$n[i]       <- fit$n
    out$k[i]       <- fit$k
    out$a[i]       <- fit$a
    out$r2[i]      <- fit$r2
    out$adj_r2[i]  <- fit$adj_r2
    out$rmse[i]    <- fit$rmse
    out$cor[i]     <- fit$cor
    out$p_k[i]     <- fit$p_k
    out$p_a[i]     <- fit$p_a
    out$p_model[i] <- fit$p_model
  }
  
  out
}




.row_max <- function(mat, na.rm = TRUE) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMaxs(mat, na.rm = na.rm)
  } else {
    apply(mat, 1, max, na.rm = na.rm)
  }
}
.row_min <- function(mat, na.rm = TRUE) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMins(mat, na.rm = na.rm)
  } else {
    apply(mat, 1, min, na.rm = na.rm)
  }
}

# NA-aware max for two matrices/vectors elementwise:
# if both NA -> NA; else max of available
.max2_na <- function(a, b) {
  aa <- a; bb <- b
  aa[is.na(aa)] <- -Inf
  bb[is.na(bb)] <- -Inf
  out <- pmax(aa, bb)
  out[is.infinite(out)] <- NA_real_
  out
}

# smoothstep on [0,1]
.smoothstep <- function(u) {
  u <- pmin(pmax(u, 0), 1)
  3*u^2 - 2*u^3
}

# New S-curve:
# f(0)=0, f(0.1)=0.1, f(0.5)=1, monotone, S-shaped, saturates at >=0.5
.f_scurve <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  out <- numeric(length(x))
  
  i1 <- x < 0.1
  if (any(i1)) {
    u <- x[i1] / 0.1
    out[i1] <- 0.1 * .smoothstep(u)
  }
  
  i2 <- (x >= 0.1) & (x < 0.5)
  if (any(i2)) {
    u <- (x[i2] - 0.1) / 0.4
    out[i2] <- 0.1 + 0.9 * .smoothstep(u)
  }
  
  i3 <- x >= 0.5
  if (any(i3)) out[i3] <- 1
  
  out
}


## Compute the final score
weighted_hmean <- function(Sv, Sw, w_v = 3, w_w = 1, eps = 1e-12) {
  ok <- is.finite(Sv) & is.finite(Sw) & Sv >= 0 & Sw >= 0
  out <- rep(NA_real_, length(Sv))
  denom <- w_v / pmax(Sv, eps) + w_w / pmax(Sw, eps)
  out[ok] <- (w_v + w_w) / denom[ok]
  out
}
## =========================
## Main: add two scores
##   1) S_warp (unchanged logic)
##   2) S_value = exp(-alpha_v * mean(S_value_dif_i) - beta_v * max(S_value_dif_i))
##      where S_value_dif_i = |G0_i - G1_i| * (1 + lambda * g_i^2) * f(G0_i / max(G0))
##      g_i from segment angle-diff delta_trend_t (normalized to [0,1]) with point max mapping
## =========================
add_scores_warp_value_v1 <- function(
    res_pairs2, Y_sig, avg_list_binary, time_order,
    baseline = c("none", "min"),
    eps = 1e-8,
    chunk_size = 50000L,
    
    ## ---- warp (same as v5)
    gamma_mean = 2, gamma_max = 1,
    w_plus = 1.2, w_minus = 1,
    
    ## ---- value score parameters
    lambda_trend = 1,
    alpha_v = 1, beta_v = 2,
    
    ## ---- guards against tiny ranges (highly recommended)
    min_range_peak_value = 0.1,
    min_range_gene_value = 0.1
) {
  baseline <- match.arg(baseline)
  
  need_cols <- c("peak", "gene", "k", "a")
  miss <- setdiff(need_cols, colnames(res_pairs2))
  if (length(miss) > 0) stop("res_pairs2 缺少必要列：", paste(miss, collapse = ", "))
  
  cols <- intersect(time_order, intersect(colnames(Y_sig), colnames(avg_list_binary)))
  
  Y6 <- as.matrix(Y_sig[, cols, drop = FALSE])
  G6 <- as.matrix(avg_list_binary[, cols, drop = FALSE])
  
  storage.mode(Y6) <- "double"
  storage.mode(G6) <- "double"
  
  peak_id <- match(res_pairs2$peak, rownames(Y6))
  gene_id <- match(res_pairs2$gene, rownames(G6))
  
  if (length(cols) < 3) stop("time_order 与矩阵共同列太少（<3）")
  Tn <- length(cols)
  
  n <- nrow(res_pairs2)
  ok0 <- if ("ok" %in% colnames(res_pairs2)) as.logical(res_pairs2$ok) else rep(TRUE, n)
  
  ok_pair <- ok0 &
    is.finite(res_pairs2$k) & is.finite(res_pairs2$a) &
    !is.na(peak_id) &
    !is.na(gene_id)
  
  ## outputs
  # value
  D_mean <- rep(NA_real_, n)
  D_max  <- rep(NA_real_, n)
  S_value <- rep(NA_real_, n)
  
  # diagnostics (optional but very useful)
  max_gi <- rep(NA_real_, n)
  max_gi_time <- rep(NA_character_, n)
  max_Sv <- rep(NA_real_, n)
  max_Sv_time <- rep(NA_character_, n)
  
  # warp
  warp_W_mean <- rep(NA_real_, n)
  warp_W_max  <- rep(NA_real_, n)
  S_warp      <- rep(NA_real_, n)
  
  idx_all <- which(ok_pair)
  if (length(idx_all) > 0) {
    for (st in seq(1, length(idx_all), by = chunk_size)) {
      ed <- min(length(idx_all), st + chunk_size - 1L)
      idx <- idx_all[st:ed]
      
      peaks_blk <- res_pairs2$peak[idx]
      genes_blk <- res_pairs2$gene[idx]
      k_blk <- res_pairs2$k[idx]
      a_blk <- res_pairs2$a[idx]
      
      P  <- Y6[peak_id[idx], , drop = FALSE]
      G0 <- G6[gene_id[idx], , drop = FALSE]
      
      if (baseline == "min") {
        P  <- P  - .row_min(P,  na.rm = TRUE)
        G0 <- G0 - .row_min(G0, na.rm = TRUE)
      }
      P[!is.finite(P)] <- NA_real_
      G0[!is.finite(G0)] <- NA_real_
      
      ## =========================
      ## (A) Build G1 = a + k*P
      ## =========================
      G1 <- sweep(P, 1, k_blk, "*")
      G1 <- sweep(G1, 1, a_blk, "+")
      
      MG0 <- .row_max(G0, na.rm = TRUE)
      MG1 <- .row_max(G1, na.rm = TRUE)
      MP <- .row_max(P, na.rm = TRUE)
      
      ## =========================
      ## (B) Trend angle diff per segment
      ##   rP_t = ΔP / range(P)
      ##   rG_t = ΔG / range(G0)
      ##   delta_trend_t = acos( (1+rP*rG)/(sqrt(1+rP^2)*sqrt(1+rG^2)) ) / pi
      ##   g_i (point) = max(adjacent delta_trend)   (your choice)
      ## =========================
      
      
      # RP <- .row_max(P,  na.rm = TRUE) - .row_min(P,  na.rm = TRUE)
      # RG <- .row_max(G0, na.rm = TRUE) - .row_min(G0, na.rm = TRUE)
      # ok_range <- is.finite(RP) & is.finite(RG) & (RP >= min_range_peak_value) & (RG >= min_range_gene_value)
      # 
      # dP <- P[, -1, drop = FALSE]  - P[, -Tn, drop = FALSE]   # n x (Tn-1)
      # dG <- G0[, -1, drop = FALSE] - G0[, -Tn, drop = FALSE]
      # 
      # rP <- sweep(dP, 1, pmax(RP, eps), "/")
      # rG <- sweep(dG, 1, pmax(RG, eps), "/")
      
      ## Use the G1 line
      RP <- .row_max(G1,  na.rm = TRUE) - .row_min(P,  na.rm = TRUE)
      RG <- .row_max(G0, na.rm = TRUE) - .row_min(G0, na.rm = TRUE)
      ok_range <- is.finite(RP) & is.finite(RG) & (RP >= min_range_peak_value) & (RG >= min_range_gene_value)
      
      dP <- G1[, -1, drop = FALSE]  - G1[, -Tn, drop = FALSE]   # n x (Tn-1)
      dG <- G0[, -1, drop = FALSE] - G0[, -Tn, drop = FALSE]
      
      rP <- sweep(dP, 1, pmax(pmax(MG0, MG1), eps), "/")
      rG <- sweep(dG, 1, pmax(pmax(MG0, MG1), eps), "/")
      
      
      
      
      
      
      # cos(theta)
      num <- 0.04 + rP * rG
      den <- sqrt(0.04 + rP^2) * sqrt(0.04 + rG^2)
      cosv <- num / pmax(den, eps)
      cosv <- pmin(pmax(cosv, -1), 1)
      
      delta_seg <- acos(cosv) / pi   # in [0,1]
      delta_seg[!is.finite(delta_seg)] <- NA_real_
      
      # point-level g_i using max
      g_point <- matrix(NA_real_, nrow = nrow(P), ncol = Tn)
      g_point[, 1]  <- delta_seg[, 1]
      g_point[, Tn] <- delta_seg[, Tn - 1]
      if (Tn > 2) {
        mid <- .max2_na(delta_seg[, 1:(Tn-2), drop = FALSE],
                        delta_seg[, 2:(Tn-1), drop = FALSE])
        g_point[, 2:(Tn-1)] <- mid
      }
      
      ## =========================
      ## (C) Value-difference per point:
      ##   Sv_i = |G0_i - G1_i| * (1 + lambda * g_i^2) * f(G0_i / max(G0))
      ## =========================
      # ---- normalize numeric difference by max(MG0, MG1) per pair
      # ---- new x_i = max(g0_i/M_G0 , g1_i/M_G1)
      
      den_num <- pmax(MG0, MG1)
      den_num[!is.finite(den_num) | den_num < eps] <- NA_real_
      
      delta_num <- abs(G0 - G1)
      delta_num <- sweep(delta_num, 1, pmax(den_num, eps), "/")
      
      
      x0 <- sweep(G0, 1, pmax(MG0, eps), "/")
      x1 <- sweep(G1, 1, pmax(MG1, eps), "/")
      
      # elementwise max with NA-safety (your helper .max2_na supports matrices too)
      x <- .max2_na(x0, x1)
      
      # clean & clamp to [0,1]
      x[!is.finite(x)] <- 0
      x <- pmin(pmax(x, 0), 1)
      
      # S-shaped weight
      f_x <- matrix(.f_scurve(as.vector(x)), nrow = nrow(x), ncol = ncol(x), byrow = FALSE)
      
      # Sv <- delta_num * (1 + lambda_trend * (g_point^2)) * f_x
      Sv <- delta_num * (1 + lambda_trend * g_point) * f_x
      Sv[!is.finite(Sv)] <- NA_real_
      
      # aggregate (mean + max)
      dmean <- rowMeans(Sv, na.rm = TRUE)
      dmax  <- .row_max(Sv, na.rm = TRUE)
      
      # apply range guard
      dmean[!ok_range] <- NA_real_
      dmax[!ok_range]  <- NA_real_
      
      D_mean[idx] <- dmean
      D_max[idx]  <- dmax
      S_value[idx] <- exp(-alpha_v * dmean - beta_v * dmax)
      
      # diagnostics: where g_i / Sv largest
      gi_max <- .row_max(g_point, na.rm = TRUE)
      max_gi[idx] <- gi_max
      arg_g <- max.col(ifelse(is.na(g_point), -Inf, g_point), ties.method = "first")
      max_gi_time[idx] <- cols[arg_g]
      
      Sv_max <- .row_max(Sv, na.rm = TRUE)
      max_Sv[idx] <- Sv_max
      arg_sv <- max.col(ifelse(is.na(Sv), -Inf, Sv), ties.method = "first")
      max_Sv_time[idx] <- cols[arg_sv]
      
      ## =========================
      ## (D) Warp score (unchanged from your v5 idea)
      ##   L1 = a + k*P ; Delta = L1 - P ; weighted_move = w_plus*Dpos + w_minus*Dneg
      ## =========================
      L1 <- G1  # same as a+kP
      m0 <- .row_max(P,  na.rm = TRUE)
      m1 <- .row_max(L1, na.rm = TRUE)
      den_amp <- pmax(m0, m1)
      den_amp[!is.finite(den_amp)] <- NA_real_
      
      Delta <- L1 - P
      Dpos <- pmax(Delta, 0)
      Dneg <- pmax(-Delta, 0)
      
      weighted_move <- w_plus * Dpos + w_minus * Dneg
      wm_mean <- rowMeans(weighted_move, na.rm = TRUE)
      wm_max  <- .row_max(weighted_move, na.rm = TRUE)
      
      warp_W_mean[idx] <- wm_mean / (den_amp + eps)
      warp_W_max[idx]  <- wm_max  / (den_amp + eps)
      
      message("Processed block: ", ed, " / ", length(idx_all))
    }
  }
  
  ok_sw <- is.finite(warp_W_mean) & is.finite(warp_W_max)
  S_warp[ok_sw] <- exp(-gamma_mean * warp_W_mean[ok_sw] - gamma_max * warp_W_max[ok_sw])
  
  out <- res_pairs2
  # value score
  out$D_mean_value <- D_mean
  out$D_max_value  <- D_max
  out$S_value      <- S_value
  # diagnostics
  out$max_gi <- max_gi
  out$max_gi_time <- max_gi_time
  out$max_Sv <- max_Sv
  out$max_Sv_time <- max_Sv_time
  
  # warp score
  out$warp_W_mean <- warp_W_mean
  out$warp_W_max  <- warp_W_max
  out$S_warp      <- S_warp
  
  # ---- Final overall score: harmonic mean of S_value and S_warp
  Sv <- out$S_value
  Sw <- out$S_warp
  
  ok_fin <- is.finite(Sv) & is.finite(Sw) & (Sv >= 0) & (Sw >= 0)
  S_final <- rep(NA_real_, nrow(out))
  denom <- Sv + Sw
  S_final[ok_fin & denom > 0] <- 2 * Sv[ok_fin & denom > 0] * Sw[ok_fin & denom > 0] / denom[ok_fin & denom > 0]
  S_final[ok_fin & denom == 0] <- 0
  
  out$S_final <- S_final
  
  
  out
}







####Script_41_get_peak_score_and_RNA_score####

## Extract the RNA expression matrix
cell_change_list <- list('sec_PO11_surface_immune_stroma' = 'sec_PO11_surface_stroma', 'sec_PO11_inner_immune_stroma' = 'sec_PO11_inner_stroma', 
                         'sec_PO11_immune_stroma' = 'sec_PO11_inner_stroma', 'sec_PO13_surface_immune_stroma' = 'sec_PO13_surface_stroma', 
                         'sec_PO13_inner_immune_stroma' = 'sec_PO13_inner_stroma')

new_Stroma_and_Gland_Seurat_for_ArchR <- subset(new_whole_sec_seurat_with_RIF_for_ArchR_True, subset = super_cell_type %in% c('stroma', 'epithelial'))

new_Stroma_and_Gland_Seurat_for_ArchR$Universal_name <- sapply(new_Stroma_and_Gland_Seurat_for_ArchR$orig.ident, function(x){
  a <- unlist(str_split(x, '_'))
  b <- paste0(a[1:3], collapse = '_')
  return(b)
})


new_Stroma_and_Gland_Seurat_for_ArchR$Universal_name <- paste0(new_Stroma_and_Gland_Seurat_for_ArchR$Universal_name, '#', new_Stroma_and_Gland_Seurat_for_ArchR$nA, 'x', new_Stroma_and_Gland_Seurat_for_ArchR$nB)
new_Stroma_and_Gland_Seurat_for_ArchR <- change_annotation(seurat_obj = new_Stroma_and_Gland_Seurat_for_ArchR, change_list = cell_change_list, old_cluster = 'period_exclusive_cell_type', new_cluster = 'period_exclusive_cell_type_2')

#### Extract the RNA matrix
new_Stroma_and_Gland_Seurat_for_ArchR <- SCTransform(new_Stroma_and_Gland_Seurat_for_ArchR)

DefaultAssay(new_Stroma_and_Gland_Seurat_for_ArchR) <- 'RNA'
new_Stroma_and_Gland_Seurat_for_ArchR <- JoinLayers(object = new_Stroma_and_Gland_Seurat_for_ArchR)
Expression_0_1_RNA_matrix <- new_Stroma_and_Gland_Seurat_for_ArchR@assays$RNA$counts
Expression_0_1_RNA_matrix@x[Expression_0_1_RNA_matrix@x > 0] <- 1
## Add the new assay
new_Stroma_and_Gland_Seurat_for_ArchR[["RNA_binary"]] <- CreateAssayObject(counts = Expression_0_1_RNA_matrix)


## New expression fractions by cell type
avg_list_binary <- AverageExpression(
  object  = new_Stroma_and_Gland_Seurat_for_ArchR,
  assays  = "RNA_binary",   
  slot    = "counts",  
  group.by = "period_super_cell_type"
)

avg_list_binary <- avg_list_binary$RNA_binary
colnames(avg_list_binary) <- gsub(pattern = '-', replacement = '_', x = colnames(avg_list_binary))






## ATAC analysis
unique(new_whole_sec_to_use_with_RIF_ArchR$period_exclusive_cell_type_4)
new_whole_sec_to_use_with_RIF_ArchR$period_exclusive_cell_type_2 <- new_whole_sec_to_use_with_RIF_ArchR$period_exclusive_cell_type_4




cells_use <- getCellNames(new_whole_sec_to_use_with_RIF_ArchR)[
  new_whole_sec_to_use_with_RIF_ArchR$super_cell_type %in% c('stroma', 'epithelial')
]


new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland <- subsetArchRProject(
  new_whole_sec_to_use_with_RIF_ArchR, cells = cells_use,
  outputDirectory = "ATAC/new_ArchR_stroma_and_gland",
  dropCells = TRUE, force = TRUE
)

saveArchRProject(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)



## Call peaks again for stromal and epithelial ArchR projects
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland <- addIterativeLSI(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_slim",   # Use a distinct name
  iterations = 2,
  dimsToUse = 2:25,
  force = TRUE
)

# 3. Compute group coverage separately for each peak call so the selected cells are present
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland <- addGroupCoverages(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland,
  groupBy = "period_super_cell_type",
  maxCells = 2000,
  minCells = 500,
  force = TRUE
)


new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland <- addReproduciblePeakSet(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland, 
  groupBy = "period_super_cell_type", 
  pathToMacs2 = '/home/liulab/anaconda3/envs/hjz_st_env/bin/macs2'
)

## Get peak scores at genomic positions
getPeakSet(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)

## Add peaks to the project
addPeakMatrix(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)

saveArchRProject(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)
gene_to_check_1 <- read.csv(file = 'gene_cluster/genes_both.csv')
gene_to_check_2 <- read.csv(file = 'gene_cluster/genes_only_time_change.csv')
moreExprInOrganoid <- readRDS("marker_gene/26_02_28_gene_more_expr_in_organoid.rds")


genes_of_interest  <- c(gene_to_check_1$gene, gene_to_check_2$gene)
genes_of_interest_include_Organoid <- unique(c(genes_of_interest, rownames(moreExprInOrganoid)))
gene_only_in_Organoid <- setdiff(genes_of_interest_include_Organoid, genes_of_interest)


## Get gene annotations and the peak set
ga <- getGeneAnnotation(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)$genes          # GRanges of gene annotations
pk <- getPeakSet(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland)                       # GRanges of the current peak set

symCol <- 'symbol'

ga.sel <- ga[mcols(ga)[[symCol]] %in% genes_of_interest_include_Organoid] 
missing_genes <- setdiff(genes_of_interest_include_Organoid, unique(mcols(ga.sel)[[symCol]]))
missing_genes_another_form <- norm_sym(missing_genes)

genes_of_interest_include_Organoid <- c(genes_of_interest_include_Organoid, missing_genes_another_form)
ga.sel <- ga[mcols(ga)[[symCol]] %in% genes_of_interest_include_Organoid]
missing_genes <- setdiff(genes_of_interest_include_Organoid, unique(mcols(ga.sel)[[symCol]]))






## Extend each gene body by 1,000 kb on either side
updown_kb <- 1000                                         

ga.ext <- GRanges(
  seqnames = seqnames(ga.sel),
  ranges   = IRanges(start = pmax(1, start(ga.sel) - updown_kb*1000L),
                     end   = end(ga.sel) + updown_kb*1000L),
  strand   = strand(ga.sel)
)
ga.ext <- trim(ga.ext)  # Trim to chromosome boundaries

## Keep overlapping peaks
peaks.keep <- subsetByOverlaps(pk, ga.ext, ignore.strand = TRUE)
peaks.keep <- sort(unique(peaks.keep))

message("原始 Peak 数量: ", length(pk))
message("筛选后 Peak 数量: ", length(peaks.keep))




## Create a new ArchR project
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addPeakSet(ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland, peakSet = peaks.keep, force = TRUE)

# Recompute PeakMatrix for peaks.keep
# Set binarize = TRUE to binarize insertion counts
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addPeakMatrix(ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered, force = TRUE, binarize = FALSE)

# Save as a new project, preferably in a new directory
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- saveArchRProject(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  outputDirectory = "ATAC/new_ArchR_stroma_and_gland_Filtered/"
)



## Recompute dimensionality reduction
new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addIterativeLSI(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  useMatrix = "PeakMatrix",
  name = "IterativeLSI_slim",   # Use a distinct name
  iterations = 2,
  dimsToUse = 2:25,
  force = TRUE
)

new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addClusters(
  input = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  reducedDims = "IterativeLSI_slim",
  name = "Clusters_slim",
  resolution = 0.6,
  force = TRUE
)


new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addCoAccessibility(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  reducedDims = "IterativeLSI_slim",
  k = 50,
  dimsToUse  = 2:25,
  corCutOff = 0.5,
  maxDist = 250000, 
)


new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered <- addPeak2GeneLinks(
  ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  reducedDims = "IterativeLSI_slim",
  useMatrix = "GeneExpressionMatrix",  # Previously added expression matrix
  corCutOff = 0.45,                    # Adjust as needed
  dimsToUse = 2:25, 
)


gm_tsse2_new <- make_group_peak_mats_cell_TSSE(
  proj       = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  group_col  = "period_super_cell_type",
  sample_col = "orig.ident",
  min_cells_per_sg = 30,
  scale_by   = 'ReadsInTSS',   # Scale by ReadsInTSS
  relative   = TRUE,       # Set TRUE to scale by ReadsInTSS / median(ReadsInTSS)
  nfrags_floor = 100,
  winsor = 0
)


# Convert to numeric
Y <- gm_tsse2_new$peak_group_mean_scaled
Y_sample <- gm_tsse2_new$Y_sg_mean_scaled

eps <- min(Y[Y > 0], na.rm = TRUE) / 100   # Use a pseudocount scaled to the data: 1/100 of the smallest positive value

Y_log <- log10(Y + eps)
Y_sample_log <- log10(Y_sample + eps)

x_nonzero <- as.numeric(Y_log[Y > 0])    # Fit N(mu, sigma^2) using nonzero values only
x_sample_nonzero <- as.numeric(Y_sample_log[Y_sample > 0])

mu    <- mean(x_nonzero)
sigma <- sd(x_nonzero)
Y_z <- (Y_log - mu) / sigma     # Z-score all entries, including zeros


mu_sample    <- mean(x_sample_nonzero)
sigma_sample <- sd(x_sample_nonzero)
Y_z_sample <- (Y_sample_log - mu_sample) / sigma_sample     # Z-score all entries, including zeros



##
group_order_no_RIF <- c("sec_PO3_stroma","sec_PO5_stroma","sec_PO7_stroma","sec_PO9_stroma","sec_PO11_stroma","sec_PO13_stroma",
                        "sec_PO3_epithelial","sec_PO5_epithelial","sec_PO7_epithelial","sec_PO9_epithelial","sec_PO11_epithelial","sec_PO13_epithelial")


## Find peaks near transcription start sites
tss_peaks <- get_tss_peaks_for_genes(
  proj    = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  genes   = genes_of_interest_include_Organoid,
  upstream = 2000,   # 2 kb upstream of the TSS
  downstream = 2000  # 2 kb downstream of the TSS
)

## Compute correlations by cell type, then select one TSS peak per gene
cor_res2 <- score_tss_peaks_by_cor_twoCT(
  tss_peaks = tss_peaks,
  Y_z       = Y_z,
  avg_expr  = avg_list_binary,
  method    = "pearson",
  use_abs_cor = FALSE,
  high_peak_cutoff = 1.8,
  high_gene_cutoff = 0.3,
  low_peak_cutoff  = -0.3,
  low_gene_cutoff  = 0.1
)


tss_peaks_scored2 <- cor_res2$all
best_tss_peak_gene2 <- cor_res2$best   # One TSS peak per gene



plot_res2 <- plot_peak_gene_scatter_twoCT(
  score_df    = best_tss_peak_gene2,
  Y_z         = Y_z,
  avg_expr    = avg_list_binary,
  cor_cutoff  = 0.8,
  use_abs_cor = FALSE
)


print(plot_res2$scatter)
View(plot_res2$data_long)

aaa <- plot_res2$data_long
View(aaa[aaa$peak_score < 0.1 & aaa$gene_expr > 0.5,])

#### Fit the data with a prior
df_fit <- plot_res2$data_long   # Long-format table used for the scatter plot

# # Fit only observations with expression above 0.1
# df_fit <- df_fit[df_fit$gene_expr > 0.1, ]

# Drop NA values
df_fit <- df_fit[is.finite(df_fit$peak_score) & is.finite(df_fit$gene_expr), ]

x_all <- df_fit$peak_score
y_all <- df_fit$gene_expr

# Clamp y to the open interval (0, 1) for a finite logit
eps <- 1e-5
y_all <- pmin(pmax(y_all, eps), 1 - eps)

xa <- 0
xb <- 1.7   # Set which point to constrain


obj_raw <- function(theta, x, y) {
  # Enforce the hard constraint by reparameterization
  p0  <- 0.05 * plogis(theta[1])          # Set the constraint value
  p18 <- 0.75 + 0.25 * plogis(theta[2])     # f(1.8) ∈ (0.8, 1) => f(1.8) >= 0.8
  
  # Derive x0 and k from (xa, p0) and (xb, p18)
  L0  <- qlogis(p0)
  L18 <- qlogis(p18)
  
  k  <- (xb - xa) / (L18 - L0)
  x0 <- xa - k * L0
  
  y_hat <- plogis((x - x0) / k)
  sum((y_hat - y)^2)
}

opt <- optim(
  par    = c(0, 0),
  fn     = obj_raw,
  x      = x_all,
  y      = y_all,
  method = "BFGS"
)

theta_hat <- opt$par

# Convert fitted theta back to f(0) and f(1.8)
p0_hat  <- 0.05 * plogis(theta_hat[1])
p18_hat <- 0.75 + 0.25 * plogis(theta_hat[2])   ## Update this value too when changing the constraint

L0  <- qlogis(p0_hat)
L18 <- qlogis(p18_hat)

k_hat  <- (xb - xa) / (L18 - L0)
x0_hat <- xa - k_hat * L0

cat("f(0)   ≈", p0_hat,  "\n")
cat("f(1.7) ≈", p18_hat, "\n")
cat("x0 =", x0_hat, "  k =", k_hat, "\n")

sigmoid_fit <- function(x) plogis((x - x0_hat) / k_hat)



x_grid <- seq(min(x_all), max(x_all), length.out = 400)
line_df <- data.frame(
  peak_score    = x_grid,
  gene_expr_hat = sigmoid_fit(x_grid)
)

library(ggplot2)
ggplot(df_fit, aes(peak_score, gene_expr, color = cell_type)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_line(data = line_df,
            aes(peak_score, gene_expr_hat),
            inherit.aes = FALSE,
            color = "black",
            size = 1.2) +
  theme_bw() +
  labs(
    x = "Peak accessibility (Y_z)",
    y = "Gene expression (avg_list_binary)",
    color = "Cell type",
    title = "Constrained sigmoid fit (f(0)∈(0,0.05), f(2)∈(0.95,1))"
  )



#### Apply the fitted function to all of Y_z

## Convert sparse matrices such as dgCMatrix to a dense matrix or operate on @x directly
Y_num <- as.matrix(Y_z)

## Allocate an output matrix of the same size
Y_sig <- matrix(
  NA_real_,
  nrow = nrow(Y_num),
  ncol = ncol(Y_num),
  dimnames = dimnames(Y_num)
)

## Transform finite values only
idx <- is.finite(Y_num)
Y_sig[idx] <- sigmoid_fit(Y_num[idx])


## Sample-level results
Y_sample_num <- as.matrix(Y_z_sample)

## Allocate an output matrix of the same size
Y_sample_sig <- matrix(
  NA_real_,
  nrow = nrow(Y_sample_num),
  ncol = ncol(Y_sample_num),
  dimnames = dimnames(Y_sample_num)
)

## Transform finite values only
idx <- is.finite(Y_sample_num)
Y_sample_sig[idx] <- sigmoid_fit(Y_sample_num[idx])




















####Script_42_get_Sync-open_scores_Using_Statistical_Methods####
pair_mat <- res_peak_gene_1Mb$mat
res1_stroma <- fit_peak_gene_pairs_ols(
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  time_order = sec_stroma_cols,
  pair_mat = res_peak_gene_1Mb$mat,
  min_points = 5
)
res1_epi <- fit_peak_gene_pairs_ols(
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  time_order = sec_epi_cols,
  pair_mat = res_peak_gene_1Mb$mat,
  min_points = 5
)
head(res1[order(-res1$r2), ], 20)


res1[res1$peak == 'chr7:38178556-38179056' & res1$gene == 'VPS41', ]




res_pairs2_scored_stroma <- add_scores_warp_value_v1(
  res_pairs2 = res1_stroma,
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  time_order = sec_stroma_cols,
  baseline = "none",
  chunk_size = 50000,
  
  # warp
  gamma_mean = 2, gamma_max = 1,
  w_plus = 1.2, w_minus = 1,
  
  # value score
  lambda_trend = 2,
  alpha_v = 1, beta_v = 2,
  
  # Keep this lower bound to stabilize r = delta/range for small ranges
  min_range_peak_value = 0.08,
  min_range_gene_value = 0.08
)

# example
res_pairs2_scored_stroma$S_final <- weighted_hmean(
  res_pairs2_scored_stroma$S_value,
  res_pairs2_scored_stroma$S_warp,
  w_v = 5, w_w = 1
)

# Inspect the two scores separately
# 1) High value agreement: high S_value
# 2) Little warping: high S_warp
res_pairs2_scored_stroma <- res_pairs2_scored_stroma[order(-res_pairs2_scored_stroma$S_final), ]
res_pairs2_scored_stroma <- res_pairs2_scored_stroma[res_pairs2_scored_stroma$k >0.45 & res_pairs2_scored_stroma$k < 2.1, ]
res_pairs2_scored_stroma <- res_pairs2_scored_stroma[!is.na(res_pairs2_scored_stroma$S_final),]
head(res_pairs2_scored_stroma[, c("peak","gene","S_value","S_warp","max_gi_time","max_Sv_time")])

plot_distribution(x = res_pairs2_scored_stroma$S_final[res_pairs2_scored_stroma$S_final > 0.0])


pdf('ATAC_plot/26_06_09_stroma_sync_SyncScore_1.pdf', height = 8, width = 9)
res_pairs2_scored_stroma <- res_pairs2_scored_stroma[order(-res_pairs2_scored_stroma$S_final), ]
for(i in seq(5000, 8000)){
  plot_sync_pair(
    peak = res_pairs2_scored_stroma$peak[i],
    gene = res_pairs2_scored_stroma$gene[i],
    res_pairs2 = res_pairs2_scored_stroma,
    Y_sig = Y_sig,
    avg_list_binary = avg_list_binary,
    time_order = sec_stroma_cols,
    scale_to_01 = FALSE,
    ylim_from_zero = TRUE
  )
}
dev.off()



res_pairs2_scored_stroma[res_pairs2_scored_stroma$peak == 'chr8:132877234-132877734' & res_pairs2_scored_stroma$gene == 'SLA', ]
res_pairs2_scored_stroma[res_pairs2_scored_stroma$peak == 'chr1:12474975-12475475' & res_pairs2_scored_stroma$gene == 'DHRS3', ]
res_pairs2_scored_stroma[res_pairs2_scored_stroma$peak == 'chr1:231590520-231591020' & res_pairs2_scored_stroma$gene == 'DISC1', ]









### Glandular epithelium


res_pairs2_scored_epi <- add_scores_warp_value_v1(
  res_pairs2 = res1_epi,
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  time_order = sec_epi_cols,
  baseline = "none",
  chunk_size = 50000,
  
  # warp
  gamma_mean = 2, gamma_max = 1,
  w_plus = 1.2, w_minus = 1,
  
  # value score
  lambda_trend = 2,
  alpha_v = 1, beta_v = 2,
  
  # Keep this lower bound to stabilize r = delta/range for small ranges
  min_range_peak_value = 0.08,
  min_range_gene_value = 0.08
)

# example
res_pairs2_scored_epi$S_final <- weighted_hmean(
  res_pairs2_scored_epi$S_value,
  res_pairs2_scored_epi$S_warp,
  w_v = 5, w_w = 1
)

# Inspect the two scores separately
# 1) High value agreement: high S_value
# 2) Little warping: high S_warp
res_pairs2_scored_epi <- res_pairs2_scored_epi[order(-res_pairs2_scored_epi$S_final), ]
res_pairs2_scored_epi <- res_pairs2_scored_epi[res_pairs2_scored_epi$k >0.45 & res_pairs2_scored_epi$k < 2.1, ]
head(res_pairs2_scored_epi[, c("peak","gene","S_value","S_warp","max_gi_time","max_Sv_time")])

plot_distribution(x = res_pairs2_scored_epi$S_final[res_pairs2_scored_epi$S_final > 0.0])


pdf('ATAC_plot/26_06_09_epi_sync_SyncScore_1.pdf', height = 8, width = 9)
res_pairs2_scored_epi <- res_pairs2_scored_epi[order(-res_pairs2_scored_epi$S_final), ]
for(i in seq(2000, 5000)){
  plot_sync_pair(
    peak = res_pairs2_scored_epi$peak[i],
    gene = res_pairs2_scored_epi$gene[i],
    res_pairs2 = res_pairs2_scored_epi,
    Y_sig = Y_sig,
    avg_list_binary = avg_list_binary,
    time_order = sec_epi_cols,
    scale_to_01 = FALSE,
    ylim_from_zero = TRUE
  )
}
dev.off()

