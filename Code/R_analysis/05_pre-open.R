#### Helper functions ####

### Build a binary peak-by-gene matrix with 1 Mb flanks
# proj: ArchRProj object
# gene_symbols: candidate gene symbols, e.g. ga.sel$symbol
# half_window: flanking bases on each side (default 1,000 kb, for 2 Mb of flanking sequence)
make_peak_gene_window_matrix <- function(proj,
                                         gene_symbols,
                                         half_window = 1000000,
                                         gene_rep = NULL,
                                         verbose = TRUE) {
  ## 1. Get peak GRanges for PeakMatrix
  peak_gr <- getPeakSet(proj)
  
  names(peak_gr) <- paste0(seqnames(peak_gr), ":", start(peak_gr), "-", end(peak_gr))
  
  ## 2. Extract the genes element from the gene annotation
  ga <- getGeneAnnotation(proj)
  gene_anno_all <- ga$genes   # This element is a GRanges object
  
  # # Drop entries with missing symbols
  # gene_anno_all <- gene_anno_all[!is.na(gene_anno_all$symbol)]
  # 
  # # Keep the first entry for each duplicate symbol by default
  # keep_idx <- !duplicated(gene_anno_all$symbol)
  # if (any(!keep_idx)) {
  #   dup_syms <- unique(gene_anno_all$symbol[!keep_idx])
  #   message("Duplicate gene symbols in the annotation. Keeping the first entry by default: ",
  #           paste(dup_syms, collapse = ", "))
  #   gene_anno_all <- gene_anno_all[keep_idx]
  # }
  # 
  
  # Use the 5'-most TSS to represent each symbol
  if(is.null(gene_rep)){
    gene_rep <- get_representative_tss_by_symbol(gene_anno_all, verbose = verbose)
  }
  
  
  ## 3. Filter candidate genes in ga.sel$symbol order
  gene_order <- gene_symbols[gene_symbols %in% gene_anno_all$symbol]
  if (length(gene_order) == 0) {
    stop("候选基因在基因注释中一个都没找到，请检查基因名是否一致（大小写等）。")
  }
  
  genes_not_found <- setdiff(gene_symbols, gene_order)
  if (verbose && length(genes_not_found) > 0) {
    message("这些候选基因在基因注释中未找到，将被忽略： ",
            paste(genes_not_found, collapse = ", "))
  }
  
  gene_anno <- gene_rep[gene_order]
  
  # gene_anno <- gene_anno_all[match(gene_order, gene_anno_all$symbol)]
  # names(gene_anno) <- gene_anno$symbol
  
  ## 4. Extend each gene on both sides by half_window (500 kb per side gives 1 Mb of flanking sequence)
  # Extend the full gene span. Use ga$TSS for TSS-centered windows.
  # genes_win <- gene_anno
  # start(genes_win) <- pmax(1L, start(gene_anno) - half_window)
  # end(genes_win)   <- end(gene_anno) + half_window
  
  
  # gene_anno contains the GRanges selected by symbol
  tss_pos <- ifelse(as.character(strand(gene_anno)) == "-", end(gene_anno), start(gene_anno))
  
  tss_gr <- GRanges(
    seqnames = seqnames(gene_anno),
    ranges   = IRanges(start = tss_pos, end = tss_pos),
    strand   = strand(gene_anno)
  )
  
  tss_gr <- resize(gene_anno, width = 1, fix = "start")  # Default: ignore.strand = FALSE
  
  
  genes_win <- tss_gr
  start(genes_win) <- pmax(1L, start(tss_gr) - half_window)
  end(genes_win)   <- end(tss_gr) + half_window
  
  
  # Clip intervals to chromosome boundaries
  # This affects only genes near chromosome ends
  genes_win <- trim(genes_win)
  
  ## 5. Find peaks overlapping gene windows
  hits <- findOverlaps(peak_gr, genes_win, ignore.strand = TRUE)
  message("findOverlaps 得到 ", length(hits), " 对 peak-gene 在窗口内重叠。")
  
  ## 6. Build a sparse binary peak-by-gene matrix
  n_peaks <- length(peak_gr)
  n_genes <- length(genes_win)
  
  PeakGene_mat <- sparseMatrix(
    i = as.integer(queryHits(hits)),
    j = as.integer(subjectHits(hits)),
    x = rep(1L, length(hits)),
    dims = c(n_peaks, n_genes),
    dimnames = list(
      peak = names(peak_gr),
      gene = names(genes_win)
    )
  )
  # Keep entries binary even if a pair is counted multiple times
  PeakGene_mat@x[] <- 1L
  
  message("矩阵维度：", nrow(PeakGene_mat), " peaks × ", ncol(PeakGene_mat), " genes")
  
  return(list(
    mat      = PeakGene_mat,  # Binary peak-by-gene matrix
    peaks_gr = peak_gr,       # Peak GRanges
    genes_gr = genes_win,     # Gene window GRanges
    tss_gr = tss_gr           # Also return the representative 1 bp TSS
  ))
}




### Build a peak-to-TSS distance matrix
# proj: ArchRProj object
# peak_gene_mat: binary peak-by-gene candidate matrix, e.g. res_peak_gene_1Mb$mat
make_peak_gene_tss_distance_matrix <- function(proj,
                                               peak_gene_mat,
                                               gene_rep = NULL,
                                               verbose = TRUE) {
  
  ## 1) Peak centers named consistently with the window matrix
  peak_gr <- getPeakSet(proj)
  names(peak_gr) <- paste0(seqnames(peak_gr), ":", start(peak_gr), "-", end(peak_gr))
  peak_center_gr <- GenomicRanges::resize(peak_gr, width = 1, fix = "center")
  names(peak_center_gr) <- names(peak_gr)
  
  peak_ids <- rownames(peak_gene_mat)
  if (is.null(peak_ids)) stop("peak_gene_mat 必须有行名（peak ID）。")
  
  peak_idx <- match(peak_ids, names(peak_center_gr))
  if (anyNA(peak_idx)) stop("peak_gene_mat 行名无法匹配 getPeakSet(proj) 的 peak 命名。")
  peak_center_gr2 <- peak_center_gr[peak_idx]
  peak_chr <- as.character(seqnames(peak_center_gr2))
  peak_pos <- start(peak_center_gr2)  # 1bp
  
  ## 2) Representative genes using the same TSS selection rule
  gene_ids <- colnames(peak_gene_mat)
  if (is.null(gene_ids)) stop("peak_gene_mat 必须有列名（gene symbol）。")
  
  ga <- getGeneAnnotation(proj)
  gene_anno_all <- ga$genes
  if (is.null(gene_rep)) {
    gene_rep <- get_representative_tss_by_symbol(gene_anno_all, verbose = verbose)
  }
  
  gene_idx <- match(gene_ids, names(gene_rep))
  keep_gene <- !is.na(gene_idx)
  if (any(!keep_gene)) {
    warning("这些基因在 gene_rep 中找不到，将被丢弃：\n",
            paste(gene_ids[!keep_gene], collapse = ", "))
  }
  
  peak_gene_mat2 <- peak_gene_mat[, keep_gene, drop = FALSE]
  gene_ids_valid <- colnames(peak_gene_mat2)
  gene_anno <- gene_rep[gene_ids_valid]
  
  if (is.null(mcols(gene_anno)$tss_pos)) stop("gene_rep 缺少 tss_pos 列。")
  
  tss_chr <- as.character(seqnames(gene_anno))
  tss_pos <- as.integer(mcols(gene_anno)$tss_pos)
  
  ## 3) Compute distances only for candidate pairs using Matrix::summary
  sm <- Matrix::summary(peak_gene_mat2)  # i/j 1-based
  i <- sm$i
  j <- sm$j
  
  same_chr <- peak_chr[i] == tss_chr[j]
  dist_vec <- abs(peak_pos[i] - tss_pos[j])
  dist_vec[!same_chr] <- NA_real_
  
  ## 4) Sparse distance matrix
  PeakGene_TSSdist <- Matrix::sparseMatrix(
    i = i, j = j, x = as.numeric(dist_vec),
    dims = dim(peak_gene_mat2),
    dimnames = list(peak = rownames(peak_gene_mat2), gene = gene_ids_valid)
  )
  
  message("TSS 距离矩阵维度：", nrow(PeakGene_TSSdist), " peaks × ", ncol(PeakGene_TSSdist), " genes")
  
  list(
    dist_mat       = PeakGene_TSSdist,
    peak_center_gr = peak_center_gr2
  )
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



run_preopen_bayes_fast <- function(Y_sig,
                                   avg_list_binary,
                                   res_mat,
                                   time_order = colnames(Y_sig),
                                   link_alt_min_prob = 0.5,
                                   pair_chunk_size = 200000L,
                                   shape_chunk_size = 100000L,
                                   shape_only_related = TRUE,
                                   mix_fit_n = 100000L,
                                   mix_maxit = 1500L,
                                   seed = 1L,
                                   atac_range_min = 0.1,
                                   expr_range_min = 0.1,
                                   atac_min_max = 0.2,
                                   atac_max_min = 0.7) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("需要安装 Matrix 包。")
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("需要安装 matrixStats 包。")
  }
  if (!requireNamespace("mixtools", quietly = TRUE)) {
    stop("需要安装 mixtools 包。")
  }
  
  stopifnot(all(time_order %in% colnames(Y_sig)))
  stopifnot(all(time_order %in% colnames(avg_list_binary)))
  
  Y_sig <- Y_sig[, time_order, drop = FALSE]
  avg_list_binary <- avg_list_binary[, time_order, drop = FALSE]
  
  common_peaks <- intersect(rownames(Y_sig), rownames(res_mat))
  common_genes <- intersect(colnames(res_mat), rownames(avg_list_binary))
  
  if (length(common_peaks) == 0L) stop("Y_sig 与 res_mat 没有共同 peak。")
  if (length(common_genes) == 0L) stop("avg_list_binary 与 res_mat 没有共同 gene。")
  
  Y_sig_sub <- as.matrix(Y_sig[common_peaks, , drop = FALSE])
  expr_sub <- as.matrix(avg_list_binary[common_genes, , drop = FALSE])
  res_sub <- res_mat[common_peaks, common_genes, drop = FALSE]
  
  message("Step 0: extracting sparse candidate indices...")
  pair_index <- .extract_candidate_pairs(res_sub)
  peak_idx <- pair_index$i
  gene_idx <- pair_index$j
  n_pairs <- length(peak_idx)
  message("Candidate pairs: ", format(n_pairs, big.mark = ","))
  
  message("Step 1: computing correlations by chunked matrix operations...")
  link <- .compute_link_stats_fast(
    peak_idx = peak_idx,
    gene_idx = gene_idx,
    Y_sig_sub = Y_sig_sub,
    expr_sub = expr_sub,
    chunk_size = pair_chunk_size,
    atac_range_min = atac_range_min,
    expr_range_min = expr_range_min,
    atac_min_max = atac_min_max,
    atac_max_min = atac_max_min
  )
  
  r_vec <- link$S_link
  keep <- is.finite(r_vec)
  if (sum(keep) < 20L) {
    stop("可用的 S_link 太少，没法拟合混合模型。")
  }
  
  r_use <- pmax(pmin(r_vec[keep], 0.999), -0.999)
  
  if (!is.null(seed)) set.seed(seed)
  if (is.finite(mix_fit_n) && length(r_use) > mix_fit_n) {
    fit_id <- sample.int(length(r_use), as.integer(mix_fit_n))
    r_fit <- r_use[fit_id]
    message(
      "Step 2: fitting 3-component mixture on a sample of ",
      format(length(r_fit), big.mark = ","),
      " correlations..."
    )
  } else {
    r_fit <- r_use
    message(
      "Step 2: fitting 3-component mixture on all ",
      format(length(r_fit), big.mark = ","),
      " correlations..."
    )
  }
  
  fit_link <- mixtools::normalmixEM(
    r_fit,
    k = 3,
    maxit = as.integer(mix_maxit),
    verb = FALSE
  )
  
  posterior_keep <- .posterior_from_normalmix(r_use, fit_link)
  
  mu_link <- fit_link$mu
  pos_comp <- which.max(mu_link)
  remaining <- setdiff(seq_along(mu_link), pos_comp)
  null_comp <- remaining[which.min(abs(mu_link[remaining]))]
  neg_comp <- setdiff(seq_along(mu_link), c(pos_comp, null_comp))
  
  p_pos <- numeric(n_pairs)
  p_unrel <- rep(1, n_pairs)
  p_pos[keep] <- posterior_keep[, pos_comp]
  p_unrel[keep] <- matrixStats::rowSums2(
    posterior_keep[, c(null_comp, neg_comp), drop = FALSE]
  )
  
  shape_idx <- if (shape_only_related) {
    which(keep & p_pos >= link_alt_min_prob)
  } else {
    seq_len(n_pairs)
  }
  
  message(
    "Step 3: computing shape statistics for ",
    format(length(shape_idx), big.mark = ","),
    " pairs..."
  )
  
  S_pre_shape <- rep(NA_real_, n_pairs)
  pre_peak_pos <- rep(NA_integer_, n_pairs)
  pre_baseline <- rep(NA_real_, n_pairs)
  first_open <- rep(NA_integer_, n_pairs)
  
  RNA_pre_max_frac <- rep(NA_real_, n_pairs)
  RNA_pre_min_frac <- rep(NA_real_, n_pairs)
  ATAC_pre_max_frac <- rep(NA_real_, n_pairs)
  ATAC_pre_min_frac <- rep(NA_real_, n_pairs)
  RNA_post_max_frac <- rep(NA_real_, n_pairs)
  RNA_post_min_frac <- rep(NA_real_, n_pairs)
  ATAC_post_max_frac <- rep(NA_real_, n_pairs)
  ATAC_post_min_frac <- rep(NA_real_, n_pairs)
  resid_mat <- matrix(NA_real_, nrow = n_pairs, ncol = ncol(Y_sig_sub))
  
  if (length(shape_idx) > 0L) {
    shape <- .compute_pre_shape_fast(
      peak_idx = peak_idx[shape_idx],
      gene_idx = gene_idx[shape_idx],
      Y_sig_sub = Y_sig_sub,
      expr_sub = expr_sub,
      chunk_size = shape_chunk_size,
      lambda_base = 3,
      lambda_negL = 1,
      lambda_negR = 1,
      lambda_ret = 3,
      lambda_left = 1,
      lambda_peak = 3,
      tau_peak = 0.08,
      early_fixed = 1:2
    )
    
    S_pre_shape[shape_idx] <- shape$S_pre_shape
    pre_peak_pos[shape_idx] <- shape$pre_peak_pos
    pre_baseline[shape_idx] <- shape$pre_baseline
    first_open[shape_idx] <- shape$first_open
    
    RNA_pre_max_frac[shape_idx] <- shape$RNA_pre_max_frac
    RNA_pre_min_frac[shape_idx] <- shape$RNA_pre_min_frac
    ATAC_pre_max_frac[shape_idx] <- shape$ATAC_pre_max_frac
    ATAC_pre_min_frac[shape_idx] <- shape$ATAC_pre_min_frac
    RNA_post_max_frac[shape_idx] <- shape$RNA_post_max_frac
    RNA_post_min_frac[shape_idx] <- shape$RNA_post_min_frac
    ATAC_post_max_frac[shape_idx] <- shape$ATAC_post_max_frac
    ATAC_post_min_frac[shape_idx] <- shape$ATAC_post_min_frac
    resid_mat[shape_idx, ] <- shape$resid_mat
  }
  
  related_idx <- which(
    p_pos >= link_alt_min_prob &
      is.finite(r_vec) &
      is.finite(S_pre_shape)
  )
  
  p_pre_cond <- numeric(n_pairs)
  p_sync_cond <- rep(1, n_pairs)
  
  pos_idx <- related_idx[S_pre_shape[related_idx] > 0]
  if (length(pos_idx) >= 20L) {
    # This is exactly the empirical CDF value, including ties.
    p_pre_cond[pos_idx] <-
      rank(S_pre_shape[pos_idx], ties.method = "max") / length(pos_idx)
    p_sync_cond[related_idx] <- 1 - p_pre_cond[related_idx]
  }
  
  resid_period <- colnames(Y_sig_sub)
  resid_period <- sub("^sec_", "", resid_period)
  resid_period <- sub("_(stroma|epithelial)$", "", resid_period)
  colnames(resid_mat) <- paste0("resid_", resid_period)
  
  pair_df <- data.frame(
    peak = rownames(Y_sig_sub)[peak_idx],
    gene = rownames(expr_sub)[gene_idx],
    S_link = r_vec,
    z_link = link$z_link,
    P_rel_link = p_pos,
    P_unrel_link = p_unrel,
    S_pre_shape = S_pre_shape,
    pre_peak_pos = pre_peak_pos,
    pre_baseline = pre_baseline,
    first_open = first_open,
    RNA_pre_max_frac = RNA_pre_max_frac,
    RNA_pre_min_frac = RNA_pre_min_frac,
    ATAC_pre_max_frac = ATAC_pre_max_frac,
    ATAC_pre_min_frac = ATAC_pre_min_frac,
    RNA_post_max_frac = RNA_post_max_frac,
    RNA_post_min_frac = RNA_post_min_frac,
    ATAC_post_max_frac = ATAC_post_max_frac,
    ATAC_post_min_frac = ATAC_post_min_frac,
    P_pre_given_link = p_pre_cond,
    P_sync_given_link = p_sync_cond,
    P_pre_overall = p_pos * p_pre_cond,
    P_sync_overall = p_pos * p_sync_cond,
    P_unrel_overall = p_unrel,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  pair_df <- cbind(
    pair_df,
    as.data.frame(resid_mat, check.names = FALSE)
  )
  
  message("Done.")
  list(
    pair_stats = pair_df,
    fit_link = fit_link,
    fit_pre = NULL,
    n_pairs = n_pairs,
    n_shape_scored = length(shape_idx),
    mixture_fit_n = length(r_fit)
  )
}








####Script_51_get_pre-open_score_by_statistic_method####

#### Build the binary and distance matrices

## This step is slow. Run the relevant part of make_peak_gene_window_matrix first.
gene_anno_all <- getGeneAnnotation(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)$genes
gene_rep <- get_representative_tss_by_symbol(gene_anno_all, verbose = T)



res_peak_gene_1Mb <- make_peak_gene_window_matrix(
  proj         = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  gene_symbols = ga.sel$symbol,   # Candidate gene list
  half_window  = 1000000,           # 1,000 kb upstream and downstream, for 2 Mb of flanking sequence
  gene_rep = gene_rep
)




PeakGene_TSSdist_1Mb <- make_peak_gene_tss_distance_matrix(
  proj          = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  peak_gene_mat = res_peak_gene_1Mb$mat,
  gene_rep = res_peak_gene_1Mb$gene_rep
)



sec_stroma_cols <- c("sec_PO3_stroma","sec_PO5_stroma","sec_PO7_stroma","sec_PO9_stroma","sec_PO11_stroma","sec_PO13_stroma")
sec_epi_cols <- c("sec_PO3_epithelial","sec_PO5_epithelial","sec_PO7_epithelial","sec_PO9_epithelial","sec_PO11_epithelial","sec_PO13_epithelial")

# Stromal correlations using sec_*_stroma columns
PeakGene_cor_1Mb_stroma <- make_peak_gene_cor_matrix_subset(
  Y_z             = Y_sig,
  avg_list_binary = avg_list_binary,
  peak_gene_mat   = res_peak_gene_1Mb$mat,
  use_cols        = sec_stroma_cols,
  method          = "pearson"   # Set this option to use Spearman correlation
)

# Epithelial correlations using sec_*_epithelial columns
PeakGene_cor_1Mb_epithelial <- make_peak_gene_cor_matrix_subset(
  Y_z             = Y_sig,
  avg_list_binary = avg_list_binary,
  peak_gene_mat   = res_peak_gene_1Mb$mat,
  use_cols        = sec_epi_cols,
  method          = "pearson"
)




# New constraints
residual_modification_f <- make_S_f_2anchors(0.2, 0.1, 0.8, 0.95)
curve(residual_modification_f(x), 0, 1, n = 1000, xlab="x", ylab="f(x)")




#### Apply the function

##stroma
time_order <- sec_stroma_cols
res_bayes_stroma_1Mb <- run_preopen_bayes(
  Y_sig            = Y_sig,
  avg_list_binary  = avg_list_binary,
  res_mat          = res_peak_gene_1Mb$mat,
  time_order       = time_order,
  link_alt_min_prob = 0.5   # Adjust later if needed
)

res_bayes_stroma_1Mb <- run_preopen_bayes_fast(
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  res_mat = res_peak_gene_1Mb$mat,
  time_order = time_order,
  
  link_alt_min_prob = 0.5,
  
  # Chunk size
  pair_chunk_size = 200000L,
  shape_chunk_size = 100000L,
  
  # Compute shape for all candidate pairs because pre-opening can occur with low correlation
  shape_only_related = FALSE,
  
  # Sample up to 100,000 correlations to fit the mixture model
  mix_fit_n = 100000L,
  
  seed = 1L
)


##epithelial
time_order <- sec_epi_cols
res_bayes_epi_1Mb <- run_preopen_bayes(
  Y_sig            = Y_sig,
  avg_list_binary  = avg_list_binary,
  res_mat          = res_peak_gene_1Mb$mat,
  time_order       = time_order,
  link_alt_min_prob = 0.5   # Adjust later if needed
)

res_bayes_epi_1Mb <- run_preopen_bayes_fast(
  Y_sig = Y_sig,
  avg_list_binary = avg_list_binary,
  res_mat = res_peak_gene_1Mb$mat,
  time_order = time_order,
  
  link_alt_min_prob = 0.5,
  
  # Chunk size
  pair_chunk_size = 200000L,
  shape_chunk_size = 100000L,
  
  # Compute shape for all candidate pairs because pre-opening can occur with low correlation
  shape_only_related = FALSE,
  
  # Sample up to 100,000 correlations to fit the mixture model
  mix_fit_n = 100000L,
  
  seed = 1L
)

# 
# res_bayes_stroma_1Mb <- readRDS('Final_sample_and_doc/26_07_13_res_bayes_stroma_1Mb.rds')
# res_bayes_epi_1Mb <- readRDS('Final_sample_and_doc/26_07_13_res_bayes_epi_1Mb.rds')
# 

## Plot triangle scores
dist_mat <- PeakGene_TSSdist_1Mb$dist_mat
stopifnot(!is.null(rownames(PeakGene_TSSdist_1Mb$dist_mat)), !is.null(colnames(PeakGene_TSSdist_1Mb$dist_mat)))

## Stroma
pair_df_stroma <- res_bayes_stroma_1Mb$pair_stats
dim(pair_df_stroma)

stopifnot(all(c("gene","peak") %in% colnames(pair_df_stroma)))


## Add distances
# Convert gene and peak IDs from factors to character strings if needed
pair_df_stroma$peak <- as.character(pair_df_stroma$peak)
pair_df_stroma$gene <- as.character(pair_df_stroma$gene)
# 1) Map peak and gene names to integer row and column indices
ri <- match(pair_df_stroma$peak, rownames(dist_mat))
ci <- match(pair_df_stroma$gene, colnames(dist_mat))
# 2) Look up paired entries, using NA if either ID is unmatched
out <- rep(NA_real_, nrow(pair_df_stroma))
ok  <- !is.na(ri) & !is.na(ci)
out[ok] <- dist_mat[cbind(ri[ok], ci[ok])]
pair_df_stroma$distance_to_TSS <- as.numeric(out)


max(pair_df_stroma$distance_to_TSS)

## Plot S_pre_shape (triangle scores). Exclude organoid genes and leave S_link unrestricted.
aaa <- pair_df_stroma[!is.na(pair_df_stroma$S_link) & pair_df_stroma$S_pre_shape > -0.4 &
                        pair_df_stroma$RNA_pre_max_frac < 0.5 &
                        pair_df_stroma$RNA_pre_max_frac < 0.5 &
                        pair_df_stroma$ATAC_post_min_frac > 0.6, ]
aaa <- pair_df_stroma[pair_df_stroma$gene %notin% gene_only_in_Organoid & pair_df_stroma$S_pre_shape > -0.4, ]




plot_distribution(aaa$S_pre_shape)


tri_Model_result_stroma_1Mb <- pair_df_stroma[!is.na(pair_df_stroma$S_link) & pair_df_stroma$S_pre_shape > -0.125 &
                                                pair_df_stroma$RNA_pre_max_frac < 0.5 &
                                                pair_df_stroma$ATAC_post_min_frac > 0.6 &
                                                pair_df_stroma$RNA_post_min_frac < 0.5,]
tri_Model_result_stroma_1Mb <- tri_Model_result_stroma_1Mb[order(tri_Model_result_stroma_1Mb$S_pre_shape, decreasing = TRUE), ]
dim(tri_Model_result_stroma_1Mb)
# plot_distribution(x = log10(pair_df_stroma$distance_to_TSS), bins = 100)
plot_distribution(x = log10(tri_Model_result_stroma_1Mb$distance_to_TSS))





pdf('ATAC_plot/25_05_21_tri_Model_result_stroma_1Mb.pdf')
for(i in seq(500)){
  tryCatch({
    plot_4_line_using_peak_and_gene(peak_to_use = tri_Model_result_stroma_1Mb$peak[i], gene_to_check = tri_Model_result_stroma_1Mb$gene[i], group_order = c(sec_stroma_cols, sec_epi_cols))
    
  },
  error = function(e){
    message(
      "出错：i = ", i,
      ", peak = ", tri_Model_result_stroma_1Mb$peak[i],
      ", gene = ", tri_Model_result_stroma_1Mb$gene[i],
      "\n错误信息: ", e$message, "\n"
    )
  })
  
}
dev.off()





## Glandular epithelium
pair_df_epi <- res_bayes_epi_1Mb$pair_stats

stopifnot(all(c("gene","peak") %in% colnames(pair_df_epi)))

# Convert gene and peak IDs from factors to character strings if needed
pair_df_epi$peak <- as.character(pair_df_epi$peak)
pair_df_epi$gene <- as.character(pair_df_epi$gene)
# 1) Map peak and gene names to integer row and column indices
ri <- match(pair_df_epi$peak, rownames(dist_mat))
ci <- match(pair_df_epi$gene, colnames(dist_mat))
# 2) Look up paired entries, using NA if either ID is unmatched
out <- rep(NA_real_, nrow(pair_df_epi))
ok  <- !is.na(ri) & !is.na(ci)
out[ok] <- dist_mat[cbind(ri[ok], ci[ok])]
pair_df_epi$distance_to_TSS <- as.numeric(out)

max(pair_df_epi$distance_to_TSS)

# Inspect results
# head(pair_df_epi)
# plot_distribution(pair_df_epi$S_pre_shape[pair_df_epi$S_pre_shape > -0.4], bins = 100)
aaa <- pair_df_epi[pair_df_epi$gene %notin% gene_only_in_Organoid & pair_df_epi$S_pre_shape > -0.5, ]
plot_distribution(x = aaa$S_pre_shape)


tri_Model_result_epi_1Mb <- pair_df_epi[!is.na(pair_df_epi$S_link) & pair_df_epi$S_pre_shape > -0.125 &
                                          pair_df_epi$RNA_pre_max_frac < 0.5 &
                                          pair_df_epi$ATAC_post_min_frac > 0.6 &
                                          pair_df_epi$RNA_post_min_frac < 0.5,]
# tri_Model_result_epi <- pair_df_epi[pair_df_epi$P_unrel_link < 0.5 & !is.na(pair_df_epi$S_pre) & pair_df_epi$S_pre_shape > -0.12, ]

tri_Model_result_epi_1Mb <- tri_Model_result_epi_1Mb[order(tri_Model_result_epi_1Mb$S_pre_shape, decreasing = TRUE), ]
dim(tri_Model_result_epi_1Mb)

tri_Model_result_epi_500kb <- tri_Model_result_epi_1Mb[tri_Model_result_epi_1Mb$distance_to_TSS <= 500000, ]

