# =============================================================================
# make_fixtures.R — Create small synthetic inputs for tests and headless validation
#
# Fixtures mirror the real input structure:
#   avg: dgCMatrix (genes by cell type and stage)
#   Y_sig: dense matrix (peaks by cell type and stage)
#   dist: list(dist_mat = dgCMatrix)
#   res: list with a binary dgCMatrix in mat
#   tri_stroma / tri_epi / sync_stroma / sync_epi: data.frame
# =============================================================================

make_fixtures <- function(out_dir, n_peaks = 30, n_genes = 20,
                          seed = 42L, with_missing_gene = TRUE) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  set.seed(seed)

  stages <- c("PO3", "PO5", "PO7", "PO9", "PO11", "PO13")
  sec_stroma_cols <- paste0("sec_", stages, "_stroma")
  sec_epi_cols    <- paste0("sec_", stages, "_epithelial")
  extra_cols <- c("RIF_PO7_E10_241226_Add_stroma", "RIF_PO7_E10_241226_Add_epithelial")
  all_cols <- c(extra_cols, sec_stroma_cols, sec_epi_cols)

  peaks <- sprintf("chr1:%d-%d", seq(1000, by = 500, length.out = n_peaks),
                   seq(1500, by = 500, length.out = n_peaks))
  genes <- sprintf("GENE%02d", seq_len(n_genes))

  # Y_sig: dense matrix with values in [0, 1]
  Y_sig <- matrix(runif(n_peaks * length(all_cols), 0, 1),
                  nrow = n_peaks, dimnames = list(peaks, all_cols))

  # avg: sparse dgCMatrix with cell fractions in [0, 1]
  avg_m <- matrix(0, nrow = n_genes, ncol = length(all_cols),
                  dimnames = list(genes, all_cols))
  for (g in genes) {
    for (c in all_cols) {
      if (runif(1) < 0.55) avg_m[g, c] <- runif(1)
    }
  }
  avg <- Matrix::Matrix(avg_m, sparse = TRUE)

  # res: sparse binary matrix with about 35% ones, with distances in dist at those entries
  res_m <- matrix(0L, nrow = n_peaks, ncol = n_genes,
                  dimnames = list(peaks, genes))
  dist_m <- matrix(0, nrow = n_peaks, ncol = n_genes,
                   dimnames = list(peaks, genes))
  for (p in peaks) {
    for (g in genes) {
      if (runif(1) < 0.35) {
        res_m[p, g] <- 1L
        dist_m[p, g] <- sample(c(-100000, -50000, 0, 50000, 100000), 1)
      }
    }
  }
  # Include a zero-distance candidate in dist to test the boundary
  if (any(res_m == 1L) && !any(dist_m[res_m == 1L] == 0)) {
    idx <- which(res_m == 1L, arr.ind = TRUE)[1, ]
    dist_m[idx[1], idx[2]] <- 0
  }
  res_mat <- Matrix::Matrix(res_m, sparse = TRUE)
  dist_mat <- Matrix::Matrix(dist_m, sparse = TRUE)

  # Statistics tables
  make_df <- function(pairs, s_col, s_vals) {
    data.frame(
      peak = pairs$peak,
      gene = pairs$gene,
      S = s_vals,
      stringsAsFactors = FALSE
    )
  }
  cand_ones <- which(res_m == 1L, arr.ind = TRUE)
  tri_pairs <- data.frame(
    peak = peaks[cand_ones[, 1]],
    gene = genes[cand_ones[, 2]],
    stringsAsFactors = FALSE
  )
  tri_pairs <- tri_pairs[seq_len(min(12, nrow(tri_pairs))), , drop = FALSE]

  tri_stroma <- data.frame(
    peak = tri_pairs$peak,
    gene = tri_pairs$gene,
    S_pre_shape = runif(nrow(tri_pairs), 0, 0.5),
    stringsAsFactors = FALSE
  )
  tri_epi <- data.frame(
    peak = tri_pairs$peak[seq_len(min(8, nrow(tri_pairs)))],
    gene = tri_pairs$gene[seq_len(min(8, nrow(tri_pairs)))],
    S_pre_shape = runif(min(8, nrow(tri_pairs)), 0, 0.5),
    stringsAsFactors = FALSE
  )
  # Add a candidate with a score exactly at the threshold to test strict exclusion
  if (nrow(tri_stroma) >= 2) {
    tri_stroma$S_pre_shape[2] <- 0.2
  }

  sync_stroma <- data.frame(
    peak = tri_pairs$peak,
    gene = tri_pairs$gene,
    S_final = runif(nrow(tri_pairs), 0, 0.9),
    stringsAsFactors = FALSE
  )
  sync_epi <- data.frame(
    peak = tri_pairs$peak[seq_len(min(8, nrow(tri_pairs)))],
    gene = tri_pairs$gene[seq_len(min(8, nrow(tri_pairs)))],
    S_final = runif(min(8, nrow(tri_pairs)), 0, 0.9),
    stringsAsFactors = FALSE
  )
  if (nrow(sync_stroma) >= 3) {
    sync_stroma$S_final[3] <- 0.7  # A score equal to the threshold must be excluded
  }

  # Optional gene missing from avg to test exclusion
  if (with_missing_gene && nrow(tri_stroma) >= 4) {
    tri_stroma$gene[4] <- "GENE_NOT_IN_AVG"
  }
  # Peak missing from res to test exclusion
  if (with_missing_gene && nrow(sync_epi) >= 4) {
    sync_epi$peak[4] <- "chr99:1-100"
  }

  # Output files
  saveRDS(avg, file.path(out_dir, "avg.rds"))
  saveRDS(Y_sig, file.path(out_dir, "Y_sig.rds"))
  saveRDS(list(dist_mat = dist_mat, peak_center_gr = NULL),
          file.path(out_dir, "dist.rds"))
  saveRDS(list(mat = res_mat, peaks_gr = NULL, genes_gr = NULL, tss_gr = NULL),
          file.path(out_dir, "res.rds"))
  saveRDS(tri_stroma, file.path(out_dir, "tri_stroma.rds"))
  saveRDS(tri_epi, file.path(out_dir, "tri_epi.rds"))
  saveRDS(sync_stroma, file.path(out_dir, "sync_stroma.rds"))
  saveRDS(sync_epi, file.path(out_dir, "sync_epi.rds"))

  cfg <- list(
    paths = list(
      root_dir = out_dir,
      app_dir = ".",
      database = "labels/annotations.sqlite",
      data_cache = "data_cache",
      logs = "logs",
      labels = "labels"
    ),
    data_files = list(
      avg = "avg.rds", Y_sig = "Y_sig.rds", dist = "dist.rds", res = "res.rds",
      tri_stroma = "tri_stroma.rds", tri_epi = "tri_epi.rds",
      sync_stroma = "sync_stroma.rds", sync_epi = "sync_epi.rds"
    ),
    stages = list(
      labels = stages,
      sec_stroma_cols = sec_stroma_cols,
      sec_epi_cols = sec_epi_cols
    ),
    columns = list(
      peak = list(tri_stroma = "peak", tri_epi = "peak",
                  sync_stroma = "peak", sync_epi = "peak"),
      gene = list(tri_stroma = "gene", tri_epi = "gene",
                  sync_stroma = "gene", sync_epi = "gene"),
      S_pre_shape = list(tri_stroma = "S_pre_shape", tri_epi = "S_pre_shape"),
      S_final = list(sync_stroma = "S_final", sync_epi = "S_final")
    ),
    mandatory = list(
      tri_S_pre_shape_threshold = 0.2,
      sync_S_final_threshold = 0.7
    ),
    sampling = list(
      random_seed = 20260806,
      n_background_per_context = 5,
      sync_repress = list(
        sample_size = 100,
        max_per_context = 3,
        dynamic_range_peak = 0.1,
        dynamic_range_gene = 0.1,
        min_drop_peak = 0.1,
        min_drop_gene = 0.1,
        max_window_gap = 1
      )
    ),
    queue_priorities = list(
      mandatory = 1,
      sync_repress_heuristic = 2,
      random_background = 3
    ),
    label_rules = list(
      hard_label_min_top = 7,
      hard_label_min_margin = 2,
      hard_label_min_confidence = 6
    ),
    app = list(
      blind_default = TRUE,
      annotator_default = "tester",
      app_version = "test-1.0.0"
    ),
    qc = list(n_consistency_repeats = 2)
  )
  yaml::write_yaml(cfg, file.path(out_dir, "config.yml"))

  list(
    dir = out_dir,
    n_peaks = n_peaks,
    n_genes = n_genes,
    peaks = peaks,
    genes = genes,
    sec_stroma_cols = sec_stroma_cols,
    sec_epi_cols = sec_epi_cols
  )
}
