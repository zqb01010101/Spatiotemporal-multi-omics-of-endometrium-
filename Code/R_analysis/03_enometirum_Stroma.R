#### Helper functions ####

### Update annotations within a Seurat object using a named mapping list
change_annotation <- function(seurat_obj, change_list, old_cluster, new_cluster = NULL, retain_new_cluster = FALSE) {
  if (is.null(new_cluster)) {
    new_cluster <- old_cluster
  }
  
  meta_data <- seurat_obj@meta.data
  old_vals <- as.character(meta_data[[old_cluster]])
  
  # Reuse existing new_cluster values when available
  # Otherwise initialize new_cluster from old_cluster
  if (new_cluster %in% colnames(meta_data)) {
    new_vals <- as.character(meta_data[[new_cluster]])
  } else {
    new_vals <- old_vals
  }
  
  # Default handling of unmatched values
  if (retain_new_cluster) {
    result <- new_vals
  } else {
    result <- old_vals
  }
  
  # Find old_cluster values listed in change_list
  idx <- old_vals %in% names(change_list)
  
  # Replace matched values
  result[idx] <- unname(unlist(change_list[old_vals[idx]]))
  
  meta_data[[new_cluster]] <- result
  seurat_obj@meta.data <- meta_data
  
  return(seurat_obj)
}


## Plot spatial gene expression with fixed limits
my_spatial_gene_plot <- function(
    seurat_obj,
    compared_gene = NULL,
    sample_name = NULL,
    feature = NULL,
    assay_to_use = "SCT"
) {
  
  # Provide compared_gene, feature, or both
  if (is.null(compared_gene) && is.null(feature)) {
    stop("请提供 compared_gene 或 feature。")
  }
  
  # Get all metadata
  meta_data <- seurat_obj@meta.data
  
  # Check required metadata columns
  required_columns <- c("orig.ident", "nA", "nB")
  missing_columns <- setdiff(required_columns, colnames(meta_data))
  
  if (length(missing_columns) > 0) {
    stop(
      "seurat_obj@meta.data 中缺少以下列：",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Select cell names without creating a Seurat subset
  if (!is.null(sample_name)) {
    cells_use <- rownames(meta_data)[
      meta_data$orig.ident %in% sample_name
    ]
    
    if (length(cells_use) == 0) {
      stop(
        "没有找到 sample_name 对应的细胞：",
        paste(sample_name, collapse = ", ")
      )
    }
  } else {
    cells_use <- rownames(meta_data)
  }
  
  # Extract metadata only for the cells being plotted
  plot_data <- meta_data[cells_use, , drop = FALSE]
  
  # Convert via character to avoid using factor level codes
  plot_data$nA <- as.numeric(as.character(plot_data$nA))
  plot_data$nB <- as.numeric(as.character(plot_data$nB))
  
  # Use the sole assay automatically
  if (length(seurat_obj@assays) == 1) {
    assay_to_use <- names(seurat_obj@assays)[[1]]
  }
  
  # Plot by gene expression
  if (!is.null(compared_gene)) {
    
    # Check that the assay exists
    if (!assay_to_use %in% names(seurat_obj@assays)) {
      stop(
        "找不到 assay：", assay_to_use,
        "。当前可用 assay 为：",
        paste(names(seurat_obj@assays), collapse = ", ")
      )
    }
    
    # Check that the gene exists
    if (!compared_gene %in%
        rownames(seurat_obj@assays[[assay_to_use]]@data)) {
      stop(
        "在 assay '", assay_to_use,
        "' 的 data 中没有找到基因：",
        compared_gene
      )
    }
    
    # Extract the target gene and cells only
    gene_expression <- seurat_obj@assays[[assay_to_use]]@data[
      compared_gene,
      cells_use,
      drop = TRUE
    ]
    
    # Align by cell name
    plot_data$gene <- as.numeric(
      gene_expression[rownames(plot_data)]
    )
    
    value_name <- compared_gene
    
  } else {
    
    # Plot by a metadata feature
    if (!feature %in% colnames(meta_data)) {
      stop(
        "seurat_obj@meta.data 中没有找到 feature：",
        feature
      )
    }
    
    plot_data$gene <- as.numeric(
      as.character(plot_data[[feature]])
    )
    
    value_name <- feature
  }
  
  # Set the plot title
  if (!is.null(sample_name)) {
    sample_label <- paste(sample_name, collapse = ", ")
  } else {
    sample_label <- paste(
      unique(plot_data$orig.ident),
      collapse = ", "
    )
  }
  
  # Set the fill scale so zero always maps to gray
  max_gene <- max(plot_data$gene, na.rm = TRUE)
  
  if (!is.finite(max_gene)) {
    max_gene <- 0
  }
  
  fill_upper <- if (max_gene > 0) max_gene else 1
  
  p <- ggplot(
    plot_data,
    aes(
      x = nB,
      y = nA,
      fill = gene
    )
  ) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradient(
      low = "grey",
      high = "red",
      limits = c(0, fill_upper),
      oob = scales::squish,
      na.value = "grey"
    ) +
    labs(
      title = paste0(sample_label, "_", value_name),
      x = "nB (X Coordinate)",
      y = "nA (Y Coordinate)",
      fill = value_name
    ) +
    coord_fixed()  # + guides(fill = "none")
  
  return(p)
}


## Inspect spatial distributions interactively by selecting a sample using orig.ident
view_my_spatial_plot <- function(
    seurat_object,
    cluster = "exclusive_cell_type",
    colors = my_colors2,
    sample_name = NULL
) {
  
  # Get metadata
  meta_data <- seurat_object[[]]
  
  # Check required metadata columns
  required_columns <- c(
    "orig.ident",
    "nA",
    "nB",
    "nCount_RNA",
    cluster
  )
  
  missing_columns <- setdiff(required_columns, colnames(meta_data))
  
  if (length(missing_columns) > 0) {
    stop(
      "Missing metadata columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Select cells by sample name
  if (is.null(sample_name)) {
    
    seurat_obj <- seurat_object
    
  } else {
    
    orig_ident <- as.character(meta_data$orig.ident)
    
    keep_cells <- rownames(meta_data)[orig_ident == sample_name]
    
    if (length(keep_cells) == 0) {
      stop(
        "No cells found for sample_name: '",
        sample_name,
        "'.\nAvailable samples:\n",
        paste(unique(orig_ident), collapse = "\n")
      )
    }
    
    seurat_obj <- subset(
      seurat_object,
      cells = keep_cells
    )
  }
  
  # Extract the filtered metadata
  plot_data <- seurat_obj[[]]
  
  # Convert coordinates to numeric
  # Convert via character to avoid using factor level codes
  plot_data$nA <- as.numeric(as.character(plot_data$nA))
  plot_data$nB <- as.numeric(as.character(plot_data$nB))
  
  if (anyNA(plot_data$nA) || anyNA(plot_data$nB)) {
    warning("Some nA or nB values could not be converted to numeric.")
  }
  
  # Define plot variables before passing them to Plotly formulas
  plot_data$plot_cluster <- factor(plot_data[[cluster]])
  
  plot_data$hover_text <- paste0(
    cluster, ": ", plot_data[[cluster]],
    "<br>nCount_RNA: ", plot_data$nCount_RNA,
    "<br>orig.ident: ", plot_data$orig.ident,
    "<br>nA: ", plot_data$nA,
    "<br>nB: ", plot_data$nB
    
  )
  
  sample_title <- paste(
    unique(as.character(plot_data$orig.ident)),
    collapse = ", "
  )
  
  title_text <- paste0(sample_title, "   ", cluster)
  
  plotly::plot_ly(
    data = plot_data,
    x = ~nB,
    y = ~nA,
    type = "scatter",
    mode = "markers",
    color = ~plot_cluster,
    colors = colors,
    hoverinfo = "text",
    text = ~hover_text,
    marker = list(
      size = 7,
      symbol = "square"
    )
  ) |>
    plotly::layout(
      title = title_text,
      xaxis = list(title = "nB (X Coordinate)"),
      yaxis = list(title = "nA (Y Coordinate)")
    )
}






draw_KEGG_picture <- function(gene_set, 
                              top_n = 10, 
                              n_to_plot = NULL,
                              sort_by = "pvalue",
                              color_by_p.adjust = TRUE,
                              single_color = "steelblue",
                              y_text_size = 20) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  
  gene.df <- bitr(gene_set,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
  gene <- gene.df$ENTREZID
  
  kk <- enrichKEGG(gene = gene,
                   keyType = "kegg",
                   organism = "human",
                   qvalueCutoff = 0.05,
                   pvalueCutoff = 0.05)
  
  print(kk)
  
  hh <- as.data.frame(kk)
  # Preserve all significant enrichments before filtering by top_n
  hh_all <- hh
  
  ## Map Gene_ID to Gene_Symbol
  all_entrez <- unique(
    unlist(strsplit(hh_all$geneID, "/"))
  )
  
  entrez_to_symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = all_entrez,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  hh_all$geneSymbol <- vapply(
    strsplit(hh_all$geneID, "/"),
    function(ids) {
      symbols <- unname(entrez_to_symbol[ids])
      symbols <- symbols[!is.na(symbols)]
      paste(symbols, collapse = "/")
    },
    character(1)
  )
  
  if (nrow(hh) == 0) {
    message("没有显著富集到 KEGG 通路")
    return(NULL)
  }
  
  if (!sort_by %in% c("pvalue", "p.adjust")) {
    stop("sort_by 只能是 'pvalue' 或 'p.adjust'")
  }
  
  # Sort by significance first
  hh <- hh[order(hh[[sort_by]], decreasing = FALSE), ]
  rownames(hh) <- seq_len(nrow(hh))
  
  # Record significance ranks
  hh$rank <- seq_len(nrow(hh))
  
  # If n_to_plot is set, plot pathways at those ranks
  # Otherwise plot the top_n pathways
  if (!is.null(n_to_plot)) {
    if (!is.numeric(n_to_plot)) {
      stop("n_to_plot 必须是数值向量，例如 c(1, 4, 6)")
    }
    
    if (any(n_to_plot != as.integer(n_to_plot)) || any(n_to_plot <= 0)) {
      stop("n_to_plot 必须是正整数向量，例如 c(1, 4, 6)")
    }
    
    if (any(n_to_plot > nrow(hh))) {
      stop("n_to_plot 中存在超过富集通路数量的编号")
    }
    
    hh <- hh[n_to_plot, , drop = FALSE]
  } else {
    hh <- head(hh, top_n)
  }
  
  # Show the most significant pathways at the top
  hh$order <- factor(hh$Description, levels = rev(hh$Description))
  
  if (color_by_p.adjust) {
    p <- ggplot(hh, aes(y = order, x = Count, fill = p.adjust)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_gradient(low = "red", high = "blue") +
      labs(fill = "p.adjust")
  } else {
    p <- ggplot(hh, aes(y = order, x = Count)) +
      geom_bar(stat = "identity", width = 0.7, fill = single_color)
  }
  
  p <- p +
    labs(title = "KEGG Pathways Enrichment",
         x = "Gene numbers",
         y = "Pathways") +
    theme_bw() +
    theme(
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      axis.text.y = element_text(size = y_text_size, face = "bold")
    )
  
  
  # ==============================
  # Second plot with -log10(P value) on the x-axis
  # ==============================
  
  # Avoid Inf when very small P values underflow to zero
  hh$minus_log10_pvalue <- -log10(
    pmax(hh$pvalue, .Machine$double.xmin)
  )
  
  if (color_by_p.adjust) {
    p_pvalue <- ggplot(
      hh,
      aes(
        y = order,
        x = minus_log10_pvalue,
        fill = p.adjust
      )
    ) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_gradient(low = "red", high = "blue") +
      labs(fill = "p.adjust")
  } else {
    p_pvalue <- ggplot(
      hh,
      aes(
        y = order,
        x = minus_log10_pvalue
      )
    ) +
      geom_bar(
        stat = "identity",
        width = 0.7,
        fill = single_color
      )
  }
  
  p_pvalue <- p_pvalue +
    labs(
      title = "KEGG Pathways Enrichment",
      x = expression(-log[10](italic(P)~value)),
      y = "Pathways"
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_text(face = "bold", size = 20),
      axis.title.y = element_text(face = "bold", size = 20),
      legend.title = element_text(face = "bold", size = 20),
      axis.text.y = element_text(size = y_text_size, face = "bold")
    )
  return(list(
    plot = p,
    pvalue_plot = p_pvalue,
    enrichment = kk,
    result = hh_all,
    plotted_result = hh,
    gene_mapping = gene.df
  ))
}










####Script_31_DEG_of_diff_stage####

pro_and_sec_seurat <- change_annotation(pro_and_sec_seurat, change_list = list('inner_stroma' = 'stroma', 'surface_stroma' = 'stroma', 'immune_stroma' = 'stroma'), old_cluster = 'exclusive_cell_type_4', new_cluster = 'exclusive_cell_type_5')
pro_and_sec_stroma_and_gland <- subset(pro_and_sec_seurat, subset = exclusive_cell_type_5 %in% c('stroma', 'gland'))
pro_and_sec_stroma_and_gland <- subset(pro_and_sec_stroma_and_gland, subset = stage_2 != 'RIF_PO7')


DefaultAssay(pro_and_sec_stroma_and_gland) <- "RNA"
pro_and_sec_stroma_and_gland[['SCT']] <- NULL
# Join layers first so all layers are used
pro_and_sec_stroma_and_gland <- JoinLayers(pro_and_sec_stroma_and_gland)
# Extract counts from the counts layer
counts_data <- LayerData(pro_and_sec_stroma_and_gland, assay = "RNA", layer = "counts")
# Create a new counts-only assay
new_assay <- CreateAssay5Object(counts = counts_data)
# Replace the original assay
pro_and_sec_stroma_and_gland[["RNA"]] <- new_assay
# Verify that only the counts layer remains
Layers(pro_and_sec_stroma_and_gland, assay = "RNA")
# Split the layers again
pro_and_sec_stroma_and_gland[['RNA']] <- split(pro_and_sec_stroma_and_gland[['RNA']], f = pro_and_sec_stroma_and_gland$orig.ident)





pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_4 <- paste0(pro_and_sec_stroma_and_gland$stage_period, '_', pro_and_sec_stroma_and_gland$exclusive_cell_type_4)
pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_5 <- paste0(pro_and_sec_stroma_and_gland$stage_period, '_', pro_and_sec_stroma_and_gland$exclusive_cell_type_5)
pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_5 <- factor(pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_5, levels = c('sec_late_stroma', 'sec_mid_stroma', 'sec_early_stroma', 'pro_stage_stroma', 'sec_late_gland', 'sec_mid_gland', 'sec_early_gland', 'pro_stage_gland'))

Idents(pro_and_sec_stroma_and_gland) <- 'stage_period_exclusive_cell_type_5'
pro_and_sec_stroma_and_gland <- SCTransform(pro_and_sec_stroma_and_gland)
pro_and_sec_stroma_and_gland <- PrepSCTFindMarkers(pro_and_sec_stroma_and_gland)


stroma_and_gland_whole_stage_marker <- FindAllMarkers(pro_and_sec_stroma_and_gland, only.pos = TRUE, assay = 'SCT', min.pct = 0.1)

pro_stroma_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$cluster == 'pro_stage_stroma' & stroma_and_gland_whole_stage_marker$avg_log2FC > 1, ]
sec_early_stroma_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$cluster == 'sec_early_stroma' & stroma_and_gland_whole_stage_marker$avg_log2FC > 1, ]
sec_mid_stroma_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$cluster == 'sec_mid_stroma' & stroma_and_gland_whole_stage_marker$avg_log2FC > 1, ]
sec_late_stroma_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$cluster == 'sec_late_stroma' & stroma_and_gland_whole_stage_marker$avg_log2FC > 1, ]

pro_gland_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$pct.1 > 0.3 & stroma_and_gland_whole_stage_marker$pct.2 < 0.3 & stroma_and_gland_whole_stage_marker$cluster == 'pro_stage_gland' & stroma_and_gland_whole_stage_marker$avg_log2FC > 0.7, ]
sec_early_gland_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$pct.1 > 0.3 & stroma_and_gland_whole_stage_marker$pct.2 < 0.3 & stroma_and_gland_whole_stage_marker$cluster == 'sec_early_gland' & stroma_and_gland_whole_stage_marker$avg_log2FC > 0.7, ]
sec_mid_gland_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$pct.1 > 0.3 & stroma_and_gland_whole_stage_marker$pct.2 < 0.3 & stroma_and_gland_whole_stage_marker$cluster == 'sec_mid_gland' & stroma_and_gland_whole_stage_marker$avg_log2FC > 0.7, ]
sec_late_gland_marker <- stroma_and_gland_whole_stage_marker[stroma_and_gland_whole_stage_marker$pct.1 > 0.3 & stroma_and_gland_whole_stage_marker$pct.2 < 0.3 & stroma_and_gland_whole_stage_marker$cluster == 'sec_late_gland' & stroma_and_gland_whole_stage_marker$avg_log2FC > 0.7, ]



stroma_marker_gene_to_plot <- c(#'VIM', 'COL1A1', 
  'SFRP4', 'MTUS2', 'KCNIP4', 'PHF21B', 'SDK2', 
  'PLCXD3', 'MATN2', 'BRINP1', 'ZBTB16', 'GDF7',
  'ALDH1A3', 'C1R', 'HAND2', 'KCND3', 'SORBS1', 
  'PLAT', 'IGFBP1', 'TIMP3', 'PLAGL1', 'FERMT2')

my_Dotplot_plot(seurat_obj = pro_and_sec_stroma_and_gland, marker_gene = stroma_marker_gene_to_plot, cluster = 'stage_period_exclusive_cell_type_5',
                idents = c('pro_stage_stroma', 'sec_early_stroma', 'sec_mid_stroma', 'sec_late_stroma'))


pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_5 <- factor(pro_and_sec_stroma_and_gland$stage_period_exclusive_cell_type_5, 
                                                                          levels = c('pro_stage_stroma', 'sec_early_stroma', 'sec_mid_stroma', 'sec_late_stroma',
                                                                                     'pro_stage_gland', 'sec_early_gland', 'sec_mid_gland', 'sec_late_gland'))

all_genes <- getFeatures(
  ArchRProj = pro_and_sec_ArchR,
  useMatrix = "GeneScoreMatrix"
)

se_gene <- getMatrixFromProject(
  ArchRProj = pro_and_sec_ArchR,
  useMatrix = "GeneScoreMatrix"
)




gene_list <- stroma_marker_gene_to_plot

gene_map <- setNames(all_genes, toupper(all_genes))

genes_found <- unname(gene_map[toupper(gene_list)])
genes_found <- genes_found[!is.na(genes_found)]
genes_found <- as.character(genes_found)

genes_missing <- gene_list[!toupper(gene_list) %in% toupper(all_genes)]

genes_found
genes_missing


mat_gene_score <- assay(se_gene)

class(mat_gene_score)
dim(mat_gene_score)

gene_names <- rowData(se_gene)$name


rd <- rowData(se_gene)

if ("name" %in% colnames(rd)) {
  gene_names <- rd$name
} else if ("symbol" %in% colnames(rd)) {
  gene_names <- rd$symbol
} else if ("gene_name" %in% colnames(rd)) {
  gene_names <- rd$gene_name
} else {
  gene_names <- getFeatures(
    ArchRProj = pro_and_sec_ArchR,
    useMatrix = "GeneScoreMatrix"
  )
}

gene_names <- as.character(gene_names)



idx <- match(toupper(gene_list), toupper(gene_names))

genes_found <- gene_names[idx[!is.na(idx)]]
genes_missing <- gene_list[is.na(idx)]

genes_found
genes_missing


idx_use <- idx[!is.na(idx)]

mat_use <- mat_gene_score[idx_use, , drop = FALSE]
rownames(mat_use) <- genes_found

mat_gene_score[idx_use, , drop = FALSE]

pro_and_sec_ArchR <- change_Archr_annotation_using_list(archr_project = pro_and_sec_ArchR, change_list = list('sec_PO3' = 'sec_early', 'pro_late' = 'pro_stage', 'sec_PO11' = 'sec_late', 'sec_PO9' = 'sec_mid'), old_cluster = 'stage', new_cluster = 'stage_period', retain_new_cluster = TRUE)
pro_and_sec_ArchR$stage_period_super_cell_type <- paste0(pro_and_sec_ArchR$stage_period, '_', pro_and_sec_ArchR$super_cell_type)

pro_and_sec_ArchR <- change_Archr_annotation_using_list(archr_project = pro_and_sec_ArchR, change_list = list('inner_stroma' = 'stroma', 'surface_stroma' = 'stroma', 'immune_stroma' = 'stroma'), 
                                                        old_cluster = 'exclusive_cell_type_4', new_cluster = 'exclusive_cell_type_5')


group_col <- "stage_period_super_cell_type"
celltype_order <- c('sec_late_stroma', 'sec_mid_stroma', 'sec_early_stroma', 'pro_stage_stroma')


meta <- as.data.frame(
  getCellColData(pro_and_sec_ArchR, select = group_col)
)

common_cells <- intersect(colnames(mat_use), rownames(meta))

mat_use2 <- mat_use[, common_cells, drop = FALSE]
meta2 <- meta[common_cells, , drop = FALSE]

group_vec <- meta2[[group_col]]

avg_mat <- sapply(
  split(seq_along(group_vec), group_vec),
  function(idx) {
    Matrix::rowMeans(mat_use2[, idx, drop = FALSE])
  }
)

avg_mat <- as.matrix(avg_mat)
avg_mat <- avg_mat[, colnames(avg_mat) %in% celltype_order]

avg_mat_z <- t(scale(t(avg_mat)))
# avg_mat_z <- avg_mat
avg_mat_z[is.na(avg_mat_z)] <- 0

gene_order <- gene_list


plot_df_ATAC <- avg_mat_z %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "celltype",
    values_to = "value"
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    celltype = factor(celltype, levels = celltype_order)
  )



## RNA comparison shared by the plots below

avg_rna <- AverageExpression(
  pro_and_sec_seurat,
  group.by = "stage_super_cell_type",
  assays = "SCT",
  slot = "data"
)
avg_rna$SCT
colnames(avg_rna$SCT) <- gsub("-", "_", colnames(avg_rna$SCT))
avg_rna$SCT <- avg_rna$SCT[, celltype_order]

avg_rna_z <- t(scale(t(avg_rna$SCT)))

gene_order <- gene_list
avg_rna_z <- avg_rna_z[gene_order, ]



plot_df_RNA <- avg_rna_z %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "celltype",
    values_to = "value"
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    celltype = factor(celltype, levels = celltype_order)
  )


## Shared plotting code
plot_df <- plot_df_ATAC
plot_df <- plot_df_RNA



celltype_order <- c(
  "pro_stage_stroma",
  "sec_early_stroma",
  "sec_mid_stroma",
  "sec_late_stroma"
)

# Get the current gene order
gene_order <- if (is.factor(plot_df$gene)) {
  levels(plot_df$gene)
} else {
  unique(as.character(plot_df$gene))
}

plot_df$celltype <- factor(
  plot_df$celltype,
  levels = celltype_order
)

# Keep genes in the original data order
plot_df$gene <- factor(
  plot_df$gene,
  levels = unique(plot_df$gene)
)




max_abs <- max(abs(plot_df$value), na.rm = TRUE)

p <- ggplot(
  plot_df,
  aes(x = celltype, y = gene, fill = value)
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "#3F26FF",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_abs, max_abs),
    breaks = pretty(c(-max_abs, max_abs), n = 7),
    oob = scales::squish
  ) +
  scale_y_discrete(
    limits = rev(gene_order)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black",
      face = "plain", 
      size = 10
    ),
    axis.text.y = element_text(angle = 0,
                               color = "black",
                               face = "plain",
                               size = 10
    ),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = "Cell type",
    y = "Gene",
    fill = "Z-score",
    title = 'RNA Gene Expression'
  )

p



ggsave(
  filename = "Final_Plot/Fig_2/Fig_2A_Hotplot_of_stage_marker_gene_of_RNA.pdf",
  plot = p,
  width = 3.5,
  height = 7.5,
  units = "in"
)



ggsave(
  filename = "Final_Plot/Fig_2/Fig_2A_Hotplot_of_stage_marker_gene_of_ATAC.pdf",
  plot = p,
  width = 3.5,
  height = 7.5,
  units = "in"
)





####Script_32_DEG_between_ISCs_and_SSCs_across_different_stage####
unique(pro_and_sec_stroma_and_gland$stage_period_2_exclusive_cell_type_surface_and_inner_stroma)
Idents(pro_and_sec_stroma_and_gland) <- 'stage_period_2_exclusive_cell_type_surface_and_inner_stroma'
Idents(pro_and_sec_stroma_and_gland) <- 'stage_period_exclusive_cell_type_4'

Idents(pro_and_sec_stroma_and_gland) <- 'exclusive_cell_type_surface_and_inner_stroma'
Marker_gene_between_surface_and_inner_stroma <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'surface_stroma', ident.2 = 'inner_stroma', group.by = 'exclusive_cell_type_surface_and_inner_stroma', min.pct = 0.1)
View(Marker_gene_between_surface_and_inner_stroma)
Marker_gene_of_surface_stroma <- Marker_gene_between_surface_and_inner_stroma[Marker_gene_between_surface_and_inner_stroma$avg_log2FC > 0.7 & Marker_gene_between_surface_and_inner_stroma$p_val_adj < 0.01, ]
Marker_gene_of_inner_stroma <- Marker_gene_between_surface_and_inner_stroma[Marker_gene_between_surface_and_inner_stroma$avg_log2FC < -0.7 & Marker_gene_between_surface_and_inner_stroma$p_val_adj < 0.01, ]


Idents(pro_and_sec_stroma_and_gland) <- 'super_cell_type'
stroma_and_gland_marker <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'stroma', ident.2 = 'epithelial', group.by = 'super_cell_type', min.pct = 0.1)
stroma_marker <- rownames(stroma_and_gland_marker)[stroma_and_gland_marker$avg_log2FC > 0.7 & stroma_and_gland_marker$p_val_adj < 0.01]
gland_marker <- rownames(stroma_and_gland_marker)[stroma_and_gland_marker$avg_log2FC < -0.7 & stroma_and_gland_marker$p_val_adj < 0.01]

Idents(pro_and_sec_stroma_and_gland) <- 'stage_period_exclusive_cell_type_4'
Marker_gene_between_surface_and_inner_stroma_of_pro_stage <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'pro_stage_surface_stroma', ident.2 = 'pro_stage_inner_stroma', assay = 'SCT', min.pct = 0.1, logfc.threshold = 0.7)
Marker_gene_between_surface_and_inner_stroma_of_sec_early <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'sec_early_surface_stroma', ident.2 = 'sec_early_inner_stroma', assay = 'SCT', min.pct = 0.1, logfc.threshold = 0.7)
Marker_gene_between_surface_and_inner_stroma_of_sec_mid <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'sec_mid_surface_stroma', ident.2 = 'sec_mid_inner_stroma', assay = 'SCT', min.pct = 0.1, logfc.threshold = 0.7)
Marker_gene_between_surface_and_inner_stroma_of_sec_late <- FindMarkers(object = pro_and_sec_stroma_and_gland, ident.1 = 'sec_late_surface_stroma', ident.2 = 'sec_late_inner_stroma', assay = 'SCT', min.pct = 0.1, logfc.threshold = 0.7)
Marker_gene_between_surface_and_inner_stroma_of_pro_stage <- Marker_gene_between_surface_and_inner_stroma_of_pro_stage[order(Marker_gene_between_surface_and_inner_stroma_of_pro_stage$avg_log2FC, decreasing = TRUE), ]
Marker_gene_between_surface_and_inner_stroma_of_sec_early <- Marker_gene_between_surface_and_inner_stroma_of_sec_early[order(Marker_gene_between_surface_and_inner_stroma_of_sec_early$avg_log2FC, decreasing = TRUE), ]
Marker_gene_between_surface_and_inner_stroma_of_sec_mid <- Marker_gene_between_surface_and_inner_stroma_of_sec_mid[order(Marker_gene_between_surface_and_inner_stroma_of_sec_mid$avg_log2FC, decreasing = TRUE), ]
Marker_gene_between_surface_and_inner_stroma_of_sec_late <- Marker_gene_between_surface_and_inner_stroma_of_sec_late[order(Marker_gene_between_surface_and_inner_stroma_of_sec_late$avg_log2FC, decreasing = TRUE), ]


View(Marker_gene_between_surface_and_inner_stroma_of_pro_stage)
View(Marker_gene_between_surface_and_inner_stroma_of_sec_early)
View(Marker_gene_between_surface_and_inner_stroma_of_sec_mid)
View(Marker_gene_between_surface_and_inner_stroma_of_sec_late)


more_expr_in_pro_stage_surface_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_pro_stage[Marker_gene_between_surface_and_inner_stroma_of_pro_stage$avg_log2FC > 0.5 & Marker_gene_between_surface_and_inner_stroma_of_pro_stage$p_val_adj < 0.05, ]))
more_expr_in_sec_early_surface_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_early[Marker_gene_between_surface_and_inner_stroma_of_sec_early$avg_log2FC > 0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_early$p_val_adj < 0.05, ]))
more_expr_in_sec_mid_surface_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_mid[Marker_gene_between_surface_and_inner_stroma_of_sec_mid$avg_log2FC > 0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_mid$p_val_adj < 0.05, ]))
more_expr_in_sec_late_surface_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_late[Marker_gene_between_surface_and_inner_stroma_of_sec_late$avg_log2FC > 0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_late$p_val_adj < 0.05, ]))

more_expr_in_pro_stage_surface_stroma <- setdiff(more_expr_in_pro_stage_surface_stroma, gland_marker)
more_expr_in_sec_early_surface_stroma <- setdiff(more_expr_in_sec_early_surface_stroma, gland_marker)
more_expr_in_sec_mid_surface_stroma <- setdiff(more_expr_in_sec_mid_surface_stroma, gland_marker)
more_expr_in_sec_late_surface_stroma <- setdiff(more_expr_in_sec_late_surface_stroma, gland_marker)


more_expr_in_pro_stage_inner_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_pro_stage[Marker_gene_between_surface_and_inner_stroma_of_pro_stage$avg_log2FC < -0.5 & Marker_gene_between_surface_and_inner_stroma_of_pro_stage$p_val_adj < 0.05, ]))
more_expr_in_sec_early_inner_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_early[Marker_gene_between_surface_and_inner_stroma_of_sec_early$avg_log2FC < -0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_early$p_val_adj < 0.05, ]))
more_expr_in_sec_mid_inner_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_mid[Marker_gene_between_surface_and_inner_stroma_of_sec_mid$avg_log2FC < -0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_mid$p_val_adj < 0.05, ]))
more_expr_in_sec_late_inner_stroma <- unique(rownames(Marker_gene_between_surface_and_inner_stroma_of_sec_late[Marker_gene_between_surface_and_inner_stroma_of_sec_late$avg_log2FC < -0.5 & Marker_gene_between_surface_and_inner_stroma_of_sec_late$p_val_adj < 0.05, ]))

more_expr_in_pro_stage_inner_stroma <- setdiff(more_expr_in_pro_stage_inner_stroma, gland_marker)
more_expr_in_sec_early_inner_stroma <- setdiff(more_expr_in_sec_early_inner_stroma, gland_marker)
more_expr_in_sec_mid_inner_stroma <- setdiff(more_expr_in_sec_mid_inner_stroma, gland_marker)
more_expr_in_sec_late_inner_stroma <- setdiff(more_expr_in_sec_late_inner_stroma, gland_marker)




gene_to_plot <- 'ADAMTS16'
p_ADAMTS16_PO3 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO3_E65_240205_RNA')
p_ADAMTS16_PO7 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO7_E16_250419_Add')
p_ADAMTS16_PO11 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO11_E15_240413_RNA')


gene_to_plot <- 'SLC35F3'
p_SLC35F3_PO3 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO3_E65_240205_RNA')
p_SLC35F3_PO7 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO7_E16_250419_Add')
p_SLC35F3_PO11 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO11_E15_240413_RNA')



gene_to_plot <- 'LRRC15'
p_LRRC15_PO3 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO3_E65_240205_RNA')
p_LRRC15_PO7 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO7_E16_250419_Add')
p_LRRC15_PO11 <- my_spatial_gene_plot(seurat_obj = pro_and_sec_seurat, compared_gene = gene_to_plot, sample_name = 'sec_PO11_E15_240413_RNA')

sec_PO3_E65_spatial_plot <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'exclusive_cell_type_4', Colors = my_colors_stroma, sample_name = 'sec_PO3_E65_240205_RNA')
sec_PO7_E16_spatial_plot <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'exclusive_cell_type_4', Colors = my_colors_stroma, sample_name = 'sec_PO7_E16_250419_Add')
sec_PO11_E15_spatial_plot <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'exclusive_cell_type_4', Colors = my_colors_stroma, sample_name = 'sec_PO11_E15_240413_RNA')



plots <- list(
  sec_PO3_E65_spatial_plot, p_ADAMTS16_PO3, p_SLC35F3_PO3, p_LRRC15_PO3,
  sec_PO7_E16_spatial_plot, p_ADAMTS16_PO7, p_SLC35F3_PO7, p_LRRC15_PO7,
  sec_PO11_E15_spatial_plot, p_ADAMTS16_PO11, p_SLC35F3_PO11, p_LRRC15_PO11
)


combined_plot <- wrap_plots(
  plots,
  ncol = 4,
  nrow = 3
)

combined_plot


ggsave(
  filename = "Final_Plot/Fig_2/Fig_2B_marker_gene_of_surface_stroma_as_different_stage.pdf",
  plot = combined_plot,
  width = 16,
  height = 16,
  units = "in",
  device = cairo_pdf
)





####Script_33_Functional_Enrichment_Analysis_of_ISCs-Specific_DEGs####

p <- draw_KEGG_picture(more_expr_in_sec_mid_surface_stroma, top_n = 20, n_to_plot = c(1,2,3,9,10), color_by_p.adjust = FALSE, single_color = 'red', y_text_size = 13)
draw_KEGG_picture(more_expr_in_sec_late_surface_stroma, top_n = 5, color_by_p.adjust = FALSE, single_color = 'green', y_text_size = 13)


ggsave("Final_Plot/Fig_2/Fig_2C_KEGG_enrichment_of_marker_gene_of_sec_mid_surface_stroma.pdf",  p$pvalue_plot,  width = 8,  height = 3)






####Script_34_Differences_in_cell–cell_communication_with_TE_between_ISCs_and_SSCs_at_PO7####
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)

options(future.globals.maxSize = 20 * 1024^3)
wd='/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/'
setwd(wd)

pro_and_sec_and_Blastaea_seurat_for_cellchat <- readRDS('integration_sample/26_05_12_pro_and_sec_and_Blastaea_seurat_for_cellchat.rds')

pro_and_sec_and_Blastaea_seurat_for_cellchat$period_exclusive_cell_type_4 <- paste0(pro_and_sec_and_Blastaea_seurat_for_cellchat$stage, '_', pro_and_sec_and_Blastaea_seurat_for_cellchat$exclusive_cell_type_4)
DefaultAssay(pro_and_sec_and_Blastaea_seurat_for_cellchat) <- 'RNA'
pro_and_sec_and_Blastaea_seurat_for_cellchat <- NormalizeData(  pro_and_sec_and_Blastaea_seurat_for_cellchat, assay = "RNA")

Idents(pro_and_sec_and_Blastaea_seurat_for_cellchat) <- 'period_exclusive_cell_type_4'


Blastaea_cellchat_sample <- createCellChat(object = pro_and_sec_and_Blastaea_seurat_for_cellchat,  group.by = "period_exclusive_cell_type_4", assay = "RNA")


CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB)
Blastaea_cellchat_sample@DB <- CellChatDB.use

Blastaea_cellchat_sample <- subsetData(Blastaea_cellchat_sample)

Blastaea_cellchat_sample <- identifyOverExpressedGenes(Blastaea_cellchat_sample)
Blastaea_cellchat_sample <- identifyOverExpressedInteractions(Blastaea_cellchat_sample) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
Blastaea_cellchat_sample <- smoothData(Blastaea_cellchat_sample, adj = PPI.human) # projectData projects ligand-receptor expression in @data.signaling onto the PPI network and saves adjusted values in @data.project


Blastaea_cellchat_sample <- computeCommunProb(Blastaea_cellchat_sample, raw.use = FALSE, population.size = TRUE) # Set raw.use = TRUE to use unadjusted expression
Blastaea_cellchat_sample <- computeCommunProbPathway(Blastaea_cellchat_sample)
Blastaea_cellchat_sample <- filterCommunication(Blastaea_cellchat_sample, min.cells = 3)

saveRDS(Blastaea_cellchat_sample, 'R_cellchat/cellchat_sample/26_05_12_pro_and_sec_and_Blastaea_seurat_for_cellchat.rds')

Blastaea_cellchat_sample <- readRDS('R_cellchat/cellchat_sample/26_05_12_pro_and_sec_and_Blastaea_seurat_for_cellchat.rds')

source('R_cellchat/my_netVisual_bubble.R')


pdf(paste0('R_cellchat/cellchat_plot/', Sys.Date(), '_PO7_stroma_to_Blastaea_communication.pdf'), height = 55)
netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('TE_TE'), sources.use = c('sec_PO7_surface_stroma', 'sec_PO7_inner_stroma'))  # Hierarchy plot
dev.off()

pdf(paste0('R_cellchat/cellchat_plot/', Sys.Date(), '_Blastaea_to_PO7_stroma_to_communication.pdf'), height = 55)
netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, sources.use = c('TE_TE'), targets.use = c('sec_PO7_surface_stroma', 'sec_PO7_inner_stroma'))  # Hierarchy plot
dev.off()


ligand_list <- c('THBS2', 'THBS2', 'THBS2', 'THBS2', 'THBS2', 'SEMA6D', 'SEMA6D', 'ANGPT1')
receptor_list <- c("SDC4", "SDC1", "CD47", "CD36", "ITGA3_ITGB1", "PLXNA1", "PLXNA1_KDR", "TEK")


pairLR.use_for_surface_stroma <- data.frame(
  interaction_name = paste0(ligand_list, '_', receptor_list)
)

p1 <- netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('TE_TE'), sources.use = c('sec_PO7_surface_stroma', 'sec_PO7_inner_stroma'), 
                       pairLR.use = pairLR.use_for_surface_stroma, sort.by.source = TRUE)  # Hierarchy plot

ggsave("Final_Plot/Fig_2/Fig_2C_Cell_Communication_from_PO7_Stroma_to_TE.pdf",  p1,  width = 3.6,  height = 4.1)



pairLR.use_for_surface_stroma <- data.frame(
  interaction_name = pathway_vec <- c("THBS2_ITGA3_ITGB1", "THBS2_CD36", "THBS2_CD47", "THBS2_SDC1", "THBS2_SDC4", 
                                      "TNC_ITGA8_ITGB1", "TNC_SDC1", "TNC_SDC4")
)
p3 <- netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('TE_TE'), sources.use = c('sec_PO3_surface_stroma', 'sec_PO5_surface_stroma', 'sec_PO7_surface_stroma','sec_PO9_surface_stroma', 'sec_PO11_surface_stroma', 'sec_PO13_surface_stroma'), 
                       pairLR.use = pairLR.use_for_surface_stroma, sort.by.source = TRUE)  # Hierarchy plot


ggsave("Final_Plot/Fig_2/Fig_2D_Cell_Communication_from_different_stage_stroma_to_TE.pdf", p3, width = 5, height = 4.4)






####Script_35_KEGG_enrichment_of_immune_stroma_marker####
sec_late_stroma <- subset(pro_and_sec_seurat, subset = stage_period == 'sec_late')
sec_late_stroma <- subset(sec_late_stroma, subset = super_cell_type == 'stroma')
sec_late_stroma[['RNA']] <- split(sec_late_stroma[['RNA']], f = sec_late_stroma$orig.ident)

sec_late_stroma <- SCTransform(sec_late_stroma)

sec_late_stroma <- PrepSCTFindMarkers(sec_late_stroma)
Idents(sec_late_stroma) <- 'exclusive_cell_type_4'
DEG_of_sec_late_stroma <- FindAllMarkers(sec_late_stroma, min.pct = 0.1, only.pos = TRUE)

View(DEG_of_sec_late_stroma)

marker_gene_of_immune_stroma <- FindMarkers(object = sec_late_stroma, ident.1 = 'immune_stroma', ident.2 = 'inner_stroma', group.by = 'exclusive_cell_type_4', min.pct = 0.1)
marker_gene_of_immune_stroma <- marker_gene_of_immune_stroma[marker_gene_of_immune_stroma$avg_log2FC > 0.7 & marker_gene_of_immune_stroma$p_val_adj < 0.01,]
marker_gene_of_immune_stroma <- rownames(marker_gene_of_immune_stroma)

aaa <- draw_KEGG_picture((marker_gene_of_immune_stroma), top_n = 15, n_to_plot = c(2, 3, 7, 8, 9), color_by_p.adjust = FALSE, single_color = 'red', y_text_size = 10)
ggsave(
  filename = "Final_Plot/Fig_2/Fig_2F_KEGG_enrichment_of_immune_stroma_marker_gene.pdf",
  plot = aaa$pvalue_plot,
  width = 7,
  height = 3,
  units = "in"
)


####Script_35_




