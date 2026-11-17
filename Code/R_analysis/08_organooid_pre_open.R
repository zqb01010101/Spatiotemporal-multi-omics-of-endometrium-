
#### Helper functions ####

## Expand and combine tables
# --- 1) Convert inputs to dgCMatrix format with genes as rows and cell types as columns ---
.to_sparse_gene_cell <- function(obj, gene_col = NULL) {
  # For data frames or tibbles, genes can be in a column or in row names
  if (is.data.frame(obj)) {
    df <- obj
    if (!is.null(gene_col) && gene_col %in% colnames(df)) {
      genes <- as.character(df[[gene_col]])
      df[[gene_col]] <- NULL
      rownames(df) <- genes
    } else if ("gene" %in% colnames(df)) {
      genes <- as.character(df[["gene"]])
      df[["gene"]] <- NULL
      rownames(df) <- genes
    } else {
      if (is.null(rownames(df)) || any(rownames(df) == "")) {
        stop("data.frame 里没有 gene 列，也没有有效的 rownames(gene)。请把 gene 放到一列 gene 或 rownames。")
      }
    }
    mat <- as.matrix(df)
    storage.mode(mat) <- "numeric"
    return(Matrix(mat, sparse = TRUE))
  }
  
  # matrix
  if (is.matrix(obj)) {
    mat <- obj
    storage.mode(mat) <- "numeric"
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
      stop("matrix 必须同时具有 rownames(gene) 和 colnames(celltype)。")
    }
    return(Matrix(mat, sparse = TRUE))
  }
  
  # Matrix / dgCMatrix
  if (inherits(obj, "Matrix")) {
    if (is.null(rownames(obj)) || is.null(colnames(obj))) {
      stop("稀疏矩阵(Matrix)必须同时具有 rownames(gene) 和 colnames(celltype)。")
    }
    return(as(obj, "dgCMatrix"))
  }
  
  stop("不支持的输入类型：请提供 data.frame / matrix / Matrix(dgCMatrix)。")
}

# --- 2) Expand to genes_all by cells_all and fill missing entries with zeros ---
.expand_sparse <- function(mat, genes_all, cells_all) {
  mat <- as(mat, "dgCMatrix")
  sm <- summary(mat)  # One-based i, j, x triplets for nonzero entries
  if (nrow(sm) == 0) {
    return(Matrix(0, nrow = length(genes_all), ncol = length(cells_all),
                  sparse = TRUE, dimnames = list(genes_all, cells_all)))
  }
  i_new <- match(rownames(mat)[sm$i], genes_all)
  j_new <- match(colnames(mat)[sm$j], cells_all)
  sparseMatrix(
    i = i_new, j = j_new, x = sm$x,
    dims = c(length(genes_all), length(cells_all)),
    dimnames = list(genes_all, cells_all)
  )
}

# --- 3) Align genes, expand columns, and fill missing values with zeros ---
merge_gene_celltype_tables <- function(x, y,
                                       gene_col = NULL,
                                       dup_cols = c("suffix", "prefer_x", "prefer_y"),
                                       suffix = c("x", "y"),
                                       out = c("dgCMatrix", "data.frame")) {
  dup_cols <- match.arg(dup_cols)
  out <- match.arg(out)
  
  mx <- .to_sparse_gene_cell(x, gene_col = gene_col)
  my <- .to_sparse_gene_cell(y, gene_col = gene_col)
  
  # Align rows to the union of genes
  genes_all <- union(rownames(mx), rownames(my))
  
  # Resolve duplicate cell-type column names before combining
  cx <- colnames(mx)
  cy <- colnames(my)
  inter_cols <- intersect(cx, cy)
  
  if (dup_cols == "suffix" && length(inter_cols) > 0) {
    colnames(mx)[match(inter_cols, cx)] <- paste0(inter_cols, "__", suffix[1])
    colnames(my)[match(inter_cols, cy)] <- paste0(inter_cols, "__", suffix[2])
  } else if (dup_cols == "prefer_x" && length(inter_cols) > 0) {
    my <- my[, setdiff(cy, inter_cols), drop = FALSE]
  } else if (dup_cols == "prefer_y" && length(inter_cols) > 0) {
    mx <- mx[, setdiff(cx, inter_cols), drop = FALSE]
  }
  
  # Combine column names after dropping or renaming conflicting columns
  cells_all_x <- colnames(mx)
  cells_all_y <- colnames(my)
  
  # Expand both matrices to the same genes and fill missing rows with zeros
  mx2 <- .expand_sparse(mx, genes_all, cells_all_x)
  my2 <- .expand_sparse(my, genes_all, cells_all_y)
  
  merged <- cbind(mx2, my2)  # Combine columns and fill missing entries with zeros
  
  if (out == "data.frame") {
    # Converting to a data frame creates a dense object and can use substantial memory
    merged_df <- as.data.frame(as.matrix(merged))
    return(merged_df)
  } else {
    return(as(merged, "dgCMatrix"))
  }
}




draw_GO_picture <- function(
    gene_set,
    top_n = 10,
    select_description = NULL,
    ontology = c("BP", "CC", "MF"),
    min_gene = 1,
    desc_word_n = 15,
    return_dataframe = FALSE,
    y_text_size = 14
) {
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  
  
  # ============================================================
  # 1. Preprocess genes and map SYMBOL to ENTREZID
  # ============================================================
  
  gene_set <- unique(trimws(as.character(gene_set)))
  gene_set <- gene_set[!is.na(gene_set) & gene_set != ""]
  
  if (length(gene_set) == 0) {
    stop("gene_set 为空。")
  }
  
  gene_df <- bitr(
    gene_set,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    stop("没有基因成功从 SYMBOL 转换为 ENTREZID。")
  }
  
  gene <- unique(gene_df$ENTREZID)
  
  
  # ============================================================
  # 2. GO enrichment
  # ============================================================
  
  ontology <- toupper(ontology)
  
  if (!all(ontology %in% c("BP", "CC", "MF"))) {
    stop("ontology 只能包含 BP、CC、MF。")
  }
  
  ego_list <- lapply(ontology, function(ont) {
    
    ego <- enrichGO(
      gene = gene,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    df <- as.data.frame(ego)
    
    if (nrow(df) == 0) {
      return(NULL)
    }
    
    df$Ontology <- ont
    
    return(df)
  })
  
  names(ego_list) <- ontology
  
  go_all <- bind_rows(ego_list)
  
  if (nrow(go_all) == 0) {
    stop("没有显著 GO enrichment 结果。")
  }
  
  
  # ============================================================
  # 3. Filter by enriched gene count
  # ============================================================
  
  go_df <- go_all %>%
    filter(Count >= min_gene)
  
  if (nrow(go_df) == 0) {
    stop("min_gene 过滤后没有剩余 GO term。")
  }
  
  
  # ============================================================
  # 4. Select GO terms to plot
  # ============================================================
  
  if (!is.null(select_description)) {
    
    select_description <- unique(
      trimws(as.character(select_description))
    )
    
    not_found <- setdiff(
      select_description,
      go_df$Description
    )
    
    if (length(not_found) > 0) {
      warning(
        "以下 Description 没有匹配到：\n",
        paste(not_found, collapse = "\n")
      )
    }
    
    go_df <- go_df %>%
      filter(Description %in% select_description)
    
    if (nrow(go_df) == 0) {
      stop("select_description 中没有任何 GO term 被匹配到。")
    }
    
    # Keep the specified order
    go_df$Description_order <- match(
      go_df$Description,
      select_description
    )
    
    go_df <- go_df %>%
      arrange(Description_order)
    
  } else {
    
    # For each ontology, select the top_n terms ranked by p.adjust
    go_df <- go_df %>%
      group_by(Ontology) %>%
      arrange(p.adjust, .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }
  
  
  # ============================================================
  # 5. Convert GeneRatio to numeric
  # ============================================================
  
  ratio_to_numeric <- function(x) {
    
    x <- strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    )[[1]]
    
    if (length(x) != 2) {
      return(NA_real_)
    }
    
    numerator <- suppressWarnings(as.numeric(x[1]))
    denominator <- suppressWarnings(as.numeric(x[2]))
    
    if (
      is.na(numerator) ||
      is.na(denominator) ||
      denominator == 0
    ) {
      return(NA_real_)
    }
    
    numerator / denominator
  }
  
  go_df$GeneRatio_num <- vapply(
    go_df$GeneRatio,
    ratio_to_numeric,
    numeric(1)
  )
  
  
  # ============================================================
  # 6. Process p.adjust values
  # ============================================================
  
  go_df$p.adjust <- as.numeric(go_df$p.adjust)
  
  # Prevent zeros from breaking the log10 scale
  go_df$p.adjust_plot <- pmax(
    go_df$p.adjust,
    1e-300
  )
  
  # Keep -log10(p.adjust) for inspection
  go_df$neglog10_padj <- -log10(
    go_df$p.adjust_plot
  )
  
  
  # ============================================================
  # 7. Shorten descriptions
  # ============================================================
  
  truncate_description <- function(x, n) {
    
    words <- strsplit(
      trimws(x),
      "\\s+"
    )[[1]]
    
    paste(
      head(words, n),
      collapse = " "
    )
  }
  
  go_df$Label <- vapply(
    go_df$Description,
    truncate_description,
    character(1),
    n = desc_word_n
  )
  
  # Keep GO term labels distinct after truncation
  go_df$Label <- make.unique(go_df$Label)
  
  
  # ============================================================
  # 8. Ontology names
  # ============================================================
  
  go_df$Ontology <- factor(
    go_df$Ontology,
    levels = c("BP", "CC", "MF"),
    labels = c(
      "Biological Process",
      "Cellular Component",
      "Molecular Function"
    )
  )
  
  
  # ============================================================
  # 9. GO term display order
  # ============================================================
  
  go_df$Label <- factor(
    go_df$Label,
    levels = rev(unique(go_df$Label))
  )
  
  
  # ============================================================
  # 10. Shared plot theme
  # ============================================================
  
  my_theme <- theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_text(
        size = y_text_size
      ),
      strip.text.y = element_text(
        angle = 0,
        face = "bold"
      ),
      panel.grid.minor = element_blank()
    )
  
  
  # ============================================================
  # 11. Plot 1: GeneRatio
  # ============================================================
  
  ratio_plot <- ggplot(
    go_df,
    aes(
      x = GeneRatio_num,
      y = Label
    )
  ) +
    geom_point(
      aes(
        size = Count,
        color = neglog10_padj
      ),
      alpha = 0.9
    ) +
    facet_grid(
      Ontology ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    labs(
      x = "Gene Ratio",
      y = NULL,
      color = "-log10(p.adjust)",
      size = "Gene Number"
    ) +
    my_theme
  
  
  # ============================================================
  # 12. Plot 2: -log10(p.adjust) bar chart
  # ============================================================
  
  # Set a separate order for the significance chart
  # Place smaller adjusted P values and larger -log10 values at the top
  padj_df <- go_df %>%
    filter(!is.na(neglog10_padj)) %>%
    arrange(neglog10_padj)
  
  padj_df$Label_padj <- factor(
    as.character(padj_df$Label),
    levels = as.character(padj_df$Label)
  )
  
  
  padj_plot <- ggplot(
    padj_df,
    aes(
      x = neglog10_padj,
      y = Label_padj,
      fill = Ontology
    )
  ) +
    geom_col(
      width = 0.7
    ) +
    scale_fill_manual(
      values = c(
        "Biological Process" = "#FD8D62",
        "Cellular Component" = "#66C3A5",
        "Molecular Function" = "#8DA1CB"
      )
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0, 0.05)
      )
    ) +
    labs(
      x = "-log10(adjusted P value)",
      y = NULL,
      fill = "Ontology"
    ) +
    my_theme  
  # ============================================================
  # 13. Return results
  # ============================================================
  
  if (return_dataframe) {
    
    return(
      list(
        ratio_plot = ratio_plot,
        padj_plot = padj_plot,
        go_df = go_df,
        go_all = go_all,
        gene_mapping = gene_df
      )
    )
    
  } else {
    
    return(
      list(
        ratio_plot = ratio_plot,
        padj_plot = padj_plot
      )
    )
  }
}




####Script_81_get_organoid_seurat####

## Build the organoid Seurat object

data_dir <- "/mnt/server1/data9/zouqb/Spatial_endometrium/organoid_24_72_hour/IM-24h-20260203T025011Z-3-001/IM-24h/"

counts <- Read10X(data.dir = data_dir)  # Read10X
if (is.list(counts)) counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]

organoid_24h <- CreateSeuratObject(counts = counts, project = "IM_24h", min.cells = 3, min.features = 200)  # CreateSeuratObject
organoid_24h <- SCTransform(organoid_24h)
organoid_24h <- get_neighbors_and_clusters(organoid_24h)
View_umap(organoid_24h, cluster = 'SCT_snn_res.0.4')
my_Dotplot_plot(seurat_obj = organoid_24h, marker_gene = c('POU5F1', 'NANOG', 'SOX2', 'TDGF1', 'LEFTY1',                  ### Embryonic epiblast
                                                           'TEAD4', 'CDX2', 'GATA3', 'TFAP2C', 'KRT7', 'CLDN4',        ### Embryonic trophectoderm
                                                           'GATA6', 'SOX17', 'PDGFRA',                                 ### Embryonic primitive endoderm
                                                           'EPCAM', 'KRT8', 'KRT18', 'KRT19',                          ### Endometrial epithelium
                                                           'PAX8', 'FOXA2', 'PAEP', 'AGR2', 'TFF3',                    ### Glandular or secretory cells
                                                           'FOXJ1', 'PIFO', 'TPPP3',                                   ### Ciliated cells
                                                           'DCN', 'LUM', 'COL1A1', 'THY1'), cluster = 'SCT_snn_res.0.4')  ## Stroma




saveRDS(organoid_24h, 'integration_sample/26_02_03_organoid_24h.rds')
organoid_24h <- readRDS('integration_sample/26_02_03_organoid_24h.rds')



data_dir <- "/mnt/server1/data9/zouqb/Spatial_endometrium/organoid_24_72_hour/IM-72h-20260203T025039Z-3-001/IM-72h/"

counts <- Read10X(data.dir = data_dir)  # Read10X
if (is.list(counts)) counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]

organoid_72h <- CreateSeuratObject(counts = counts, project = "IM_72h", min.cells = 3, min.features = 200)  # CreateSeuratObject
saveRDS(organoid_72h, 'integration_sample/26_02_03_organoid_72h.rds')
organoid_72h <- readRDS('integration_sample/26_02_03_organoid_72h.rds')

whole_sec_to_use_with_RIF <- readRDS('integration_sample/25_08_08_whole_sec_to_use_with_RIF_new.rds')



### Merge
organoid_24h_cellname <- paste0(colnames(organoid_24h), '_24h')
organoid_72h_cellname <- paste0(colnames(organoid_72h), '_72h')
organoid_cellname <- c(organoid_24h_cellname, organoid_72h_cellname)
organoid_Seurat <- merge(organoid_24h, organoid_72h, add.cell.ids = c('organoid_24h', 'organoid_72h'))


## Initial clustering
organoid_Seurat <- SCTransform(organoid_Seurat)
organoid_Seurat <- get_neighbors_and_clusters(organoid_Seurat)

View_umap(organoid_Seurat, cluster = 'SCT_snn_res.0.7')
View_umap(organoid_Seurat, cluster = 'orig.ident')

saveRDS(organoid_Seurat, 'integration_sample/26_02_03_organoid_Seurat.rds')
organoid_Seurat <- readRDS('integration_sample/26_02_03_organoid_Seurat.rds')

## Correct batch effects
organoid_Seurat <- RunHarmony(organoid_Seurat, group.by.vars = c('orig.ident'))

organoid_Seurat <- RunUMAP(organoid_Seurat, reduction = "harmony", dims = 1:20, reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
organoid_Seurat <- FindNeighbors(organoid_Seurat, reduction = 'harmony', dims = 1:30, graph.name = c('harmony_orig_nn', 'harmony_orig_snn'))
organoid_Seurat <- FindClusters(organoid_Seurat, resolution = seq(0.3, 0.9, 0.1), graph.name = 'harmony_orig_snn')
View_umap(organoid_Seurat, cluster = 'orig.ident', my_umap = 'harmony_umap')
View_umap(organoid_Seurat, cluster = 'harmony_orig_snn_res.0.6', my_umap = 'harmony_umap')

saveRDS(organoid_Seurat, 'integration_sample/26_02_03_organoid_Seurat.rds')
organoid_Seurat <- readRDS('integration_sample/26_02_03_organoid_Seurat.rds')




organoid_Seurat <- PrepSCTFindMarkers(organoid_Seurat)

Idents(organoid_Seurat) <- 'harmony_orig_snn_res.0.6'
aaaa <- FindAllMarkers(object = organoid_Seurat, assay = 'SCT', min.pct = 0.1, only.pos = TRUE)
aaa = aaaa[order(aaaa$avg_log2FC, decreasing = T),]
aaa = aaa[order(aaa$cluster),]
# Select the top num most differentially expressed genes in each group
aaa %>% group_by(cluster) %>% slice_max(n=10, order_by = avg_log2FC) -> top_num

aaa <- aaa %>% 
  group_by(cluster) %>% 
  top_n(15, avg_log2FC)

pdf(paste0('Only_merge_plot/26_02_04_top_10_Heatmap.pdf'), width=12, height=16)
DoHeatmap(organoid_Seurat, label=TRUE, features=aaa$gene) + NoLegend() + theme(axis.text.y = element_text(size = 4)) 
dev.off()

write.csv(aaa, 'marker_gene/26_02_04_organoid_Seurat_marker_gene_top_15.csv', row.names = FALSE)

#### Annotate cells
VlnPlot(object = organoid_Seurat, features = 'nCount_RNA', group.by = 'harmony_orig_snn_res.0.6')

my_Dotplot_plot(seurat_obj = organoid_Seurat, marker_gene = c('POU5F1', 'NANOG', 'SOX2', 'TDGF1', 'LEFTY1', 
                                                              'TEAD4', 'CDX2', 'GATA3', 'TFAP2C', 'KRT7', 'CLDN4',
                                                              'GATA6', 'SOX17', 'PDGFRA',
                                                              'EPCAM', 'KRT8', 'KRT18', 'KRT19',
                                                              'PAX8', 'FOXA2', 'PAEP', 'AGR2', 'TFF3',
                                                              'FOXJ1', 'PIFO', 'TPPP3',
                                                              'DCN', 'LUM', 'COL1A1', 'THY1'), cluster = 'harmony_orig_snn_res.0.6')

my_Dotplot_plot(seurat_obj = organoid_Seurat, marker_gene = c('EPCAM', 'KRT8', 'CDH1', 
                                                              'VIM', 'DCN', 'COL1A1', 'LUM',
                                                              'KRT7', 'GATA3', 'TFAP2C',
                                                              'HLA_G', 'MMP2', 'STB',
                                                              'FOXJ1', 'TPPP3', 'PIFO',
                                                              'PARP', 'SPP1', 'WFDC2', 
                                                              'LGR5'), cluster = 'harmony_orig_snn_res.0.6')



#### Start with the 24 h post-implantation organoid data
PO7_seurat_stroma <- subset(Stroma_and_Gland_Seurat_for_ArchR, subset = period_super_cell_type =='sec_PO7_stroma')
PO9_seurat_stroma <- subset(Stroma_and_Gland_Seurat_for_ArchR, subset = period_super_cell_type =='sec_PO9_stroma')
PO11_seurat_stroma <- subset(Stroma_and_Gland_Seurat_for_ArchR, subset = period_super_cell_type =='sec_PO11_stroma')

organoid_24h_stroma <- subset(organoid_24h, subset = SCT_snn_res.0.4 %in% c('0','1'))
organoid_stroma <- subset(organoid_Seurat, subset = super_cell_type == 'stroma')

DefaultAssay(PO7_seurat_stroma) <- 'RNA'
DefaultAssay(PO9_seurat_stroma) <- 'RNA'
DefaultAssay(PO11_seurat_stroma) <- 'RNA'
DefaultAssay(organoid_24h_stroma) <- 'RNA'
DefaultAssay(organoid_stroma) <- 'RNA'

PO7_seurat_stroma[['SCT']] <- NULL
PO9_seurat_stroma[['SCT']] <- NULL
PO11_seurat_stroma[['SCT']] <- NULL
organoid_24h_stroma[['SCT']] <- NULL
organoid_stroma[['SCT']] <- NULL


PO7_9_11_and_organoid_stroma <- merge(PO7_seurat_stroma, y = c(PO9_seurat_stroma, PO11_seurat_stroma, organoid_stroma))
PO7_9_11_and_organoid_stroma[['RNA']] <- split(PO7_9_11_and_organoid_stroma[['RNA']], f = PO7_9_11_and_organoid_stroma$orig.ident)

PO7_9_11_and_organoid_stroma <- change_annotation(PO7_9_11_and_organoid_stroma, change_list = list('IM_24h' = 'organoid', 'IM_72h' = 'organoid', 
                                                                                                   'sec_PO7_E16_250419_Add' = 'human', 'sec_PO7_E21_240321_RNA' = 'human', 'sec_PO7_E3_230707_RNA' = 'human', 'sec_PO7_E132_250419_Add' = 'human',
                                                                                                   'sec_PO9_E11_250304_RNA' = 'human', 'sec_PO9_E20_250304_RNA' = 'human', 'sec_PO9_E7_250220_RNA' = 'human', 
                                                                                                   'sec_PO11_E15_240413_RNA' = 'human', 'sec_PO11_E13_240629_RNA' = 'human', 'sec_PO11_E1_240629_RNA' = 'human'), old_cluster = 'orig.ident', new_cluster = 'source')

PO7_9_11_and_organoid_stroma <- change_annotation(PO7_9_11_and_organoid_stroma, change_list = list('IM_24h' = 'IM_24H', 'IM_72h' = 'IM_72H', 
                                                                                                   'sec_PO7_E16_250419_Add' = 'sec_PO7', 'sec_PO7_E21_240321_RNA' = 'sec_PO7', 'sec_PO7_E3_230707_RNA' = 'sec_PO7', 'sec_PO7_E132_250419_Add' = 'sec_PO7',
                                                                                                   'sec_PO9_E11_250304_RNA' = 'sec_PO9', 'sec_PO9_E20_250304_RNA' = 'sec_PO9', 'sec_PO9_E7_250220_RNA' = 'sec_PO9', 
                                                                                                   'sec_PO11_E15_240413_RNA' = 'sec_PO11', 'sec_PO11_E13_240629_RNA' = 'sec_PO11', 'sec_PO11_E1_240629_RNA' = 'sec_PO11'), old_cluster = 'orig.ident', new_cluster = 'stage')


Idents(PO7_9_11_and_organoid_stroma) <- 'orig.ident'
DefaultAssay(PO7_9_11_and_organoid_stroma) <- 'RNA'

PO7_9_11_and_organoid_stroma <- NormalizeData(PO7_9_11_and_organoid_stroma, assay = "RNA", )
PO7_9_11_and_organoid_stroma <- FindVariableFeatures(PO7_9_11_and_organoid_stroma, assay = "RNA", )

PO7_9_11_and_organoid_stroma <- SCTransform(PO7_9_11_and_organoid_stroma)
PO7_9_11_and_organoid_stroma <- RunPCA(PO7_9_11_and_organoid_stroma)

PO7_9_11_and_organoid_stroma <- RunHarmony(PO7_9_11_and_organoid_stroma, group.by.vars = c('source', 'orig.ident'))
PO7_9_11_and_organoid_stroma <- RunUMAP(PO7_9_11_and_organoid_stroma, reduction = "harmony", dims = 1:20, reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
PO7_9_11_and_organoid_stroma <- FindNeighbors(PO7_9_11_and_organoid_stroma, reduction = 'harmony', dims = 1:30, graph.name = c('harmony_orig_nn', 'harmony_orig_snn'))
PO7_9_11_and_organoid_stroma <- FindClusters(PO7_9_11_and_organoid_stroma, resolution = seq(0.3, 0.9, 0.1), graph.name = 'harmony_orig_snn')








####Script_82_get_organoid_pre-open####

DefaultAssay(PO7_9_11_and_organoid_stroma) <- 'RNA'
# PO7_9_11_and_organoid_stroma <- JoinLayers(object = PO7_9_11_and_organoid_stroma)
Expression_0_1_RNA_matrix_Organoid <- PO7_9_11_and_organoid_stroma@assays$RNA$counts
Expression_0_1_RNA_matrix_Organoid@x[Expression_0_1_RNA_matrix_Organoid@x > 0] <- 1

## Add the new assay
PO7_9_11_and_organoid_stroma[["RNA_binary"]] <- CreateAssayObject(counts = Expression_0_1_RNA_matrix_Organoid)

avg_list_binary_Organoid_RNA <- AverageExpression(
  object  = PO7_9_11_and_organoid_stroma,
  assays  = "RNA_binary",   # Select an assay, such as "SCT"
  slot    = "counts",  # Select a slot. Usually use "data" for log-normalized values.
  group.by = "period_super_cell_type"  # Column name in meta.data
)

avg_list_binary_Organoid_RNA <- avg_list_binary_Organoid_RNA$RNA_binary
colnames(avg_list_binary_Organoid_RNA) <- gsub(pattern = '-', replacement = '_', x = colnames(avg_list_binary_Organoid_RNA))



# Build a list of cycle and organoid data
avg_merged_list  <- merge_gene_celltype_tables(
  avg_list_binary,
  avg_list_binary_Organoid_RNA,
  dup_cols = "prefer_x",
  out      = "dgCMatrix"
)

avg_merged_list_change_PO9_11_13_to_organoid <- avg_merged_list

# Replace selected columns with organoid values
avg_merged_list_change_PO9_11_13_to_organoid[, "sec_PO9_stroma"]  <- avg_merged_list_change_PO9_11_13_to_organoid[, "IM_24h_stroma"]
avg_merged_list_change_PO9_11_13_to_organoid[, "sec_PO11_stroma"] <- avg_merged_list_change_PO9_11_13_to_organoid[, "IM_72h_stroma"]
avg_merged_list_change_PO9_11_13_to_organoid[, "sec_PO13_stroma"] <- avg_merged_list_change_PO9_11_13_to_organoid[, "IM_72h_stroma"]

# Drop the original organoid columns. Map sec_PO9 to IM_24h and sec_PO11 and sec_PO13 to IM_72h.
avg_merged_list_change_PO9_11_13_to_organoid <- avg_merged_list_change_PO9_11_13_to_organoid[
  , !colnames(avg_merged_list_change_PO9_11_13_to_organoid) %in% c("IM_24h_stroma", "IM_72h_stroma"),
  drop = FALSE
]



ga <- getGeneAnnotation(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)$genes          # GRanges of gene annotations
pk <- getPeakSet(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)                       # GRanges of the current peak set



ga <- getGeneAnnotation(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)$genes
gene_rep <- get_representative_tss_by_symbol(ga, verbose = T)

res_peak_gene_1Mb_organoid <- make_peak_gene_window_matrix(
  proj         = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  gene_symbols = ga.sel$symbol,   # Candidate gene list
  half_window  = 1000000,           # 1,000 kb upstream and downstream, for 2 Mb of flanking sequence
  gene_rep = gene_rep
)

PeakGene_TSSdist_1Mb_organoid <- make_peak_gene_tss_distance_matrix(
  proj          = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
  peak_gene_mat = res_peak_gene_1Mb_organoid$mat,
  gene_rep = res_peak_gene_1Mb_organoid$gene_rep
)
PeakGene_1Mb_mat_aligned_organoid <- res_peak_gene_1Mb_organoid$mat[
  rownames(PeakGene_TSSdist_1Mb_organoid),
  colnames(PeakGene_TSSdist_1Mb_organoid),
  drop = FALSE
]




group_order <- c("RIF_PO7_stroma", "sec_PO3_stroma","sec_PO5_stroma","sec_PO7_stroma","sec_PO9_stroma","sec_PO11_stroma","sec_PO13_stroma",
                 "RIF_PO7_epithelial", "sec_PO3_epithelial","sec_PO5_epithelial","sec_PO7_epithelial","sec_PO9_epithelial","sec_PO11_epithelial","sec_PO13_epithelial")


# 
# Use the neural network to score organoid peak–gene links
#









####Script_83_basic_information_of_organoid_pre-open####


tri_Model_result_stroma_1Mb_organoid <- pre_active_of_stroma_organoid   ## Use the resulting scores

PO7_9_11_and_organoid_stroma <- JoinLayers(PO7_9_11_and_organoid_stroma)
PO7_9_11_and_organoid_stroma <- NormalizeData(PO7_9_11_and_organoid_stroma)
gene_more_expr_in_organoid_vs_PO7 <- FindMarkers(object = PO7_9_11_and_organoid_stroma, ident.1 = 'sec_PO7_stroma', ident.2 = 'IM_24h_stroma', group.by = 'period_super_cell_type', min.pct = 0.1, 
                                                 logfc.threshold = 0, assay = 'RNA')


gene_more_expr_in_organoid <- rownames(gene_more_expr_in_organoid_vs_PO7[gene_more_expr_in_organoid_vs_PO7$avg_log2FC < -0.5 & gene_more_expr_in_organoid_vs_PO7$p_val_adj < 0.01, ])
gene_more_expr_in_PO7 <- rownames(gene_more_expr_in_organoid_vs_PO7[gene_more_expr_in_organoid_vs_PO7$avg_log2FC > 0.5 & gene_more_expr_in_organoid_vs_PO7$p_val_adj < 0.01, ])


tri_Model_result_stroma_1Mb_organoid_filtered <- tri_Model_result_stroma_1Mb_organoid[tri_Model_result_stroma_1Mb_organoid$gene %in% gene_more_expr_in_organoid, ]





library(dplyr)
library(ggplot2)
library(patchwork)


judge_first_open <- function(peak){
  idx <- which(peak >= max(peak) * 0.5)[1]
  return(idx)
}


tri_Model_result_stroma_1Mb_organoid_filtered$first_open <- NA
tri_Model_result_stroma_1Mb_organoid_filtered$first_expr <- NA
peak_to_check <- tri_Model_result_stroma_1Mb_organoid_filtered$peak
gene_to_check <- tri_Model_result_stroma_1Mb_organoid_filtered$gene

for(i in seq(1,dim(tri_Model_result_stroma_1Mb_organoid_filtered)[1])){
  tri_Model_result_stroma_1Mb_organoid_filtered$first_open[i] <- judge_first_open(Y_sig[peak_to_check[i], sec_stroma_cols])
  
}

for(i in seq(1,dim(tri_Model_result_stroma_1Mb_organoid_filtered)[1])){
  tri_Model_result_stroma_1Mb_organoid_filtered$first_expr[i] <- judge_first_open(avg_merged_list_change_PO9_11_13_to_organoid[gene_to_check[i], sec_stroma_cols])
  
}

tri_Model_result_stroma_1Mb_organoid_filtered$first_expr[tri_Model_result_stroma_1Mb_organoid_filtered$gene %in% c("PTK7", "TCF7L1", "PTPRJ", "BNC2", 'COL27A1', 'SLC6A6')] <- 4
tri_Model_result_stroma_1Mb_organoid_filtered <- tri_Model_result_stroma_1Mb_organoid_filtered[tri_Model_result_stroma_1Mb_organoid_filtered$gene != 'SLC7A1',]
tri_Model_result_stroma_1Mb_organoid_filtered$first_open[tri_Model_result_stroma_1Mb_organoid_filtered$first_open == '4'] <- 3



my_Dotplot_plot(seurat_obj = PO7_9_11_and_organoid_stroma, marker_gene = unique(tri_Model_result_stroma_1Mb_organoid_filtered[tri_Model_result_stroma_1Mb_organoid_filtered$first_expr == '1',]$gene),
                cluster = 'stage')

my_Dotplot_plot(seurat_obj = new_Stroma_and_Gland_Seurat_for_ArchR, marker_gene = unique(tri_Model_result_stroma_1Mb_organoid_filtered[tri_Model_result_stroma_1Mb_organoid_filtered$first_expr == '1',]$gene),
                cluster = 'period_super_cell_type')


stroma_first_open_unique <- tri_Model_result_stroma_1Mb_organoid_filtered %>%
  distinct(peak, first_open)

stroma_first_expr_unique <- tri_Model_result_stroma_1Mb_organoid_filtered %>%
  distinct(gene, first_expr)




make_count_df <- function(epi_vec, stroma_vec, metric_name) {
  
  epi_df <- as.data.frame(table(factor(epi_vec, levels = 1:6)))
  stroma_df <- as.data.frame(table(factor(stroma_vec, levels = 1:6)))
  
  colnames(epi_df) <- c("pos", "count")
  colnames(stroma_df) <- c("pos", "count")
  
  epi_df$cell_type <- "epi"
  stroma_df$cell_type <- "stroma"
  
  bind_rows(epi_df, stroma_df) %>%
    mutate(
      metric = metric_name,
      pos = factor(pos, levels = 1:6),
      cell_type = factor(cell_type, levels = c("epi", "stroma"))
    )
}

df_first_open <- make_count_df(
  epi_first_open_unique$first_open,
  stroma_first_open_unique$first_open,
  "first_open"
)



df_first_expr <- make_count_df(
  epi_first_expr_unique$first_expr,
  stroma_first_expr_unique$first_expr,
  "first_expr"
)


p_first_open_stroma <- df_first_open %>%
  filter(cell_type == "stroma") %>%
  ggplot(aes(x = pos, y = count)) +
  geom_col(
    width = 0.7,
    fill = '#D95F5F'
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_open in stroma",
    x = "first_open",
    y = "Count"
  ) +
  theme_classic(base_size = 14)


p_first_expr_stroma <- df_first_expr %>%
  filter(cell_type == "stroma") %>%
  ggplot(aes(x = pos, y = count)) +
  geom_col(
    width = 0.7,
    fill = '#1E90FF'
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_expr in stroma",
    x = "first_expr",
    y = "Count"
  ) +
  theme_classic(base_size = 14)

p_first_open_stroma
p_first_expr_stroma




p_peak_no_x <- p_first_open_stroma +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_combined <- p_peak_no_x / p_first_expr_stroma +
  plot_layout(
    heights = c(1, 1)
  )

p_combined





ggsave(
  filename = "Final_Plot/Fig_4/Fig_4A_organoid_first_open_and_expr.pdf",
  plot = p_combined,
  width = 6,
  height = 4.5,
  units = "in"
)









####Script_84_organoid_pre-open_example####
View(tri_Model_result_stroma_1Mb_organoid_filtered[tri_Model_result_stroma_1Mb_organoid_filtered$gene %in% setdiff(unique(tri_Model_result_stroma_1Mb_organoid_filtered$gene), unique(tri_Model_result_stroma_1Mb$gene)) 
                                                   & tri_Model_result_stroma_1Mb_organoid_filtered$distance <50000 & tri_Model_result_stroma_1Mb_organoid_filtered$peak %notin% tri_Model_result_stroma_1Mb$peak 
                                                   & tri_Model_result_stroma_1Mb_organoid_filtered$gene %notin% tri_Model_result_stroma_1Mb$gene, ])

plot_peak_gene_with_organoid(
  peak_to_use = 'chr22:39212636-39213136',
  gene_to_check = 'PDGFB',
  group_order = c(sec_stroma_cols)
)


gene_to_plot <- c( "PDGFB", "DUSP1", "CTSD", "PLAT")
pdf(paste0('Final_Plot/Fig_4/Fig_4B_organoid_pre_open_peak_gene.pdf'), height = 8, width = 20)
for(i in gene_to_plot){
  gene_to_check <- i
  
  # pdf(paste0('ATAC_plot/', gene_to_check, '_track2.pdf'), height = 15, width = 20)
  p <- plotBrowserTrack(
    ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
    groupBy = "period_super_cell_type",
    useGroups = sec_stroma_cols,
    geneSymbol = gene_to_check,
    upstream = 50000,
    downstream = 50000, 
    maxCells = 10000,
    ylim = c(0, 1)
  )
  
  .draw_any(p)
}
dev.off()

p1 <- VlnPlot(object = new_Stroma_and_Gland_Seurat_for_ArchR, features = 'PDGFB', pt.size = 0, assay = 'RNA', group.by = 'period_super_cell_type'
              , idents = sec_stroma_cols)+
  scale_x_discrete(
    limits = c("sec_PO3_stroma", "sec_PO5_stroma", "sec_PO7_stroma", "sec_PO9_stroma", "sec_PO11_stroma", "sec_PO13_stroma")
  )
p2 <- VlnPlot(object = PO7_9_11_and_organoid_stroma, features = 'PDGFB', pt.size = 0, assay = 'RNA', group.by = 'stage', idents = c("IM_24H", "IM_72H"))+
  scale_x_discrete(
    limits = c( "IM_24H", "IM_72H")
  )







####Script_85_classify_pre-open_gene_and_get_enrichment_result####

pre_active_gene_of_sec <- unique(tri_Model_result_stroma_1Mb$gene)
pre_active_gene_of_organoid <- unique(tri_Model_result_stroma_1Mb_organoid_filtered$gene)
whole_pre_active_gene <- unique(c(pre_active_gene_of_sec, pre_active_gene_of_organoid))

DefaultAssay(PO7_9_11_and_organoid_stroma) <- 'RNA'
PO7_9_11_and_organoid_stroma <- JoinLayers(PO7_9_11_and_organoid_stroma, assay = 'RNA')
PO7_9_11_and_organoid_stroma <- NormalizeData(PO7_9_11_and_organoid_stroma)
## Find genes pre-open only in the nonpregnant condition. These are in pre_active_gene_of_sec (genes_A).
pre_active_gene_of_organoid_gene_A <- FindMarkers(object = PO7_9_11_and_organoid_stroma, ident.1 = 'sec_PO11', ident.2 = 'IM_24H', 
                                                  group.by = 'stage', assay = 'RNA', features = pre_active_gene_of_sec, min.pct = 0.1)

pre_active_gene_of_organoid_gene_A$pct_diff <- pre_active_gene_of_organoid_gene_A$pct.1 - pre_active_gene_of_organoid_gene_A$pct.2

pre_active_gene_more_expr_in_PO11 <- pre_active_gene_of_organoid_gene_A[pre_active_gene_of_organoid_gene_A$pct_diff > 0.2 & 
                                                                          pre_active_gene_of_organoid_gene_A$avg_log2FC > 0.5 & 
                                                                          pre_active_gene_of_organoid_gene_A$pct.2 < 0.3  &
                                                                          pre_active_gene_of_organoid_gene_A$pct.1 > 0.3  &
                                                                          pre_active_gene_of_organoid_gene_A$p_val_adj < 0.05, ]




genes_A <- (rownames(pre_active_gene_more_expr_in_PO11))
genes_B <- intersect(setdiff(pre_active_gene_of_sec, genes_A), rownames(pre_active_gene_of_organoid_gene_A)[pre_active_gene_of_organoid_gene_A$pct.1 > 0.3])
genes_C <- setdiff(whole_pre_active_gene, c(genes_A, genes_B))






if(1){
  gene_list <- list(
    Class_A = genes_A,
    Class_B = genes_B,
    Class_C = genes_C
  )
  genes_use <- unique(unlist(gene_list))
  # =========================
  # 2. Extract the expression matrix
  # =========================
  obj <- PO7_9_11_and_organoid_stroma
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
  genes_present <- genes_use[genes_use %in% rownames(expr_mat)]
  expr_mat <- expr_mat[genes_present, , drop = FALSE]
  # =========================
  # 3. Get group assignments
  # =========================
  stage_vec <- obj@meta.data$stage
  # stage_levels <- c("IM_24H", "IM_72H", "sec_PO11", "sec_PO9", "sec_PO7")
  stage_levels <- c("sec_PO7", "IM_24H", "IM_72H", "sec_PO11")
  # =========================
  # 4. For each stage and gene, calculate the fraction of cells expressing the gene
  # =========================
  pct_exp <- sapply(stage_levels, function(st) {
    cells_use <- colnames(obj)[stage_vec == st]
    Matrix::rowMeans(expr_mat[, cells_use, drop = FALSE] > 0)
  })
  pct_exp <- as.matrix(pct_exp)
  colnames(pct_exp) <- stage_levels
  avg_exp <- pct_exp
  # =========================
  # 5. Prepare the matrix
  # =========================
  mat <- avg_exp
  mat[is.na(mat)] <- 0
  # Keep genes present in the matrix and order them by group
  gene_list_present <- lapply(gene_list, function(g) g[g %in% rownames(mat)])
  # Drop groups with no genes in the expression matrix
  gene_list_present <- gene_list_present[lengths(gene_list_present) > 0]
  gene_order <- unlist(gene_list_present, use.names = FALSE)
  mat <- mat[gene_order, , drop = FALSE]
  # =========================
  # 6. Build row annotations
  # =========================
  group_labels <- rep(names(gene_list_present), times = lengths(gene_list_present))
  ann_row <- data.frame(
    Group = factor(group_labels, levels = names(gene_list_present))
  )
  rownames(ann_row) <- gene_order
  
  
  # Insert gaps between groups
  gaps_row <- cumsum(lengths(gene_list_present))
  gaps_row <- gaps_row[-length(gaps_row)]  # Omit a trailing gap
  # Optional colors
  ann_colors <- list(
    Group = c(
      Class_A = "#E41A1C",
      Class_B = "#377EB8",
      Class_C = "#4DAF4A"
    )
  )
  # =========================
  # 7. Draw the heatmap
  # =========================
  p <- pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    gaps_row = gaps_row,
    border_color = NA,
    fontsize_row = 7,
    fontsize_col = 10
  )
}


pdf(
  "Final_Plot/Fig_4/Fig_4C_heatmap_of_pre_open_gene.pdf",
  width = 8,
  height = 7,
  useDingbats = FALSE
)

grid::grid.newpage()
grid::grid.draw(p$gtable)

dev.off()



draw_GO_picture(gene_set = genes_A, select_description = c('immune response-regulating signaling pathway', 'positive regulation of leukocyte activation', 'positive regulation of leukocyte proliferation', 'cell adhesion mediated by integrin', 'cell-matrix adhesion'))
draw_GO_picture(gene_set = genes_B, select_description = c('extracellular matrix organization', 'regulation of cell-substrate adhesion', 'tissue remodeling', 'fibroblast migration', 'collagen metabolic process'))
draw_GO_picture(gene_set = genes_C, select_description = c('cytoplasmic translation', 'cellular oxidant detoxification', 'regulation of apoptotic signaling pathway', 'intrinsic apoptotic signaling pathway', 'focal adhesion'))






####Script_86_Cellchat_and_Nitchnet_analysis####

Blastaea_and_endometrium_and_organoid_stroma <- readRDS('integration_sample/26_03_16_Blastaea_and_endometrium_and_organoid_stroma.rds')
# pre_active_gene_of_organoid <- readRDS('marker_gene/26_03_16_organoid_preactive_gene.rds')
pre_active_gene_of_organoid <- unique(tri_Model_result_stroma_1Mb_organoid_filtered$gene)



organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}






## CellChat analysis

options(future.globals.maxSize = 20 * 1024^3)

Idents(seurat_obj) <- 'period_super_cell_type'

Blastaea_and_organoid_cellchat_sample <- createCellChat(object = seurat_obj,  group.by = "period_super_cell_type", assay = "SCT")



# Load the ligand–receptor database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
# Select a database subset
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB)
Blastaea_and_organoid_cellchat_sample@DB <- CellChatDB.use

## Extract signaling genes into data.signaling
Blastaea_and_organoid_cellchat_sample <- subsetData(Blastaea_and_organoid_cellchat_sample)

# Identify overexpressed ligands and receptors in each group, similar to Seurat FindMarkers
Blastaea_and_organoid_cellchat_sample <- identifyOverExpressedGenes(Blastaea_and_organoid_cellchat_sample)
Blastaea_and_organoid_cellchat_sample <- identifyOverExpressedInteractions(Blastaea_and_organoid_cellchat_sample) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
# The previous step stores results in cellchat@LR$LRsig
Blastaea_and_organoid_cellchat_sample <- smoothData(Blastaea_and_organoid_cellchat_sample, adj = PPI.human) # projectData projects ligand-receptor expression in @data.signaling onto the PPI network and saves adjusted values in @data.project


# Infer interaction probabilities from expression. This is the slowest step.
Blastaea_and_organoid_cellchat_sample <- computeCommunProb(Blastaea_and_organoid_cellchat_sample, raw.use = FALSE, population.size = TRUE) # Set raw.use = TRUE to use unadjusted expression
Blastaea_and_organoid_cellchat_sample <- computeCommunProbPathway(Blastaea_and_organoid_cellchat_sample)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
Blastaea_and_organoid_cellchat_sample <- filterCommunication(Blastaea_and_organoid_cellchat_sample, min.cells = 3)


pairLR.use_for_surface_stroma <- data.frame(
  interaction_name = c("SEMA3E_PLXND1", "SEMA3D_NRP2_PLXNA2", 'SEMA3D_NRP1_PLXNA2', 'NRXN1_NLGN1', 'NRXN1_CLSTN2', 'NRXN1_CLSTN1')
)


pdf(paste0('Final_Plot/Fig_4/Fig_4E_the_communication_only_TE_can_trans.pdf'), height = 5, width = 7)
netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('sec_PO7_stroma'), sources.use = c('TE_TE', 'sec_PO7_stroma', 'sec_PO7_epithelial', 'sec_PO7_endothelial', 'sec_PO7_immune', 'sec_PO7_luminal'), 
                 pairLR.use = pairLR.use_for_surface_stroma, sort.by.source = TRUE)  # Hierarchy plot
dev.off()



pdf(paste0('R_cellchat/cellchat_plot/', Sys.Date(), '_which_one_is_TE_have_but_other_dont_have_original.pdf'), height = 50, width = 8)
netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('sec_PO7_stroma'), sources.use = c('TE_TE', 'sec_PO7_stroma', 'sec_PO7_epithelial', 'sec_PO7_endothelial', 'sec_PO7_immune', 'sec_PO7_luminal'), 
                 sort.by.source = TRUE)  # Hierarchy plot
dev.off()







## NicheNet analysis
library(Seurat)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(nichenetr)
library(CellChat)

wd='/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/'
setwd(wd)


## Load reference inputs
seurat_obj <- readRDS('integration_sample/26_03_16_Blastaea_and_endometrium_and_organoid_stroma.rds')
PO7_9_11_and_organoid_stroma <- readRDS('integration_sample/26_02_26_PO7_9_11_and_organoid_stroma_using_whole_sec_to_use_with_RIF.rds')
pre_active_gene_of_organoid <- readRDS('marker_gene/26_07_27_organoid_preactive_gene.rds')

pre_active_gene_of_organoid <- readRDS('Final_sample_and_doc/Fig_4/26_09_07_pre_active_gene_of_organoid.rds')
pre_active_gene_of_organoid <- unique(pre_active_gene_of_organoid$gene)
organism <- "human"


lr_network <- readRDS('NitchNet/human_lr_network.rds')
ligand_target_matrix <- readRDS('NitchNet/human_ligand_target_matrix.rds')
weighted_networks <- readRDS('NitchNet/human_weighted_networks.rds')

lr_network <- lr_network %>% distinct(from, to)


## Define sender and receiver groups
Idents(seurat_obj) <- "period_super_cell_type"
Idents(PO7_9_11_and_organoid_stroma) <- "stage"
sender_celltype <- "TE"
# receiver_celltypes <- c("IM_24h_stroma", "IM_72h_stroma")
receiver_celltypes <- c("IM_24h_stroma")
receiver_cells <- WhichCells(seurat_obj, idents = receiver_celltypes)
sender_cells   <- WhichCells(seurat_obj, idents = sender_celltype)


## Define expressed receiver genes and expressed sender ligands
get_expressed_genes_custom <- function(celltype_names, seurat_obj, pct = 0.10, assay = "RNA") {
  cells <- WhichCells(seurat_obj, idents = celltype_names)
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = "data")[, cells, drop = FALSE]
  expressed_genes <- rownames(expr_mat)[Matrix::rowSums(expr_mat > 0) / ncol(expr_mat) >= pct]
  return(expressed_genes)
}


seurat_obj <- JoinLayers(object = seurat_obj, assay = 'RNA')
seurat_obj <- NormalizeData(seurat_obj, assay = 'RNA')

PO7_9_11_and_organoid_stroma <- JoinLayers(object = PO7_9_11_and_organoid_stroma, assay = 'RNA')
PO7_9_11_and_organoid_stroma <- NormalizeData(PO7_9_11_and_organoid_stroma, assay = 'RNA')

# expressed_genes_receiver <- get_expressed_genes_custom(
#   celltype_names = receiver_celltypes,
#   seurat_obj = seurat_obj,
#   pct = 0.10,
#   assay = "SCT"
# )
# 
# expressed_genes_sender <- get_expressed_genes_custom(
#   celltype_names = sender_celltype,
#   seurat_obj = seurat_obj,
#   pct = 0.10,
#   assay = "SCT"
# )
# 
expressed_genes_receiver <- get_expressed_genes_custom(
  celltype_names = receiver_celltypes,
  seurat_obj = seurat_obj,
  pct = 0.10,
  assay = "RNA"
)


expressed_genes_sender <- get_expressed_genes_custom(
  celltype_names = sender_celltype,
  seurat_obj = seurat_obj,
  pct = 0.10,
  assay = "RNA"
)

# 


## Find ligands expressed by senders and receptors expressed by receivers
lr_network <- lr_network %>% distinct(from, to)

potential_ligands <- lr_network %>%
  filter(from %in% expressed_genes_sender,
         to %in% expressed_genes_receiver) %>%
  pull(from) %>%
  unique()


## Define the target gene set
moreExprInOrganoid <- FindMarkers(object = PO7_9_11_and_organoid_stroma, ident.1 = 'IM_24H', ident.2 = 'sec_PO7', group.by = 'stage', min.pct = 0.1, assay = 'RNA')
moreExprInOrganoid <- moreExprInOrganoid[moreExprInOrganoid$avg_log2FC > 0.7, ]

more_expr_in_PO11_vs_PO7_of_stroma <- FindMarkers(object = PO7_9_11_and_organoid_stroma, ident.1 = 'sec_PO11', ident.2 = 'sec_PO7', group.by = 'stage', min.pct = 0.1)
more_expr_in_PO11_vs_PO7_of_stroma <- more_expr_in_PO11_vs_PO7_of_stroma[more_expr_in_PO11_vs_PO7_of_stroma$avg_log2FC > 0.7,]

geneset_oi <- intersect(pre_active_gene_of_organoid, rownames(moreExprInOrganoid))

# geneset_oi <- intersect(rownames(moreExprInOrganoid), expressed_genes_receiver)
# geneset_oi <- intersect(rownames(moreExprInOrganoid), expressed_genes_receiver)
background_expressed_genes <- expressed_genes_receiver


# geneset_oi <- intersect(rownames(more_expr_in_PO11_vs_PO7_of_stroma), expressed_genes_receiver)
# background_expressed_genes <- get_expressed_genes_custom(
#   celltype_names = 'sec_PO11',
#   seurat_obj = PO7_9_11_and_organoid_stroma,
#   pct = 0.10,
#   assay = "SCT"
# )


### Calculate ligand activity
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected))


head(ligand_activities, 20)


top_ligands <- ligand_activities %>%
  slice_max(order_by = aupr_corrected, n = 30) %>%
  pull(test_ligand)


top_lr  <- lr_network %>%
  filter(from %in% top_ligands,
         to %in% expressed_genes_receiver)
top_lr


### Identify candidate pre_active target genes for each ligand
top_ligands <- ligand_activities %>%
  slice_max(
    order_by = aupr_corrected,
    n = 30,
    with_ties = FALSE
  ) %>%
  pull(test_ligand)

active_ligand_target_links_df <- top_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 100
  ) %>%
  bind_rows() %>%
  tidyr::drop_na() %>%
  group_by(ligand) %>%
  arrange(desc(weight), .by_group = TRUE) %>%
  mutate(
    target_rank_among_retained = row_number()
  ) %>%
  ungroup()


active_ligand_target_links_df



## Run steps 4–10 first

top_ligand_activity <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%
  slice_head(n = 30) %>%
  mutate(
    test_ligand = factor(test_ligand, levels = rev(test_ligand))
  ) 

top_ligands <- ligand_activities %>%
  slice_max(order_by = aupr_corrected, n = 30) %>%
  pull(test_ligand)



top_lr <- lr_network %>%
  filter(from %in% top_ligands,
         to %in% expressed_genes_receiver)

ligands_keep <- c("SEMA3D", "SEMA3E", "NRXN1")

top_lr <- lr_network %>%
  filter(
    from %in% ligands_keep,
    to %in% expressed_genes_receiver
  )

top_lr_plot_df <- top_lr %>%
  mutate(
    from = factor(from, levels = rev(top_ligands)),
    to = factor(to, levels = unique(to))
  )



p_top_lr <- ggplot(top_lr_plot_df, aes(x = to, y = from)) +
  geom_point(size = 3, color = "#0072B2") +
  theme_classic(base_size = 13) +
  labs(
    x = "Receptor expressed in receiver stroma",
    y = "Top ligand expressed in TE",
    title = "Predicted ligand-receptor links"
    # subtitle = "Filtered by ligand expression in TE and receptor expression in stroma"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold")
  )



pdf('Final_Plot/Fig_S4/Fig_S4D_NitchNet_NRXN_SEMA3_pathway_about_pre_active_gene.pdf', width = 4, height = 3)
print(p_top_lr)
dev.off()







