

#### Helper functions ####

## Get the initial nearest-neighbor graph built from PCA
get_neighbors_and_clusters <- function(seurat_obj){
  seurat_obj <- RunPCA(seurat_obj, assay = "SCT", verbose = T)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:k1)
  seurat_obj <- FindClusters(seurat_obj, verbose = T, resolution = seq(0.3, 0.7, 0.1))
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:k1)
  return(seurat_obj)
}

my_Dotplot_plot <- function(seurat_obj, marker_gene, cluster='SCT_snn_res.0.5', title='Marker_gene', idents = NULL){
  if(!is.na(idents)){
    Idents(seurat_obj) <- idents
  }
  p <- DotPlot(seurat_obj, 
               features = marker_gene,  # Use selected marker genes
               group.by = cluster,
               idents = idents) +  # Group by cluster
    ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6))+  # Rotate x-axis labels
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('red','#336699','#66CC66','#FFCC33'))
  return(p)
}





# Project back into spatial coordinates
my_spatial_plot <- function(
    seurat_object,
    cluster = "SCT_snn_res.0.7",
    Colors = NULL,
    sample_name = NULL
) {
  meta <- seurat_object@meta.data
  
  required_cols <- c("orig.ident", "nA", "nB", cluster)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("meta.data 中缺少这些列: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(sample_name)) {
    plot_data <- meta[, required_cols, drop = FALSE]
  } else {
    plot_data <- meta[meta$orig.ident %in% sample_name, required_cols, drop = FALSE]
  }
  
  plot_data$nA <- as.numeric(as.character(plot_data$nA))
  plot_data$nB <- as.numeric(as.character(plot_data$nB))
  plot_data$cluster_label <- as.factor(plot_data[[cluster]])
  
  if (is.null(Colors)) {
    my_colors <- my_colors2
  } else {
    my_colors <- Colors
  }
  
  fixed_points <- data.frame(
    nA = c( 0, 99, 99),
    nB = c( 99, 99, 0)
  )
  
  title_text <- paste(unique(plot_data$orig.ident), collapse = ", ")
  
  p <- ggplot(plot_data, aes(x = nB, y = nA, fill = cluster_label)) +
    geom_tile() +
    scale_fill_manual(values = my_colors) +
    theme_minimal() +
    labs(
      title = title_text,
      x = "nB",
      y = "nA",
      fill = cluster
    ) +
    coord_fixed() +
    geom_point(
      data = fixed_points,
      aes(x = nB, y = nA),
      inherit.aes = FALSE,
      color = "red",
      shape = 16,
      size = 3
    )
  
  return(p)
}



#### Clustering ####
pro_and_sec_seurat <- get_neighbors_and_clusters(pro_and_sec_seurat)
pro_and_sec_seurat_no_RIF <- subset(pro_and_sec_seurat, subset = stage != 'RIF_PO7')

pro_and_sec_seurat_no_RIF <- RunHarmony(pro_and_sec_seurat_no_RIF, group.by.vars = c('exper_method', 'orig.ident'))
pro_and_sec_seurat_no_RIF <- RunUMAP(pro_and_sec_seurat_no_RIF, reduction = "harmony", dims = 1:20, reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
pro_and_sec_seurat_no_RIF <- FindNeighbors(pro_and_sec_seurat_no_RIF, reduction = 'harmony', dims = 1:30, graph.name = c('harmony_orig_nn', 'harmony_orig_snn'))
pro_and_sec_seurat_no_RIF <- FindClusters(pro_and_sec_seurat_no_RIF, resolution = seq(0.3, 0.7, 0.1), graph.name = 'harmony_orig_snn')



## Plot known markers using celltype.ls defined above
pdf(paste0('Cell_type_list_seurat.pdf'))
for (i in celltype.ls){
  p <- my_Dotplot_plot(seurat_obj, marker_gene = unique(get(i)), title = i, cluster = 'harmony_orig_snn_res.0.5')
  print(p)
}
dev.off()


## Plot UMAP
p1 <- DimPlot(
  pro_and_sec_seurat_no_RIF,
  reduction = "harmony_umap",
  group.by = "cell_type_to_plot",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  ) +
  ggtitle("UMAP of Normal Sample") +
  theme(
    # Legend text on the right
    legend.text = element_text(size = 15),
    
    # Adjust the legend title if present
    legend.title = element_text(size = 13),
    legend.key.height = grid::unit(0.65, "cm")
  )



ggsave(
  filename = "Final_Plot/Fig_1/Fig_1B_Harmony_UMAP_of_pro_and_sec_seurat.pdf",
  plot = p1,
  width = 12,
  height = 7.5,
  units = "in"
)


stroma_and_epi_cell <- subset(pro_and_sec_seurat_no_RIF, subset = super_cell_type %in% c('stroma', 'epithelial'))
stroma_and_epi_cell <- SCTransform(object = stroma_and_epi_cell)


p2 <- DimPlot(
  stroma_and_epi_cell,
  reduction = "umap",
  group.by = "stage_period_2_exclusive_cell_type_4",
  label = FALSE,
  repel = TRUE,
  raster = FALSE
)

ggsave(
  filename = "Final_Plot/Fig_1/Fig_1B_UMAP_of_pro_and_sec_seurat_of_stroma_and_gland.pdf",
  plot = p2,
  width = 11,
  height = 7.5,
  units = "in"
)


## Plot spatial projections
p1 <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'cell_type_to_plot', sample_name = 'pro_late_E19_240718_RNA', Colors = my_colors_for_UMAP)
p2 <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'cell_type_to_plot', sample_name = 'sec_PO3_E26_241024_RNA', Colors = my_colors_for_UMAP)
p3 <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'cell_type_to_plot', sample_name = 'sec_PO7_E16_250419_Add', Colors = my_colors_for_UMAP)
p4 <- my_spatial_plot(seurat_object = pro_and_sec_seurat, cluster = 'cell_type_to_plot', sample_name = 'sec_PO11_E15_240413_RNA', Colors = my_colors_for_UMAP)

library(patchwork)
combined_plot <- p1 / p2 / p3 / p4
ggsave(
  "Final_Plot/Fig_1/Fig_1D_Celltype_distribution_in_space.pdf",
  combined_plot,
  width = 8,
  height = 20,
  units = "in"
)


marker_gene_to_plot <- c('CAPS', 'DNAH12', 'MS4A8',
                         'LGR5', 'XACT', 'GDA',
                         'TRPM6', 'HSD17B2', 'DISC1FP1',
                         'TMEM196', 'SLC22A3', 'PLCXD3',
                         'VWF', 'NR5A2', 'ANO2',
                         'MCAM', 'GUCY1A2', 'SULF1',
                         'PTPRC', 'NCAM1', 'KCNQ5',
                         'CSF1R', 'PLXDC2', 'CMAHP')


p1 <- my_Dotplot_plot(pro_and_sec_seurat_no_RIF, marker_gene = marker_gene_to_plot, cluster = 'exclusive_cell_type_factor_2', idents = rev(c('myeloid', 'lymphoid', 'Perivascular_cell', 'blood_endothelial', 'stroma', 'gland', 'luminal', 'ciliated')))

ggsave(
  "Final_Plot/Fig_1/Fig_1E_Celltype_marker_Dotplot.pdf",
  p1,
  width = 11,
  height = 3.5,
  units = "in"
)

