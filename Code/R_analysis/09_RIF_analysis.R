#### Helper functions ####


### Plot cell-type proportions as bars
my_bar_plot <- function(
    seurat_obj,
    clu_1 = "stage",
    clu_2 = "cell_type",
    clu_1_order = NULL,
    clu_2_order = NULL,
    my_colors = my_colors2
) {
  
  meta_data <- seurat_obj@meta.data
  
  # Set the x-axis order
  if (!is.null(clu_1_order)) {
    meta_data[[clu_1]] <- factor(meta_data[[clu_1]], levels = clu_1_order)
  }
  
  # Set the cell-type order
  if (!is.null(clu_2_order)) {
    meta_data[[clu_2]] <- factor(meta_data[[clu_2]], levels = clu_2_order)
  }
  
  # Calculate cluster proportions within each stage
  cluster_stage_counts <- meta_data %>%
    group_by(across(all_of(c(clu_1, clu_2)))) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(across(all_of(clu_1))) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  p <- ggplot(
    cluster_stage_counts,
    aes(x = .data[[clu_1]], y = percentage, fill = .data[[clu_2]])
  ) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      y = "Percentage",
      x = clu_1,
      title = "Cluster Distribution by Stage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Use clu_2_order for the legend order when supplied
  if (!is.null(clu_2_order)) {
    p <- p +
      scale_fill_manual(
        values = my_colors,
        breaks = clu_2_order,
        limits = clu_2_order
      )
  } else {
    p <- p +
      scale_fill_manual(values = my_colors)
  }
  
  return(p)
}










####Script_91_get_sec_PO7_and_RIF_PO7_seurat####
sec_PO7_and_RIF_seurat <- subset(pro_and_sec_seurat, subset = stage_2 %in% c('sec_PO7', 'RIF_PO7'))


####
# Set the default assay to RNA
DefaultAssay(sec_PO7_and_RIF_seurat) <- "RNA"
sec_PO7_and_RIF_seurat[['SCT']] <- NULL
# Join layers first so all layers are used
sec_PO7_and_RIF_seurat <- JoinLayers(sec_PO7_and_RIF_seurat)
# Extract counts from the counts layer
counts_data <- LayerData(sec_PO7_and_RIF_seurat, assay = "RNA", layer = "counts")
# Create a new counts-only assay
new_assay <- CreateAssay5Object(counts = counts_data)
# Replace the original assay
sec_PO7_and_RIF_seurat[["RNA"]] <- new_assay
# Verify that only the counts layer remains
Layers(sec_PO7_and_RIF_seurat, assay = "RNA")
# Split the layers again
sec_PO7_and_RIF_seurat[['RNA']] <- split(sec_PO7_and_RIF_seurat[['RNA']], f = sec_PO7_and_RIF_seurat$orig.ident)


sec_PO7_and_RIF_seurat <- SCTransform(sec_PO7_and_RIF_seurat)
sec_PO7_and_RIF_seurat <- RunPCA(sec_PO7_and_RIF_seurat)
sec_PO7_and_RIF_seurat <- RunHarmony(sec_PO7_and_RIF_seurat, group.by.vars = c('id', 'exper_method'))
sec_PO7_and_RIF_seurat <- RunUMAP(sec_PO7_and_RIF_seurat, reduction = "harmony", dims = 1:20, reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
sec_PO7_and_RIF_seurat <- FindNeighbors(sec_PO7_and_RIF_seurat, reduction = 'harmony', dims = 1:30, graph.name = c('harmony_orig_nn', 'harmony_orig_snn'))
sec_PO7_and_RIF_seurat <- FindClusters(sec_PO7_and_RIF_seurat, resolution = seq(0.3, 0.9, 0.1), graph.name = 'harmony_orig_snn')



p1 <- DimPlot(
  sec_PO7_and_RIF_seurat,
  reduction = "harmony_umap",
  group.by = "stage", 
  label = TRUE,
  repel = TRUE
)

p2 <- DimPlot(
  sec_PO7_and_RIF_seurat,
  reduction = "harmony_umap",
  group.by = "exclusive_cell_type_4",
  label = TRUE,
  repel = TRUE
)



RIF_and_normal_sample_order <- c(
  "sec_PO7",
  "RIF_PO7_E6_230815_AddRNA",
  "RIF_PO7_E35_231127_AddRNA",
  "RIF_PO7_E59_240312_AddRNA",
  "RIF_PO7_E15_240321_AddRNA",
  "RIF_PO7_E10_241226_Add",
  "RIF_PO7_E20_RNA_250925"
)
cell_type_order <- c(
  "luminal",
  "ciliated",
  "gland",
  "inner_stroma",
  "surface_stroma",
  "lymphatic_endothelial",
  "blood_endothelial",
  "Perivascular_cell",
  "lymphoid",
  "myeloid"
)

p <- my_bar_plot(seurat_obj = sec_PO7_and_RIF_seurat, clu_1 = 'stage', clu_2 = 'exclusive_cell_type_4',  clu_1_order = RIF_and_normal_sample_order, clu_2_order = cell_type_order)
ggsave(
  filename = "Final_Plot/Fig_5/Fig_5A_celltype_distribution_of_sec_and_RIF.pdf",
  plot = p,
  width = 5,
  height = 9,
  units = "in"
)


p1 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'sec_PO7_E16_250419_Add')
p2 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'RIF_PO7_E35_231127_AddRNA')

p3 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'sec_PO7_E16_250419_Add', Colors = my_colors_stroma_3)
p4 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'RIF_PO7_E35_231127_AddRNA', my_colors_stroma_3)

p5 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'sec_PO7_E16_250419_Add', Colors = my_colors_stroma_2)
p6 <- my_spatial_plot(seurat_object = sec_PO7_and_RIF_seurat, cluster = 'exclusive_cell_type_4', sample_name = 'RIF_PO7_E35_231127_AddRNA', my_colors_stroma_2)

ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_sec_PO7.pdf", plot = p1, width = 7, height = 7, units = "in")
ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_RIF_PO7_E35_231127_AddRNA.pdf", plot = p2, width = 7, height = 7, units = "in")

ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_sec_PO7_surface_stroma.pdf", plot = p3, width = 7, height = 7, units = "in")
ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_RIF_PO7_E35_231127_AddRNA_surface_stroma.pdf", plot = p4, width = 7, height = 7, units = "in")

ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_sec_PO7_luminal.pdf", plot = p5, width = 7, height = 7, units = "in")
ggsave(filename = "Final_Plot/Fig_5/Fig_5B_spatial_distribution_of_RIF_PO7_E35_231127_AddRNA_luminal.pdf", plot = p6, width = 7, height = 7, units = "in")








####Script_92_immune-dysfunctional_RIF####

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
sec_PO7_and_RIF_seurat <- readRDS('Final_sample_and_doc/Fig_5/26_06_22_sec_PO7_and_RIF_seurat.rds')





aaa <- table(sec_PO7_and_RIF_seurat$exclusive_cell_type_4, sec_PO7_and_RIF_seurat$orig.ident)
saveRDS(aaa, 'Final_sample_and_doc/Supplementary_table_and_information/Fig_5B_celltype_number_of_normal_and_RIF_PO7.rds')

lymphoid_ratio <- sec_PO7_and_RIF_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    total_cells = n(),
    lymphoid_cells = sum(exclusive_cell_type_4 == "lymphoid", na.rm = TRUE),
    lymphoid_ratio = lymphoid_cells / total_cells,
    lymphoid_percent = lymphoid_ratio * 100,
    .groups = "drop"
  )

normal <- c(0.0571, 0.055, 0.0731, 0.0873)
abnormal <- 0.138

# 2. Apply the arcsine square-root transformation
normal_arcsin <- asin(sqrt(normal))
abnormal_arcsin <- asin(sqrt(abnormal))


# 3. Crawford–Howell modified t-test
n <- length(normal_arcsin)

t_stat <- (abnormal_arcsin - mean(normal_arcsin)) /
  (sd(normal_arcsin) * sqrt((n + 1) / n))

p_value <- 2 * pt(-abs(t_stat), df = n - 1)

p_label <- paste0("p = ", sprintf("%.4f", p_value))

# 4. Prepare the bar chart data
bar_df <- data.frame(
  group = factor(
    c("Normal", "Abnormal"),
    levels = c("Normal", "Abnormal")
  ),
  mean_ratio = c(mean(normal), abnormal),
  sd_ratio = c(sd(normal), NA_real_)
)

# 5. Prepare sample-level points
point_df <- data.frame(
  group = factor(
    c(rep("Normal", length(normal)), "Abnormal"),
    levels = c("Normal", "Abnormal")
  ),
  ratio = c(normal, abnormal)
)

# 6. Draw the plot
p <- ggplot(bar_df, aes(x = group, y = mean_ratio,
                        fill = group)) +
  
  # Bars
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.4
  ) +
  
  # Normal group error bars (mean ± SD)
  geom_errorbar(
    data = filter(bar_df, group == "Normal"),
    aes(
      ymin = mean_ratio - sd_ratio,
      ymax = mean_ratio + sd_ratio
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual sample points
  geom_jitter(
    data = point_df,
    aes(x = group, y = ratio),
    inherit.aes = FALSE,
    width = 0.08,
    height = 0,
    size = 2.5,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8
  ) +
  
  # Significance bars
  annotate("segment",
           x = 1, xend = 2,
           y = 0.18, yend = 0.18) +
  annotate("segment",
           x = 1, xend = 1,
           y = 0.176, yend = 0.18) +
  annotate("segment",
           x = 2, xend = 2,
           y = 0.176, yend = 0.18) +
  
  # P values
  annotate(
    "text",
    x = 1.5, y = 0.188,
    label = p_label,
    size = 5
  ) +
  
  # Colors
  scale_fill_manual(
    values = c(
      "Normal" = "#76B7B2",
      "Abnormal" = "#FF61C3"
    )
  ) +
  
  # Percentage scale on the y-axis
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 0.21),
    expand = expansion(mult = c(0, 0))
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  labs(
    y = "Lymphoid cell proportion"
  )

print(p)

# 7. Save a PDF for editing in Adobe Illustrator
ggsave(
  "Final_Plot/Fig_5/Fig_5B_lymphoid_celltype_proportion_test.pdf",
  plot = p,
  width = 4,
  height = 4,
  device = cairo_pdf
)











pro_and_sec_and_Blastaea_seurat_for_cellchat$period_super_cell_type <- paste0(pro_and_sec_and_Blastaea_seurat_for_cellchat$stage, '_', pro_and_sec_and_Blastaea_seurat_for_cellchat$exclusive_cell_type_4)
DefaultAssay(pro_and_sec_and_Blastaea_seurat_for_cellchat) <- 'SCT'

pro_and_sec_and_Blastaea_seurat_for_cellchat <- JoinLayers(
  object = pro_and_sec_and_Blastaea_seurat_for_cellchat,
  assay = "RNA"
)
pro_and_sec_and_Blastaea_seurat_for_cellchat <- NormalizeData(  pro_and_sec_and_Blastaea_seurat_for_cellchat, assay = "RNA")

pro_and_sec_and_Blastaea_seurat_for_cellchat$period_super_cell_type <- paste0(pro_and_sec_and_Blastaea_seurat_for_cellchat$stage, '_', pro_and_sec_and_Blastaea_seurat_for_cellchat$super_cell_type)
Idents(pro_and_sec_and_Blastaea_seurat_for_cellchat) <- 'period_super_cell_type'

Blastaea_cellchat_sample <- createCellChat(object = pro_and_sec_and_Blastaea_seurat_for_cellchat,  group.by = "period_super_cell_type", assay = "SCT")


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

saveRDS(Blastaea_cellchat_sample, 'R_cellchat/cellchat_sample/26_08_04_pro_and_sec_and_Blastaea_seurat_for_cellchat_using_SCT_period_super_cell_type.rds')





pdf(paste0('Final_Plot/Fig_5/Fig_5C_Blastaea_to_PO7_and_RIF_stroma_to_communication.pdf'), height = 6)
netVisual_bubble(Blastaea_cellchat_sample, remove.isolate = FALSE, targets.use = c('TE_TE'), sources.use = c('sec_PO7_surface_stroma', 'sec_PO11_surface_stroma', 'RIF_PO7_E6_230815_AddRNA_surface_stroma', 'RIF_PO7_E35_231127_AddRNA_surface_stroma', 'RIF_PO7_E59_240312_AddRNA_surface_stroma', 'RIF_PO7_E15_240321_AddRNA_surface_stroma', 'RIF_PO7_E10_241226_Add_surface_stroma', 'RIF_PO7_E20_RNA_250925_surface_stroma'),
                 signaling = c('PARs'), sort.by.source = TRUE)  # Hierarchy plot
dev.off()


diff_gene_between_RIF_PO7_E6_and_sec_PO7 <- FindMarkers(object = sec_PO7_and_RIF_seurat, ident.1 = 'sec_PO7_stroma', ident.2 = 'RIF_PO7_E6_230815_AddRNA_stroma', group.by = 'period_super_cell_type', min.pct = 0.1)

more_expr_in_RIF_PO7_E6_vs_PO7 <- rownames(diff_gene_between_RIF_PO7_E6_and_sec_PO7)[diff_gene_between_RIF_PO7_E6_and_sec_PO7$avg_log2FC < -0.7 & diff_gene_between_RIF_PO7_E6_and_sec_PO7$p_val_adj < 0.01]


bbb <- draw_GO_picture(more_expr_in_RIF_PO7_E6_vs_PO7, select_description = c('antigen receptor-mediated signaling pathway', 'positive regulation of leukocyte activation', 'immune response-activating cell surface receptor signaling pathway', 'positive regulation of leukocyte cell-cell adhesion', 'regulation of leukocyte migration'))

ggsave(
  filename = "Final_Plot/Fig_5/Fig_5D_GO_enrichment_of_diff_gene.pdf",
  plot = bbb$padj_plot,
  width = 11,
  height = 3,
  units = "in"
)




## Test immune cell counts with the Crawford–Howell test
aaa <- table(sec_PO7_and_RIF_seurat$exclusive_cell_type_4, sec_PO7_and_RIF_seurat$orig.ident)
saveRDS(aaa, 'Final_sample_and_doc/Supplementary_table_and_information/Fig_5B_celltype_number_of_normal_and_RIF_PO7.rds')

lymphoid_ratio <- sec_PO7_and_RIF_seurat@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    total_cells = n(),
    lymphoid_cells = sum(exclusive_cell_type_4 == "lymphoid", na.rm = TRUE),
    lymphoid_ratio = lymphoid_cells / total_cells,
    lymphoid_percent = lymphoid_ratio * 100,
    .groups = "drop"
  )

normal <- c(0.0571, 0.055, 0.0731, 0.0873)
abnormal <- 0.138

# 2. Apply the arcsine square-root transformation
normal_arcsin <- asin(sqrt(normal))
abnormal_arcsin <- asin(sqrt(abnormal))


# 3. Crawford–Howell modified t-test
n <- length(normal_arcsin)

t_stat <- (abnormal_arcsin - mean(normal_arcsin)) /
  (sd(normal_arcsin) * sqrt((n + 1) / n))

p_value <- 2 * pt(-abs(t_stat), df = n - 1)

p_label <- paste0("p = ", sprintf("%.4f", p_value))

# 4. Prepare the bar chart data
bar_df <- data.frame(
  group = factor(
    c("Normal", "Abnormal"),
    levels = c("Normal", "Abnormal")
  ),
  mean_ratio = c(mean(normal), abnormal),
  sd_ratio = c(sd(normal), NA_real_)
)

# 5. Prepare sample-level points
point_df <- data.frame(
  group = factor(
    c(rep("Normal", length(normal)), "Abnormal"),
    levels = c("Normal", "Abnormal")
  ),
  ratio = c(normal, abnormal)
)

# 6. Draw the plot
p <- ggplot(bar_df, aes(x = group, y = mean_ratio,
                        fill = group)) +
  
  # Bars
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.4
  ) +
  
  # Normal group error bars (mean ± SD)
  geom_errorbar(
    data = filter(bar_df, group == "Normal"),
    aes(
      ymin = mean_ratio - sd_ratio,
      ymax = mean_ratio + sd_ratio
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  
  # Individual sample points
  geom_jitter(
    data = point_df,
    aes(x = group, y = ratio),
    inherit.aes = FALSE,
    width = 0.08,
    height = 0,
    size = 2.5,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 0.8
  ) +
  
  # Significance bars
  annotate("segment",
           x = 1, xend = 2,
           y = 0.18, yend = 0.18) +
  annotate("segment",
           x = 1, xend = 1,
           y = 0.176, yend = 0.18) +
  annotate("segment",
           x = 2, xend = 2,
           y = 0.176, yend = 0.18) +
  
  # P values
  annotate(
    "text",
    x = 1.5, y = 0.188,
    label = p_label,
    size = 5
  ) +
  
  # Colors
  scale_fill_manual(
    values = c(
      "Normal" = "#76B7B2",
      "Abnormal" = "#FF61C3"
    )
  ) +
  
  # Percentage scale on the y-axis
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 0.21),
    expand = expansion(mult = c(0, 0))
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  labs(
    y = "Lymphoid cell proportion"
  )

print(p)

# 7. Save a PDF for editing in Adobe Illustrator
ggsave(
  "Final_Plot/Fig_5/Fig_5B_lymphoid_celltype_proportion_test.pdf",
  plot = p,
  width = 4,
  height = 4,
  device = cairo_pdf
)







####Script_93_implantation_windows_shift_RIF####

pre_active_peak_of_organoid <- unique(pre_active_of_stroma_organoid$peak)
pre_active_peak_of_sec <- unique(pre_active_of_stroma$peak)
Sync_active_peak_of_sec_stroma <- unique(Sync_of_stroma$peak)
Sync_active_peak_of_sec_stroma_of_High_in_sec_late <- unique(Sync_of_stroma[Sync_of_stroma$gene %in% High_in_sec_late_gene_of_stroma, ]$peak)
Sync_active_peak_of_sec_stroma_of_High_in_sec_early <- unique(Sync_of_stroma[Sync_of_stroma$gene %in% High_in_sec_early_gene_of_stroma, ]$peak)


whole_pre_active_peak <- unique(c(pre_active_peak_of_sec, pre_active_peak_of_organoid))
length(whole_pre_active_peak)

pre_active_of_sec_open_in_PO7 <- pre_active_of_stroma[pre_active_of_stroma$first_open == 3, ]
pre_active_of_sec_open_in_PO9 <- pre_active_of_stroma[pre_active_of_stroma$first_open == 4, ]

pre_active_peak_of_sec_open_in_PO7 <- unique(pre_active_of_sec_open_in_PO7$peak)
pre_active_peak_of_sec_open_in_PO9 <- unique(pre_active_of_sec_open_in_PO9$peak)



plot_distribution(Y_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('RIF_PO7_E10_241226_Add_stroma')] - Y_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('sec_PO7_stroma')])
plot_distribution(Y_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('RIF_PO7_E10_241226_Add_stroma')] - Y_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('sec_PO9_stroma')])

plot_distribution(Y_sample_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('sec_PO7_E16_250419_Add||sec_PO7_stroma')] - Y_sig[Sync_active_peak_of_sec_stroma_of_High_in_sec_late, c('sec_PO7_stroma')])



diff_PO7 <- 
  Y_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "RIF_PO7_E10_241226_Add_stroma"
  ] -
  Y_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "sec_PO7_stroma"
  ]


diff_PO9 <- 
  Y_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "RIF_PO7_E10_241226_Add_stroma"
  ] -
  Y_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "sec_PO9_stroma"
  ]


diff_sec_PO7 <- 
  Y_sample_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "sec_PO7_E16_250419_Add||sec_PO7_stroma"
  ] -
  Y_sig[
    Sync_active_peak_of_sec_stroma_of_High_in_sec_late,
    "sec_PO7_stroma"
  ]




group_colors <- c(
  "RIF_PO7 - sec_PO7" = "#E64B35",
  "RIF_PO7 - sec_PO9" = "#4DBBD5",
  "PO7_sample - sec_PO7" = "#00A087"
)


# =========================
# ECDF 1
# =========================

df_compare <- data.frame(
  value = c(diff_PO7, diff_PO9),
  group = c(
    rep("RIF_PO7 - sec_PO7", length(diff_PO7)),
    rep("RIF_PO7 - sec_PO9", length(diff_PO9))
  )
)

df_compare <- df_compare[is.finite(df_compare$value), ]


p_ecdf_1 <- ggplot(
  df_compare,
  aes(x = value, color = group)
) +
  stat_ecdf(
    linewidth = 1.2,
    geom = "step"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = "Difference",
    y = "Cumulative proportion",
    color = NULL
  ) +
  theme_minimal(base_size = 14)

p_ecdf_1


df_compare <- data.frame(
  value = c(diff_PO7, diff_sec_PO7),
  group = c(
    rep("RIF_PO7 - sec_PO7", length(diff_PO7)),
    rep("PO7_sample - sec_PO7", length(diff_sec_PO7))
  )
)

df_compare <- df_compare[is.finite(df_compare$value), ]



p_ecdf_2 <- ggplot(
  df_compare,
  aes(x = value, color = group)
) +
  stat_ecdf(
    linewidth = 1.2,
    geom = "step"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = "Difference",
    y = "Cumulative proportion",
    color = NULL
  ) +
  theme_minimal(base_size = 14)

p_ecdf_2


ggsave(
  filename = "Final_Plot/Fig_5/Fig_5F_Sync_open_in_sec_late_diff_of_RIFvsPO7_and_RIFvsPO9.pdf",
  plot = p_ecdf_1,
  width = 8,
  height = 3,
  units = "in"
)

ggsave(
  filename = "Final_Plot/Fig_5/Fig_5F_Sync_open_in_sec_late_diff_of_RIFvsPO7_and_PO7sampleVSPO7.pdf",
  plot = p_ecdf_2,
  width = 8,
  height = 3,
  units = "in"
)


diff_PO7_of_pre_active_peak <- 
  Y_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "RIF_PO7_E10_241226_Add_stroma"
  ] -
  Y_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "sec_PO7_stroma"
  ]



diff_PO9_of_pre_active_peak <- 
  Y_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "RIF_PO7_E10_241226_Add_stroma"
  ] -
  Y_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "sec_PO9_stroma"
  ]


diff_sec_PO7_of_pre_active_peak <- 
  Y_sample_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "sec_PO7_E16_250419_Add||sec_PO7_stroma"
  ] -
  Y_sig[
    pre_active_peak_of_sec_open_in_PO9,
    "sec_PO7_stroma"
  ]


group_colors <- c(
  "RIF_PO7 - sec_PO7" = "yellow",
  "PO7_sample - sec_PO7" = "pink"
)



df_compare <- data.frame(
  value = c(diff_PO7_of_pre_active_peak, diff_sec_PO7_of_pre_active_peak),
  group = c(
    rep("RIF_PO7 - sec_PO7", length(diff_PO7_of_pre_active_peak)),
    rep("PO7_sample - sec_PO7", length(diff_sec_PO7_of_pre_active_peak))
  )
)



# Remove NA and infinite values
df_compare <- df_compare[is.finite(df_compare$value), ]



p_ecdf_3 <- ggplot(
  df_compare,
  aes(x = value, color = group)
) +
  stat_ecdf(
    linewidth = 1.2,
    geom = "step"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  labs(
    x = "Difference",
    y = "Cumulative proportion",
    color = NULL
  ) +
  theme_minimal(base_size = 14)

p_ecdf_3




ggsave(
  filename = "Final_Plot/Fig_5/Fig_5F_pre_active_peak_of_sec_open_in_PO9_diff_of_RIFvsPO7_and_PO7sampleVSPO7.pdf",
  plot = p_ecdf_3,
  width = 7,
  height = 3,
  units = "in"
)




