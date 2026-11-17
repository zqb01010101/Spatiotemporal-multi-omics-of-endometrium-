
####Script_71_show_example_case####
pro_and_sec_stroma_and_gland <- subset(pro_and_sec_seurat, subset = super_cell_type %in% c('stroma', 'gland'))

plot_order <- c(
  "sec_PO3_stroma",
  "sec_PO3_epithelial",
  "sec_PO5_stroma",
  "sec_PO5_epithelial",
  "sec_PO7_stroma",
  "sec_PO7_epithelial",
  "sec_PO9_stroma",
  "sec_PO9_epithelial",
  "sec_PO11_stroma",
  "sec_PO11_epithelial",
  "sec_PO13_stroma",
  "sec_PO13_epithelial"
)

groups <- c(
  "sec_PO3_stroma",
  "sec_PO5_stroma",
  "sec_PO7_stroma",
  "sec_PO9_stroma",
  "sec_PO11_stroma",
  "sec_PO13_stroma"
)


pro_and_sec_stroma_and_gland$period_super_cell_type <- factor(
  pro_and_sec_stroma_and_gland$period_super_cell_type,
  levels = plot_order
)


gene_to_plot <- c('DOCK8', 'ITGAX', 'GLRX', 'ITGAD', 'IKZF1', 'CLEC12A', 'ALPL', 
                  'IGFBP1')

gene_to_plot <- c('RDH10', )

VlnPlot(object = pro_and_sec_stroma_and_gland, features = 'BACH1', group.by = 'period_super_cell_type', assay = 'SCT', pt.size = 0, 
        idents = c('sec_PO3_stroma', 'sec_PO5_stroma', 'sec_PO7_stroma', 'sec_PO9_stroma', 'sec_PO11_stroma', 'sec_PO13_stroma'), cols = rep("#4C78A8", 6)) +
  scale_x_discrete(limits = groups) + 
  NoLegend()




sec_stroma_cols <- c("sec_PO3_stroma","sec_PO5_stroma","sec_PO7_stroma","sec_PO9_stroma","sec_PO11_stroma","sec_PO13_stroma")

pdf(paste0('Final_Plot/Fig_3/immune_gene_pre_active_peak.pdf'), height = 8, width = 20)
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


gene_to_plot <- c('CD247', 'ARHGAP20', 'IGFBP1')

p1 <- VlnPlot(object = pro_and_sec_stroma_and_gland, features = 'ARHGAP20', group.by = 'period_super_cell_type', assay = 'RNA', pt.size = 0, 
              idents = c('sec_PO3_stroma', 'sec_PO5_stroma', 'sec_PO7_stroma', 'sec_PO9_stroma', 'sec_PO11_stroma', 'sec_PO13_stroma'), cols = rep("#4C78A8", 6)) +
  scale_x_discrete(limits = groups) + 
  NoLegend()

p2 <- VlnPlot(object = pro_and_sec_stroma_and_gland, features = 'IGFBP1', group.by = 'period_super_cell_type', assay = 'RNA', pt.size = 0, 
              idents = c('sec_PO3_stroma', 'sec_PO5_stroma', 'sec_PO7_stroma', 'sec_PO9_stroma', 'sec_PO11_stroma', 'sec_PO13_stroma'), cols = rep("#4C78A8", 6)) +
  scale_x_discrete(limits = groups) + 
  NoLegend()


ggsave(
  filename = "Final_Plot/Fig_3/Fig_3A_ARHGA20_RNA_Expression.pdf",
  plot = p1,
  width = 6,
  height = 3,
  units = "in"
)


ggsave(
  filename = "Final_Plot/Fig_3/Fig_3A_IGFBP1_RNA_Expression.pdf",
  plot = p2,
  width = 6,
  height = 3,
  units = "in"
)



pdf(paste0('Final_Plot/Fig_3/Fig_3A_ARHGA20_ATAC_Open.pdf'), height = 7, width = 7)
for(i in 1){
  gene_to_check <- 'ARHGAP20'
  
  p <- plotBrowserTrack(
    ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
    groupBy = "period_super_cell_type",
    useGroups = sec_stroma_cols,
    geneSymbol = gene_to_check,
    upstream = 15000,
    downstream = 15000, 
    maxCells = 10000,
    ylim = c(0, 1)
  )
  .draw_any(p)
}
dev.off()


pdf(paste0('Final_Plot/Fig_3/Fig_3A_IGFBP1_ATAC_Open.pdf'), height = 7, width = 7)
for(i in 1){
  gene_to_check <- 'IGFBP1'
  p <- plotBrowserTrack(
    ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
    groupBy = "period_super_cell_type",
    useGroups = sec_stroma_cols,
    geneSymbol = gene_to_check,
    upstream = 15000,
    downstream = 15000, 
    maxCells = 10000,
    ylim = c(0, 1)
  )
  .draw_any(p)
}
dev.off()






####Script_72_basic_information_of_Sync-opening####

Sync_peak_gene_stroma <- Sync_of_stroma
Sync_peak_gene_epi <- Sync_of_epi



print('------------stroma的情况------------')
dim(Sync_peak_gene_stroma)
length(unique((Sync_peak_gene_stroma$gene)))  ## Gene count
length(unique((Sync_peak_gene_stroma$peak)))  ## Peak count

print('------------epi的情况----------')
dim(Sync_peak_gene_epi)
length(unique((Sync_peak_gene_epi$gene)))  ## Gene count
length(unique((Sync_peak_gene_epi$peak)))  ## Peak count




judge_first_open <- function(peak){
  idx <- which(peak >= max(peak) * 0.5)[1]
  return(idx)
}


Sync_peak_gene_stroma$first_open <- NA
Sync_peak_gene_stroma$first_expr <- NA
peak_to_check <- Sync_peak_gene_stroma$peak
gene_to_check <- Sync_peak_gene_stroma$gene

for(i in seq(1,dim(Sync_peak_gene_stroma)[1])){
  Sync_peak_gene_stroma$first_open[i] <- judge_first_open(Y_sig[peak_to_check[i], sec_stroma_cols])
  
}

for(i in seq(1,dim(Sync_peak_gene_stroma)[1])){
  Sync_peak_gene_stroma$first_expr[i] <- judge_first_open(avg_list_binary[gene_to_check[i], sec_stroma_cols])
  
}


Sync_peak_gene_epi$first_open <- NA
Sync_peak_gene_epi$first_expr <- NA
peak_to_check <- Sync_peak_gene_epi$peak
gene_to_check <- Sync_peak_gene_epi$gene

for(i in seq(1,dim(Sync_peak_gene_epi)[1])){
  Sync_peak_gene_epi$first_open[i] <- judge_first_open(Y_sig[peak_to_check[i], sec_epi_cols])
  
}

for(i in seq(1,dim(Sync_peak_gene_epi)[1])){
  Sync_peak_gene_epi$first_expr[i] <- judge_first_open(avg_list_binary[gene_to_check[i], sec_epi_cols])
  
}




epi_first_open_unique <- Sync_peak_gene_epi %>%
  distinct(peak, first_open)

stroma_first_open_unique <- Sync_peak_gene_stroma %>%
  distinct(peak, first_open)

epi_first_expr_unique <- Sync_peak_gene_epi %>%
  distinct(gene, first_expr)

stroma_first_expr_unique <- Sync_peak_gene_stroma %>%
  distinct(gene, first_expr)



library(dplyr)
library(ggplot2)
library(patchwork)

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

p_first_open <- ggplot(df_first_open, aes(x = pos, y = count, fill = cell_type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_open",
    x = "first_open",
    y = "Count",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14)





p_first_expr <- ggplot(df_first_expr, aes(x = pos, y = count, fill = cell_type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_expr",
    x = "first_expr",
    y = "Count",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14)

p_first_open
p_first_expr



ggsave(
  filename = "Final_Plot/Fig_3/Fig_3B_Sync_open_first_open.pdf",
  plot = p_first_open,
  width = 7,
  height = 3,
  units = "in"
)

ggsave(
  filename = "Final_Plot/Fig_3/Fig_3B_Sync_open_first_expr.pdf",
  plot = p_first_expr,
  width = 7,
  height = 3,
  units = "in"
)











####Script_73_classify_Sync_open_gene####
Sync_gene_expr_matrix <- avg_list_binary[unique(Sync_gene_set_stroma), sec_stroma_cols]
Sync_gene_expr_matrix <- avg_list_binary[unique(Sync_gene_set_epi), sec_epi_cols]
mat <- as.matrix(Sync_gene_expr_matrix)
# Scale columns before clustering so large-valued columns do not dominate distances
# mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat


## K-means clustering
set.seed(123)

library(cluster)
k_range <- 2:25
sil_scores <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  km <- kmeans(mat_scaled, centers = k, nstart = 20)
  sil <- silhouette(km$cluster, dist(mat_scaled))
  sil_scores[i] <- mean(sil[, 3])
}

data.frame(k = k_range, silhouette = sil_scores)

plot(k_range, sil_scores, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Mean silhouette width")


# km_res <- kmeans(mat_scaled, centers = 14, nstart = 50)  ## Stroma
km_res <- kmeans(mat_scaled, centers = 16, nstart = 50)  ## Epithelium


# Cluster assignment for each row
cluster_kmeans <- km_res$cluster

table(cluster_kmeans)
head(cluster_kmeans)

result_df <- data.frame(
  element = rownames(mat),
  cluster = cluster_kmeans
)


## Plot temporal patterns
clusters <- cluster_kmeans
# clusters <- res$kmeans
length(clusters)
nrow(Sync_gene_expr_matrix)
library(Matrix)

cluster_ids <- sort(unique(clusters))

avg_expr <- do.call(rbind, lapply(cluster_ids, function(cl) {
  Matrix::colMeans(Sync_gene_expr_matrix[clusters == cl, , drop = FALSE])
}))

rownames(avg_expr) <- paste0("Cluster_", cluster_ids)

avg_expr

cluster_colors <- rainbow(nrow(avg_expr))

for (i in 1:nrow(avg_expr)) {
  plot(
    x = 1:ncol(avg_expr),
    y = avg_expr[i, ],
    type = "b",
    lwd = 2,
    pch = 19,
    col = cluster_colors[i],
    xaxt = "n",
    xlab = "Stage",
    ylab = "Mean expression",
    main = rownames(avg_expr)[i],
    ylim = range(avg_expr)   # Use a shared y-axis range for comparison
  )
  
  axis(1, at = 1:ncol(avg_expr), labels = colnames(avg_expr))
}

# stroma_cluster_kmeans <- cluster_kmeans  
epi_cluster_kmeans <- cluster_kmeans





### Process separately below

Idents(new_Stroma_and_Gland_Seurat_for_ArchR) <- 'period_super_cell_type'
my_Dotplot_plot(seurat_obj = new_Stroma_and_Gland_Seurat_for_ArchR, marker_gene = names(cluster_kmeans[cluster_kmeans == 3]), cluster = 'period_super_cell_type',
                idents = c('sec_PO3_stroma', 'sec_PO5_stroma', 'sec_PO7_stroma', 'sec_PO9_stroma', 'sec_PO11_stroma', 'sec_PO13_stroma'))



stroma_cluster_label <- c(
  "1" = "Low_in_sec_mid",
  "2" = "High_in_sec_late",
  "3" = "High_in_sec_late",    ###Low_expression/Low_in_sec_mid
  "4" = "High_in_sec_late",
  "5" = "Low_in_sec_mid",
  "6" = "High_in_sec_early",    ###Low_expression/Low_in_sec_mid
  "7" = "High_in_sec_late",
  "8" = "Low_in_sec_mid",
  "9" = "High_in_sec_late",
  "10" = "High_in_sec_early",
  "11" = "High_in_sec_early",
  "12" = "High_in_sec_early",
  "13" = "High_in_sec_late",
  "14" = "High_in_sec_late"
)
stroma_cluster_label <- setNames(
  unname(stroma_cluster_label[as.character(stroma_cluster_kmeans)]),
  names(stroma_cluster_kmeans)
)




my_Dotplot_plot(seurat_obj = new_Stroma_and_Gland_Seurat_for_ArchR, marker_gene = names(cluster_kmeans[cluster_kmeans == 1]), cluster = 'period_super_cell_type',
                idents = c('sec_PO3_epithelial', 'sec_PO5_epithelial', 'sec_PO7_epithelial', 'sec_PO9_epithelial', 'sec_PO11_epithelial', 'sec_PO13_epithelial'))


epi_cluster_label <- c(
  "1" = "Low_in_sec_mid",
  "2" = "High_in_sec_early",
  "3" = "High_in_sec_early",
  "4" = "High_in_sec_late",
  "5" = "High_in_sec_late",
  "6" = "High_in_sec_late",
  "7" = "Low_in_sec_mid",   ##?
  "8" = "High_in_sec_late",
  "9" = "Low_in_sec_mid",
  "10" = "High_in_sec_early",
  "11" = "High_in_sec_early",
  "12" = "High_in_sec_late",
  "13" = "High_in_sec_late",
  "14" = "High_in_sec_late",
  "15" = "High_in_sec_early",
  "16" = "High_in_sec_early"
)
epi_cluster_label <- setNames(
  unname(epi_cluster_label[as.character(epi_cluster_kmeans)]),
  names(epi_cluster_kmeans)
)





Sync_gene_expr_matrix <- avg_list_binary[unique(Sync_gene_set_stroma), sec_stroma_cols]
Sync_gene_expr_matrix <- avg_list_binary[unique(Sync_gene_set_epi), sec_epi_cols]


clusters <- stroma_cluster_label
clusters <- epi_cluster_label
# clusters <- res$kmeans
length(clusters)
nrow(Sync_gene_expr_matrix)
library(Matrix)

cluster_ids <- sort(unique(clusters))

avg_expr <- do.call(rbind, lapply(cluster_ids, function(cl) {
  Matrix::colMeans(Sync_gene_expr_matrix[clusters == cl, , drop = FALSE])
}))

rownames(avg_expr) <- paste0(cluster_ids)

avg_expr

cluster_colors <- rainbow(nrow(avg_expr))

matplot(
  x = 1:ncol(avg_expr),
  y = t(avg_expr),
  type = "l",
  lty = 1,
  lwd = 2,
  col = cluster_colors,
  xaxt = "n",
  xlab = "Stage",
  ylab = "Mean expression",
  main = "Mean expression pattern of each cluster"
)

axis(1, at = 1:ncol(avg_expr), labels = colnames(avg_expr))

legend(
  "topright",
  legend = rownames(avg_expr),
  col = cluster_colors,
  lty = 1,
  lwd = 2,
  cex = 0.8
)



pdf('Final_Plot/Fig_3/Fig_3C_three_situation_shape.pdf', width = 4, height = 4)
for (i in 1:nrow(avg_expr)) {
  plot(
    x = 1:ncol(avg_expr),
    y = avg_expr[i, ],
    type = "b",
    lwd = 4,
    pch = 19,
    col = cluster_colors[i],
    xaxt = "n",
    xlab = "Stage",
    ylab = "Mean expression",
    main = rownames(avg_expr)[i],
    ylim = range(avg_expr)   # Use a shared y-axis range for comparison
  )
  
  axis(1, at = 1:ncol(avg_expr), labels = colnames(avg_expr))
}
dev.off()










####Script_74_TF_enrichment_of_Sync-open_peak####

High_in_sec_late_peak_of_Epi <- unique(Sync_of_epi[Sync_of_epi$gene %in% names(epi_cluster_label[epi_cluster_label == 'High_in_sec_late']), ]$peak)
High_in_sec_early_peak_of_Epi <- unique(Sync_of_epi[Sync_of_epi$gene %in% names(epi_cluster_label[epi_cluster_label == 'High_in_sec_early']), ]$peak)
Low_in_sec_mid_peak_of_Epi <- unique(Sync_of_epi[Sync_of_epi$gene %in% names(epi_cluster_label[epi_cluster_label == 'Low_in_sec_mid']), ]$peak)

High_in_sec_late_peak_of_stroma <- unique(Sync_of_stroma[Sync_of_stroma$gene %in% names(stroma_cluster_label[stroma_cluster_label == 'High_in_sec_late']), ]$peak)
High_in_sec_early_peak_of_stroma <- unique(Sync_of_stroma[Sync_of_stroma$gene %in% names(stroma_cluster_label[stroma_cluster_label == 'High_in_sec_early']), ]$peak)
Low_in_sec_mid_peak_of_stroma <- unique(Sync_of_stroma[Sync_of_stroma$gene %in% names(stroma_cluster_label[stroma_cluster_label == 'Low_in_sec_mid']), ]$peak)
peaks_to_gr <- function(x){
  x <- unique(x)
  parts <- do.call(rbind, strsplit(x, "[:-]"))
  chr <- parts[,1]
  chr <- ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
  GRanges(seqnames = chr, ranges = IRanges(as.integer(parts[,2]), as.integer(parts[,3])))
}

## Match peaks to the existing peak set
peak_to_enrich <- list(High_in_sec_late_peak_of_Epi, High_in_sec_early_peak_of_Epi, Low_in_sec_mid_peak_of_Epi, High_in_sec_late_peak_of_stroma, High_in_sec_early_peak_of_stroma, Low_in_sec_mid_peak_of_stroma)
names(peak_to_enrich) <- c('High_in_sec_late_peak_of_Epi', 'High_in_sec_early_peak_of_Epi', 'Low_in_sec_mid_peak_of_Epi', 'High_in_sec_late_peak_of_stroma', 'High_in_sec_early_peak_of_stroma', 'Low_in_sec_mid_peak_of_stroma')

pdf('try.pdf')
for(peak_vec in peak_to_enrich){
  if(1){
    peak_df <- data.frame(
      peak = peak_vec,
      chr = sub(":.*", "", peak_vec),
      start = as.integer(sub(".*:(\\d+)-.*", "\\1", peak_vec)),
      end = as.integer(sub(".*-", "", peak_vec))
    )
    
    peak_gr <- GRanges(
      seqnames = peak_df$chr,
      ranges = IRanges(
        start = peak_df$start,
        end = peak_df$end
      )
    )
    
    peak_gr
    length(peak_gr)
    
    
    pSet <- getPeakSet(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)
    
    hits <- findOverlaps(
      query = peak_gr,
      subject = pSet
    )
    
    length(unique(queryHits(hits)))
    length(peak_gr)
    
    
    ##
    # motif_matches <- getMatches(
    #   ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
    #   name = "Motif"
    # )
    # 
    # bgdPeaks <- getBgdPeaks(new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered)
    
    
    
    motif_enrich <- customEnrichment(
      ranges = peak_gr,
      matches = motif_matches,
      bgdPeaks = bgdPeaks
    )
    
    
    motif_enrich_df <- as.data.frame(motif_enrich)
    
    motif_enrich_df <- motif_enrich_df %>%
      arrange(dplyr::desc(mlog10Padj)) %>%
      mutate(
        padj = 10^(-mlog10Padj),
        pvalue = 10^(-mlog10p)
      )
    
    motif_enrich <- motif_enrich[motif_enrich$mlog10Padj > 3, ]
    
    head(motif_enrich_df, 30)
    # saveRDS(motif_enrich_df, '/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/Final_sample_and_doc/Fig_3/Low_in_sec_mid_peak_of_Epi_motif_enrich_df.rds')
    # saveRDS(motif_enrich_df, '/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/Final_sample_and_doc/Fig_3/Low_in_sec_mid_peak_of_Epi_motif_enrich_df_Motif_version.rds')
    # 
    # View(motif_enrich_df)
    
    motif_enrich_df <- motif_enrich_df[order(motif_enrich_df$mlog10Padj, decreasing = TRUE), ]
    head(motif_enrich_df$feature)
    # View(motif_enrich_df)
    
    
    library(ggplot2)
    library(ggrepel)
    plot_df <- motif_enrich_df %>%
      arrange(desc(mlog10Padj)) %>%
      mutate(rank = row_number())
    
    p <- ggplot(plot_df, aes(x = Enrichment, y = mlog10Padj)) +
      geom_point(size = 1) +
      ggrepel::geom_text_repel(
        data = plot_df %>% filter(rank <= 10), 
        aes(label = feature),
        size = 7,
        
        # Set the label offset
        nudge_x = 0.4,       # Move right
        nudge_y = -0.4,      # Move up
        
        # Increase the distance between labels and points
        point.padding = 1,
        box.padding = 0.8,
        
        # Always draw connector lines from labels to points
        min.segment.length = 0,
        segment.color = "grey40",
        segment.size = 0.4,
        
        max.overlaps = Inf
      ) +
      theme_classic() +
      xlab("Enrichment") +
      ylab("-log10 P value") +
      ggtitle("Motif enrichment of High_in_sec_late_peak_of_Epi")
    
    print(p)
    
  }
}
dev.off()







####Script_75_basic_information_of_pre-open####
judge_first_open <- function(peak){
  idx <- which(peak >= max(peak) * 0.5)[1]
  return(idx)
}


pre_active_of_stroma$first_open <- NA
pre_active_of_stroma$first_expr <- NA
peak_to_check <- pre_active_of_stroma$peak
gene_to_check <- pre_active_of_stroma$gene

for(i in seq(1,dim(pre_active_of_stroma)[1])){
  pre_active_of_stroma$first_open[i] <- judge_first_open(Y_sig[peak_to_check[i], sec_stroma_cols])
  
}

for(i in seq(1,dim(pre_active_of_stroma)[1])){
  pre_active_of_stroma$first_expr[i] <- judge_first_open(avg_list_binary[gene_to_check[i], sec_stroma_cols])
  
}

pre_active_of_stroma$first_expr[pre_active_of_stroma$gene %in% c("CLSTN1", "COL7A1", "SAMHD1", "MEG3", 'LAMC3', 'SPIRE1')] <- 5
pre_active_of_stroma$first_expr[pre_active_of_stroma$gene %in% c("FREM1")] <- 6


pre_active_of_epi$first_open <- NA
pre_active_of_epi$first_expr <- NA
peak_to_check <- pre_active_of_epi$peak
gene_to_check <- pre_active_of_epi$gene

for(i in seq(1,dim(pre_active_of_epi)[1])){
  pre_active_of_epi$first_open[i] <- judge_first_open(Y_sig[peak_to_check[i], sec_epi_cols])
  
}

for(i in seq(1,dim(pre_active_of_epi)[1])){
  pre_active_of_epi$first_expr[i] <- judge_first_open(avg_list_binary[gene_to_check[i], sec_epi_cols])
  
}

pre_active_of_epi$first_expr[pre_active_of_epi$gene %in% c("PSAP", "ABTB2", "TLN2")] <- 5




tri_Model_result_stroma_1Mb <- pre_active_of_stroma
tri_Model_result_epi_1Mb <- pre_active_of_epi



epi_first_open_unique <- tri_Model_result_epi_1Mb %>%
  distinct(peak, first_open)

stroma_first_open_unique <- tri_Model_result_stroma_1Mb %>%
  distinct(peak, first_open)

epi_first_expr_unique <- tri_Model_result_epi_1Mb %>%
  distinct(gene, first_expr)

stroma_first_expr_unique <- tri_Model_result_stroma_1Mb %>%
  distinct(gene, first_expr)



library(dplyr)
library(ggplot2)
library(patchwork)

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

p_first_open <- ggplot(df_first_open, aes(x = pos, y = count, fill = cell_type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_open",
    x = "first_open",
    y = "Count",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14)





p_first_expr <- ggplot(df_first_expr, aes(x = pos, y = count, fill = cell_type)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  labs(
    title = "Distribution of first_expr",
    x = "first_expr",
    y = "Count",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 14)

p_first_open
p_first_expr


ggsave(
  filename = "Final_Plot/Fig_3/Fig_3F_pre_active_first_open.pdf",
  plot = p_first_open,
  width = 6,
  height = 3,
  units = "in"
)

ggsave(
  filename = "Final_Plot/Fig_3/Fig_3F_pre_active_first_expr.pdf",
  plot = p_first_expr,
  width = 6,
  height = 3,
  units = "in"
)







####Script_76_immune_gene_with_pre_open####

## Select immune genes expressed in immune-like stromal cells
Idents(pro_and_sec_seurat) <- 'stage_period_2_exclusive_cell_type_4'
immune_stroma_marker_gene <- FindMarkers(object = pro_and_sec_seurat, ident.1 = 'sec_late_immune_stroma', ident.2 = 'sec_late_surface_stroma', group.by = 'stage_period_2_exclusive_cell_type_4', min.pct = 0.1)
immune_marker_gene <- FindMarkers(object = pro_and_sec_seurat, ident.1 = 'sec_late_lymphoid', ident.2 = 'sec_late_surface_stroma', group.by = 'stage_period_2_exclusive_cell_type_4', min.pct = 0.1)

immune_stroma_marker_gene <- rownames(immune_stroma_marker_gene)[immune_stroma_marker_gene$avg_log2FC > 0.7 & immune_stroma_marker_gene$p_val_adj < 0.01]
immune_marker_gene <- rownames(immune_marker_gene)[immune_marker_gene$avg_log2FC > 0.7 & immune_marker_gene$p_val_adj < 0.01]

immune_marker_gene_in_immune_stroma <- intersect(immune_stroma_marker_gene, immune_marker_gene)


## Count how many of these genes show pre-opening
table(immune_stroma_marker_gene %in% unique(tri_Model_result_stroma_1Mb$gene))
table(immune_stroma_marker_gene %in% unique(tri_Model_result_stroma_1Mb$gene))

pre_active_immune_stroma_marker_gene <- immune_stroma_marker_gene[immune_stroma_marker_gene %in% unique(tri_Model_result_stroma_1Mb$gene)]
tri_Model_result_stroma_1Mb$distance_to_TSS <- tri_Model_result_stroma_1Mb$distance
peak_to_show <- tri_Model_result_stroma_1Mb[tri_Model_result_stroma_1Mb$gene %in% pre_active_immune_stroma_marker_gene & tri_Model_result_stroma_1Mb$distance_to_TSS < 55000,]



# gene_to_plot <- c( "GNLY", "LCP2", "IL2RB",
#                    "SEMA4D", "NCAM1", "ITGAX", "GNPTAB",
#                    "PECAM1", "SRGAP3", "ARHGEF6", "UBASH3B",
#                    "AFAP1L2", "ITGAD", "LNPEP", "LINC00299", "PTK2B", "RIN3",
#                    "ANKRD44", "CD38", "SLA", "MBP", "RUNX3", "GAS7", "CCDC88C")
# 
gene_to_plot <- c( "GNLY", "LCP2", "IL2RB", "SEMA4D", "NCAM1", "ITGAX", "UBASH3B", "ITGAD", "CD38", "SLA", "RUNX3")
gene_to_plot <- c( "ITGAD", "LCP2", "IL2RB", "SEMA4D", "NCAM1")



sec_stroma_cols <- c("sec_PO3_stroma","sec_PO5_stroma","sec_PO7_stroma","sec_PO9_stroma","sec_PO11_stroma","sec_PO13_stroma")
sec_epi_cols <- c("sec_PO3_epithelial","sec_PO5_epithelial","sec_PO7_epithelial","sec_PO9_epithelial","sec_PO11_epithelial","sec_PO13_epithelial")




immune_gene_pre_open_peak_to_plot <- c('chr16:31386141-31386641', 'chr5:170307476-170307976', 'chr22:37162447-37162947', 'chr9:89502034-89502534', 'chr11:112950266-112950766')
gene_to_plot <- c( "ITGAD", "LCP2", "IL2RB", "SEMA4D", "NCAM1")


for(i in seq(1,5)){
  peak_to_use <- immune_gene_pre_open_peak_to_plot[i]
  
  pdf(paste0('Final_Plot/Fig_3/Fig_3G_', gene_to_plot[i], '_pre_open.pdf'), height = 10, width = 4)
  
  extend_range <- 2000
  peak_part <- split_peak_name(peak_to_use)
  gene_ranges <- GRanges(seqnames = peak_part$chrom, ranges = IRanges(start = peak_part$start - floor(extend_range*1.5), end = peak_part$end + floor(extend_range*1.5)))
  p <- plotBrowserTrack(
    ArchRProj = new_whole_sec_to_use_with_RIF_ArhcR_Stroma_and_Gland_Filtered,
    groupBy = "period_super_cell_type",
    useGroups = sec_stroma_cols,
    region = gene_ranges,
    maxCells = 10000
  )
  .draw_any(p)
  
  dev.off()
  
}



## Plot expression of the selected genes

avg_list_binary[gene_to_plot, sec_stroma_cols]

plot_df <- avg_list_binary[gene_to_plot, sec_stroma_cols] |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  pivot_longer(
    cols = -gene,
    names_to = "stage",
    values_to = "expression"
  ) |>
  mutate(
    # Remove the "_stroma" suffix from column names
    stage = sub("_stroma$", "", stage),
    
    # Set the x-axis order
    stage = factor(
      stage,
      levels = c(
        "sec_PO3", "sec_PO5", "sec_PO7",
        "sec_PO9", "sec_PO11", "sec_PO13"
      )
    )
  )


p <- ggplot(
  plot_df,
  aes(
    x = stage,
    y = expression,
    color = gene,
    group = gene
  )
) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 2.5) +
  labs(
    x = "Stage",
    y = "Expression",
    color = "Gene"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p

ggsave(
  filename = "Final_Plot/Fig_3/Fig_3G_immune_gene_expr.pdf",
  plot = p,
  width = 6,
  height = 2,
  units = "in"
)









####Script_77_diff_of_gene_and_peak_between_Sync_and_pre-open####
library(VennDiagram)
library(grid)

Sync_peak_gene_stroma_High_in_sec_late <- Sync_peak_gene_stroma[Sync_peak_gene_stroma$gene %in% High_in_sec_late_gene_of_stroma, ]


A <- unique(tri_Model_result_stroma_1Mb$peak)  # Count unique values. Remove unique() to count all entries.
B <- unique(Sync_peak_gene_stroma_High_in_sec_late$peak)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))


pdf(
  file = "Final_Plot/Fig_3/Fig_3H_peak_Venn_of_diff_between_pre_and_Sync.pdf",
  width = 4,
  height = 4,
  useDingbats = FALSE
)

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)
dev.off()




pdf(
  file = "Final_Plot/Fig_3/Fig_3H_Gene_Venn_of_diff_between_pre_and_Sync.pdf",
  width = 4,
  height = 4,
  useDingbats = FALSE
)

A <- unique(tri_Model_result_stroma_1Mb$gene)  # Count unique values. Remove unique() to count all entries.
B <- unique(Sync_peak_gene_stroma_High_in_sec_late$gene)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#F0E68C", "#8A2BE2"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)

dev.off()





####Script_78_diff_of_gene_and_peak_between_epi_and_stroma_in_Sync-open####

pdf('Final_Plot/Fig_S3/Fig_S3A_peak_number_diff_of_Sync_open.pdf')
A <- unique(Sync_of_stroma$peak)  # Count unique values. Remove unique() to count all entries.
B <- unique(Sync_of_epi$peak)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)
dev.off()


## Gene summary
pdf('Final_Plot/Fig_S3/Fig_S3A_gene_number_diff_of_Sync_open.pdf')
A <- unique(Sync_of_stroma$gene)  # Count unique values. Remove unique() to count all entries.
B <- unique(Sync_of_epi$gene)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)

dev.off()






####Script_79_diff_of_gene_and_peak_between_epi_and_stroma_in_pre-open####


pdf('Final_Plot/Fig_S3/Fig_S3C_peak_diff_between_epi_and_stroma_of_pre_open.pdf')
A <- unique(tri_Model_result_stroma_1Mb$peak)  # Count unique values. Remove unique() to count all entries.
B <- unique(tri_Model_result_epi_1Mb$peak)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)
dev.off()


## Gene summary
pdf('Final_Plot/Fig_S3/Fig_S3C_gene_diff_between_epi_and_stroma_of_pre_open.pdf')
A <- unique(tri_Model_result_stroma_1Mb$gene)  # Count unique values. Remove unique() to count all entries.
B <- unique(tri_Model_result_epi_1Mb$gene)

nA  <- length(A)
nB  <- length(B)
nAB <- length(intersect(A, B))

grid.newpage()
draw.pairwise.venn(
  area1 = nA,
  area2 = nB,
  cross.area = nAB,
  category = c("A", "B"),
  scaled = TRUE,               # Scale circle sizes by value
  fill = c("#4DBBD5", "#E64B35"),
  alpha = 0.5,
  cat.cex = 1.2,
  cex = 1.2
)

dev.off()









####Script_710_pre-open_peak_gene_link_in_stroma_and_gland####

aaa <- dim(tri_Model_result_stroma_1Mb)[1]
bbb <- dim(tri_Model_result_epi_1Mb)[1]


# Build the plotting data
plot_data <- data.frame(
  cell_type = factor(
    c("stroma", "epi"),
    levels = c("stroma", "epi")
  ),
  number = c(aaa, bbb)
)

# Draw a bar chart

p <- ggplot(plot_data, aes(x = cell_type, y = number, fill = cell_type)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = number),
    vjust = -0.4,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "stroma" = "#4DBBD5",
      "epi" = "#E64B35"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Number",
    title = "number of pre-active peak-gene link"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0.5,
    ),
    axis.text.x = element_text(
      size = 16,
    )
  )


ggsave(
  filename = "Final_Plot/Fig_S3/Fig_S3C_pre_open_paek_gene_link_number_of_epi_and_stroma.pdf",
  plot = p,
  width = 3.2,
  height = 4,
  units = "in"
)




peak_number_distribution_of_stroma <- as.integer(table(tri_Model_result_stroma_1Mb$gene))
peak_number_distribution_of_epi <- as.integer(table(tri_Model_result_epi_1Mb$gene))

library(ggplot2)

count_df <- data.frame(
  category = rep(c("<3", ">=3"), each = 2),
  sample = rep(c("stroma", "epi"), times = 2),
  count = c(
    sum(peak_number_distribution_of_stroma < 3, na.rm = TRUE),
    sum(peak_number_distribution_of_epi < 3, na.rm = TRUE),
    sum(peak_number_distribution_of_stroma >= 3, na.rm = TRUE),
    sum(peak_number_distribution_of_epi >= 3, na.rm = TRUE)
  )
)

count_df$x_pos <- ifelse(count_df$category == "<3", 2.2, 3)
count_df$sample <- factor(
  count_df$sample,
  levels = c("stroma", "epi")
)
p <- ggplot(count_df, aes(x = x_pos, y = count, fill = sample)) +
  geom_col(
    position = position_dodge(width = 0.55),
    width = 0.55
  ) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.55),
    vjust = -0.3
  ) +
  scale_fill_manual(
    values = c(
      "stroma" = "#4DBBD5",
      "epi" = "#E64B35"
    )
  ) +
  scale_x_continuous(
    breaks = c(2.2, 3),
    labels = c("<3", ">=3")
  ) +
  theme_classic() +
  labs(
    x = "Peak number category",
    y = "Number of elements",
    fill = "Sample",
    title = "Peak number distribution"
  )


ggsave(
  filename = "Final_Plot/Fig_S3/Fig_S3D_pre_open_gene_number_if_more_than_3.pdf",
  plot = p,
  width = 4,
  height = 4,
  units = "in"
)






####Script_711_GO_enrichment_analysis_of_pre-opening_stromal_genes####
draw_KEGG_picture(gene_set = unique(tri_Model_result_stroma_1Mb$gene), top_n = 13)
draw_KEGG_picture(gene_set = unique(tri_Model_result_stroma_1Mb$gene), top_n = 20, y_text_size = 10, n_to_plot = c(1, 2, 3, 4, 6), color_by_p.adjust = FALSE, single_color = 'red')
aaa <- draw_GO_picture(gene_set = unique(tri_Model_result_stroma_1Mb$gene), top_n = 50, return_dataframe = TRUE)

bbb <- draw_GO_picture(gene_set = unique(tri_Model_result_stroma_1Mb$gene), 
                       select_description = c('extracellular matrix organization', 'cell-matrix adhesion', 'focal adhesion assembly', 'actin filament organization','leukocyte migration'))


ggsave(
  filename = "Final_Plot/Fig_S3/Fig_S3E_pre_open_gene_GO_enrichment.pdf",
  plot = bbb$padj_plot,
  width = 9,
  height = 3,
  units = "in"
)






####Script_712_number_of_immune_gene_with_pre-opening####
library(ggplot2)

# Count TRUE and FALSE values
link_count <- table(
  immune_stroma_marker_gene %in% unique(tri_Model_result_stroma_1Mb$gene)
)

# Build the plotting data
plot_data <- data.frame(
  category = factor(
    c("Have pre-gene link", "Don't have pre-gene link"),
    levels = c("Have pre-gene link", "Don't have pre-gene link")
  ),
  count = c(
    unname(link_count["TRUE"]),
    unname(link_count["FALSE"])
  )
)

# Draw a bar chart
p <- ggplot(plot_data, aes(x = category, y = count)) +
  geom_col(
    width = 0.65,
    fill = "#4C78A8"
  ) +
  geom_text(
    aes(label = count),
    vjust = -0.4,
    size = 4
  ) +
  labs(
    x = NULL,
    y = "Number of marker genes"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  theme_classic(base_size = 13)


ggsave(
  filename = "Final_Plot/Fig_S3/Fig_S3F_number_have_pre_open_link_of_immune_gene.pdf",
  plot = p,
  width = 3,
  height = 3,
  units = "in"
)






