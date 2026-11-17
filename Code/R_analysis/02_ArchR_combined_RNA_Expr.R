


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





####01_get_ArchR_project####

read_data_from_files <- function(data_name) {
  primary_file <- get_file_path(data_name)
  if (file.exists(primary_file)){
    file_path <- primary_file
  }else{
    secondary_file <- secondary_path(data_name)
    file_path <- secondary_file
  }
  if (file.exists(file_path)) {
    # Choose a reader for the file format, such as read.table or read.csv
    count.m <- read.table(file_path)  # Read the TSV file
    idx <- grep("ENSG", colnames(count.m))
    count.m <- count.m[,setdiff(1:length(colnames(count.m)), idx)]      # Exclude genes with ENSG identifiers
    count.m <- as.matrix(count.m)
    colnames(count.m) <- gsub("_", "-",colnames(count.m))
    tmp <- CreateSeuratObject(counts= t(count.m), project= data_name, assay ="RNA")
    tmp$channel = c("10um")
    tmp$stage =c("implantation")
    tmp$coord = paste0( data_name,"_",row.names(tmp@meta.data))
    tmp$nA = as.integer(str_split_fixed(tmp$coord,"_|x",7)[,6])
    tmp$nB = as.integer(str_split_fixed(tmp$coord,"_|x",7)[,7])
    print(paste0("已完成:", data_name))
    return(tmp)
  } else {
    warning(paste("文件不存在：", file_path))
    return(NULL)
  }
}


### Map barcoded cell names to coordinates and return a data frame #######
get_change_name <- function(ArchR_obj, ID_file='/mnt/server2/data5/Spatial_endometrium/core/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt'){
  cell_names <- ArchR_obj$cellNames
  # Read ID_file
  id_file <- read.table(ID_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("Barcode", "nA", "nB"))
  # Extract barcodes and sample names from the cell names
  cell_barcodes <- sub(".*#", "", cell_names)
  cell_orig_ident <- sub("#.*", "", cell_names)
  # Build a data frame for matching
  cell_data <- data.frame(CellName = cell_names, orig.ident = cell_orig_ident, Barcode = cell_barcodes, stringsAsFactors = FALSE)
  # Join cell names to ID_file
  merged_data <- cell_data %>%
    left_join(id_file, by = c("Barcode" = "Barcode"))
  # Replace each barcode with its matched coordinates
  merged_data <- merged_data %>%
    mutate(NewCellName = ifelse(!is.na(nA) & !is.na(nB),
                                paste0(gsub("#.*", "", CellName), "#", nA, "_", nB),
                                CellName))  # Keep the original name if no coordinates match
  return(merged_data)
}



getBarcodeFromCellName <- function(cell_name, id_file="/mnt/server2/data5/Spatial_endometrium/core/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt") {
  # Extract x and y coordinates
  ID_file <- read.table(id_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("Barcode", "nA", "nB"))
  coords <- strsplit(cell_name, "x")[[1]]
  x_coord <- as.integer(coords[1])  # X coordinate on the left
  y_coord <- as.integer(coords[2])  # Y coordinate on the right
  
  # Look up the matching barcode
  barcode <- ID_file$Barcode[ID_file$nA == x_coord & ID_file$nB == y_coord]
  
  # Return the matching barcode
  return(barcode)
}


get_CellName_From_Barcode <- function(cell_name, id_file="/mnt/server2/data5/Spatial_endometrium/core/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt") {
  # Extract x and y coordinates
  ID_file <- read.table(id_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = c("Barcode", "nA", "nB"))
  coords <- strsplit(cell_name, "x")[[1]]
  x_coord <- as.integer(coords[1])  # X coordinate on the left
  y_coord <- as.integer(coords[2])  # Y coordinate on the right
  
  # Look up the matching barcode
  barcode <- ID_file$Barcode[ID_file$nA == x_coord & ID_file$nB == y_coord]
  
  # Return the matching barcode
  return(barcode)
}

### Transfer matching Seurat metadata to ArchR
transfer_Seurat_metadata_to_ArchR_cellColData <- function(ArchR_obj, seurat_obj){
  seurat_meta <- seurat_obj@meta.data
  seurat_meta$Universal_name <- as.character(seurat_meta$Universal_name)
  archr_cells <- getCellColData(projHeme2)
  archr_cells$cellName <- rownames(archr_cells)
  archr_cells$Universal_name <- projHeme2$Universal_name
  meta_merge <- merge(archr_cells, seurat_meta, by = "Universal_name", all.x = TRUE)
  
}



#### Convert the TSV output to ArchR objects ####

## Core settings
wd='/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/'
setwd(wd)
file_list_ATAC_name=paste0('/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/File_list_ATAC_all_to_use.txt')
file_ATAC_names <- readLines(file_list_ATAC_name)
file_ATAC_names <- trimws(file_ATAC_names)
file_ATAC_names <- file_ATAC_names[file_ATAC_names != ""]
file_ATAC_list <- read.table(file_list_ATAC_name)

file_list_RNA_name=paste0('/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/File_list_ATAC_paired_RNA.txt')
file__RNA_names <- readLines(file_list_RNA_name)
file_RNA_names <- trimws(file__RNA_names)
file_RNA_names <- file_RNA_names[file_RNA_names != ""]
file_RNA_list <- read.table(file_list_RNA_name)


### To run samples individually, update the parameters below and run the function body.
ATAC_name <- file_ATAC_list[[1]][6]
RNA_name <- file_RNA_list[[1]][6]
id_file <- id_file_2
save_dir = 'ATAC/ArchR_project_with_RIF_True/'
seurat_obj <- whole_sec_seurat_with_RIF


### Create ArchR Arrow files
addArchRThreads(threads = 10)
addArchRGenome("hg38")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")          ### Check whether omitting the three "w" characters causes an error.



#### Seurat ####
pro_and_sec_seurat_for_ArchR <- get_neighbors_and_clusters(pro_and_sec_seurat)
pro_and_sec_seurat_for_ArchR <- subset(pro_and_sec_seurat, subset = stage != 'RIF_PO7')

pro_and_sec_seurat_for_ArchR <- RunHarmony(pro_and_sec_seurat_for_ArchR, group.by.vars = c('exper_method', 'orig.ident'))
pro_and_sec_seurat_for_ArchR <- RunUMAP(pro_and_sec_seurat_for_ArchR, reduction = "harmony", dims = 1:20, reduction.name = 'harmony_umap', reduction.key = 'harmony_umap')
pro_and_sec_seurat_for_ArchR <- FindNeighbors(pro_and_sec_seurat_for_ArchR, reduction = 'harmony', dims = 1:30, graph.name = c('harmony_orig_nn', 'harmony_orig_snn'))
pro_and_sec_seurat_for_ArchR <- FindClusters(pro_and_sec_seurat_for_ArchR, resolution = seq(0.3, 0.7, 0.1), graph.name = 'harmony_orig_snn')



## Plot known markers using celltype.ls defined above
pdf(paste0('Cell_type_list_seurat.pdf'))
for (i in celltype.ls){
  p <- my_Dotplot_plot(seurat_obj, marker_gene = unique(get(i)), title = i, cluster = 'harmony_orig_snn_res.0.5')
  print(p)
}
dev.off()














#### ArchR ####

## Set local paths, convert the ATAC pipeline TSV.GZ files to Arrow files, and create an ArchR project.

##eg
ArrowFiles_pro_late_E19 <- createArrowFiles(
  inputFiles = 'pro_late_E19_ATAC_Add.sorted.barcoded.rmdup.fragment.tsv.gz',
  sampleNames = 'pro_late_E19_240718',
  # outputNames = paste0('/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/ATAC/ArrowFiles/', ATAC_name),
  minTSS = 1, # Don't set this too high because you can always increase later
  minFrags = 500, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

pro_and_sec_with_RIF_ArchR <- ArchRProject(
  ArrowFiles = all_ArrowFiles,
  outputDirectory = "ArchR_pro_and_sec_True",
  copyArrows = TRUE
)



#### Transfer Seurat annotations and the gene expression matrix
id_file_1 <- 'combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt'
id_file_2 <- 'combine_barcode.round2round1_index1_index2.v7_big_r1_25to97_r2_9to98.20241209.txt.IDFILE'

addSpatialCoordToArchR <- function(ArchRProj, sample_name, ID_file) {
  
  ## 1. Check inputs
  if (!inherits(ArchRProj, "ArchRProject")) {
    stop("ArchRProj 必须是一个 ArchRProject 对象。")
  }
  
  if (!sample_name %in% ArchRProj$Sample) {
    stop("sample_name 不在 ArchRProj$Sample 中：", sample_name)
  }
  
  ## 2. Read ID_file
  if (is.character(ID_file) && length(ID_file) == 1 && file.exists(ID_file)) {
    id_df <- read.table(
      ID_file,
      header = FALSE,
      sep = "",
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  } else if (is.data.frame(ID_file)) {
    id_df <- ID_file
  } else {
    stop("ID_file 必须是一个存在的文件路径，或者一个 data.frame。")
  }
  
  if (ncol(id_df) < 3) {
    stop("ID_file 至少需要 3 列：barcode, nA, nB")
  }
  
  id_df <- id_df[, 1:3, drop = FALSE]
  colnames(id_df) <- c("barcode", "nA", "nB")
  
  ## Handle an optional header row
  first_row <- tolower(trimws(as.character(unlist(id_df[1, ]))))
  if (identical(first_row, c("barcode", "na", "nb"))) {
    id_df <- id_df[-1, , drop = FALSE]
  }
  
  ## 3. Format ID_file
  id_df$barcode <- trimws(as.character(id_df$barcode))
  id_df$nA <- suppressWarnings(as.numeric(as.character(id_df$nA)))
  id_df$nB <- suppressWarnings(as.numeric(as.character(id_df$nB)))
  
  if (anyDuplicated(id_df$barcode)) {
    dup_bc <- unique(id_df$barcode[duplicated(id_df$barcode)])
    stop(
      "ID_file 中 barcode 存在重复，无法一一对应。重复 barcode 示例：",
      paste(head(dup_bc, 10), collapse = ", ")
    )
  }
  
  ## 4. Select cells from the specified sample
  all_cells <- ArchRProj$cellNames
  all_samples <- ArchRProj$Sample
  
  sample_idx <- which(all_samples == sample_name)
  sample_cells <- all_cells[sample_idx]
  
  ## 5. Extract the barcode after "#" in cellNames
  sample_barcodes <- sub(".*#", "", sample_cells)
  
  ## 5. Match barcodes
  match_idx <- match(sample_barcodes, id_df$barcode)
  
  matched_nA <- id_df$nA[match_idx]
  matched_nB <- id_df$nB[match_idx]
  
  ## 6. Preserve columns populated for other samples
  meta_cols <- colnames(ArchRProj@cellColData)
  
  if ("nA" %in% meta_cols) {
    nA_vec <- as.numeric(ArchRProj@cellColData[["nA"]])
  } else {
    nA_vec <- rep(NA_real_, length(all_cells))
  }
  
  if ("nB" %in% meta_cols) {
    nB_vec <- as.numeric(ArchRProj@cellColData[["nB"]])
  } else {
    nB_vec <- rep(NA_real_, length(all_cells))
  }
  
  if ("spatial_barcode" %in% meta_cols) {
    barcode_vec <- as.character(ArchRProj@cellColData[["spatial_barcode"]])
  } else {
    barcode_vec <- rep(NA_character_, length(all_cells))
  }
  
  if ("spatial_coord" %in% meta_cols) {
    coord_vec <- as.character(ArchRProj@cellColData[["spatial_coord"]])
  } else {
    coord_vec <- rep(NA_character_, length(all_cells))
  }
  
  ## 7. Fill in the results for the current sample
  nA_vec[sample_idx] <- matched_nA
  nB_vec[sample_idx] <- matched_nB
  barcode_vec[sample_idx] <- sample_barcodes
  coord_vec[sample_idx] <- ifelse(
    is.na(matched_nA) | is.na(matched_nB),
    NA,
    paste0(matched_nA, "_", matched_nB)
  )
  
  ## 8. Write back to the ArchR project
  ArchRProj$nA <- nA_vec
  ArchRProj$nB <- nB_vec
  ArchRProj$spatial_barcode <- barcode_vec
  ArchRProj$spatial_coord <- coord_vec
  
  ## 9. Report matching results
  n_matched <- sum(!is.na(match_idx))
  n_total <- length(sample_cells)
  
  message("Sample: ", sample_name)
  message("该 sample 细胞数: ", n_total)
  message("成功匹配 barcode 数: ", n_matched)
  message("未匹配 barcode 数: ", n_total - n_matched)
  
  return(ArchRProj)
}



ATAC_sample_ID_file <- read.table(file = 'ATAC_Sample_ID_file.txt')
for(i in seq(dim(ATAC_sample_ID_file)[1])){
  pro_and_sec_with_RIF_ArchR <- addSpatialCoordToArchR(
    ArchRProj   = pro_and_sec_with_RIF_ArchR,
    sample_name = ATAC_sample_ID_file$V1[i],
    ID_file     = paste0('/mnt/server1/data9/zouqb/ID_file/', ATAC_sample_ID_file$V2[i])
  )
}



unique(pro_and_sec_with_RIF_ArchR$Sample)

pro_and_sec_with_RIF_ArchR$sample_id <- sapply(pro_and_sec_with_RIF_ArchR$Sample, function(x){
  aaa <- unlist(str_split(x, pattern = '_'))[1:3]
  paste0(aaa, collapse = '_')
})

pro_and_sec_with_RIF_ArchR$Universal_name <- paste0(pro_and_sec_with_RIF_ArchR$sample_id, '#', pro_and_sec_with_RIF_ArchR$nA, 'x', pro_and_sec_with_RIF_ArchR$nB)

pro_and_sec_seurat_for_ArchR$Universal_name <- paste0(pro_and_sec_seurat_for_ArchR$sample_id, '#', pro_and_sec_seurat_for_ArchR$nA, 'x', pro_and_sec_seurat_for_ArchR$nB)



### Add cluster metadata and keep shared cells before integrating RNA data.
seurat_meta <- pro_and_sec_seurat_for_ArchR@meta.data
seurat_meta$Universal_name <- as.character(seurat_meta$Universal_name)

archr_meta <- as.data.frame(getCellColData(pro_and_sec_with_RIF_ArchR))
archr_meta$cellName <- rownames(archr_meta)
archr_meta$Universal_name <- as.character(pro_and_sec_with_RIF_ArchR$Universal_name)


# Keep Universal_name values shared by both objects
merged_meta <- merge(
  x = archr_meta,
  y = seurat_meta,
  by = "Universal_name",
  all = FALSE  # inner join
)
# 3. Find shared Universal_name values
shared_universal_names <- unique(merged_meta$Universal_name)

# Subset Seurat by cell names in rownames(meta.data)
cells_to_keep_seurat <- rownames(seurat_meta)[seurat_meta$Universal_name %in% shared_universal_names]
pro_and_sec_seurat_for_ArchR <- subset(pro_and_sec_seurat_for_ArchR, cells = cells_to_keep_seurat)

# Subset ArchR by pro_and_sec_with_RIF_ArchR$Universal_name
cells_to_keep_archr <- pro_and_sec_with_RIF_ArchR$cellNames[pro_and_sec_with_RIF_ArchR$Universal_name %in% shared_universal_names]
pro_and_sec_with_RIF_ArchR <- subsetArchRProject(
  ArchRProj = pro_and_sec_with_RIF_ArchR,
  cells = cells_to_keep_archr,
  outputDirectory = '/mnt/server1/data9/zouqb/ttttttry/ATAC_RNA/multi_data/ATAC/ArchR_Final_tmp/', # Set the output path
  dropCells = TRUE,
  force = T
)


## Select metadata fields to add
archr_fields <- colnames(getCellColData(pro_and_sec_with_RIF_ArchR))
fields_to_add <- setdiff(colnames(seurat_meta), c("Universal_name", archr_fields))
print(fields_to_add)
print(merged_meta)

for (field in fields_to_add) {
  print(field)
  vec <- merged_meta[[field]]
  if (is.factor(vec)) {
    vec <- as.character(vec)
  }
  pro_and_sec_with_RIF_ArchR <- addCellColData(
    ArchRProj = pro_and_sec_with_RIF_ArchR,
    data = vec,
    cells = merged_meta$cellName,
    name = field
  )
}


print(my_spatial_plot_ATAC(ArchR_obj = pro_and_sec_with_RIF_ArchR, cluster = 'exclusive_cell_type', sample_to_use = 'RIF_PO7_E10_241223'))



## Add the GeneExpression matrix
DefaultAssay(pro_and_sec_seurat_for_ArchR) <- "RNA"

pro_and_sec_seurat_for_ArchR <- JoinLayers(pro_and_sec_seurat_for_ArchR)

rna_counts <- Seurat::GetAssayData(
  object = pro_and_sec_seurat_for_ArchR,
  assay  = "RNA",
  slot   = "counts"
)

rna_counts <- Matrix::Matrix(rna_counts, sparse = TRUE)

message("Seurat RNA counts dim: ",
        paste(dim(rna_counts), collapse = " x "))


# -----------------------------
# 2) Map Seurat cell names to ArchR cellNames
#    Reuse the existing Universal_name
# -----------------------------

seurat_meta2 <- pro_and_sec_seurat_for_ArchR@meta.data
seurat_meta2$seurat_cell <- rownames(seurat_meta2)
seurat_meta2$Universal_name <- as.character(seurat_meta2$Universal_name)

archr_meta2 <- as.data.frame(
  getCellColData(pro_and_sec_with_RIF_ArchR),
  stringsAsFactors = FALSE
)
archr_meta2$archr_cell <- rownames(archr_meta2)
archr_meta2$Universal_name <- as.character(archr_meta2$Universal_name)

# Check inputs
if (!"Universal_name" %in% colnames(seurat_meta2)) {
  stop("Seurat meta.data 中没有 Universal_name，请先运行前面的 Universal_name 构建代码。")
}

if (!"Universal_name" %in% colnames(archr_meta2)) {
  stop("ArchR cellColData 中没有 Universal_name，请先运行前面的 Universal_name 构建代码。")
}

if (any(is.na(seurat_meta2$Universal_name))) {
  stop("Seurat 中存在 NA 的 Universal_name，请检查 sample_id/nA/nB。")
}

if (any(is.na(archr_meta2$Universal_name))) {
  stop("ArchR 中存在 NA 的 Universal_name，请检查 sample_id/nA/nB。")
}

# Universal_name must be unique for one-to-one matching
if (anyDuplicated(seurat_meta2$Universal_name)) {
  dup_use <- unique(seurat_meta2$Universal_name[duplicated(seurat_meta2$Universal_name)])
  stop(
    "Seurat 中 Universal_name 存在重复，无法唯一映射。示例：",
    paste(head(dup_use, 10), collapse = ", ")
  )
}

if (anyDuplicated(archr_meta2$Universal_name)) {
  dup_use <- unique(archr_meta2$Universal_name[duplicated(archr_meta2$Universal_name)])
  stop(
    "ArchR 中 Universal_name 存在重复，无法唯一映射。示例：",
    paste(head(dup_use, 10), collapse = ", ")
  )
}

cell_map <- merge(
  x = seurat_meta2[, c("Universal_name", "seurat_cell"), drop = FALSE],
  y = archr_meta2[, c("Universal_name", "archr_cell"), drop = FALSE],
  by = "Universal_name",
  all = FALSE
)

message("Seurat-ArchR matched cells: ", nrow(cell_map))

if (nrow(cell_map) == 0) {
  stop("Seurat 和 ArchR 没有任何细胞通过 Universal_name 匹配上。请检查 sample_id/nA/nB 格式。")
}


# -----------------------------
# 3) Rename RNA count columns with the matched ArchR cellNames
# -----------------------------

archr_cell_for_seurat <- cell_map$archr_cell[
  match(colnames(rna_counts), cell_map$seurat_cell)
]

keep_cols <- !is.na(archr_cell_for_seurat)

rna_counts2 <- rna_counts[, keep_cols, drop = FALSE]
colnames(rna_counts2) <- archr_cell_for_seurat[keep_cols]

message("RNA counts after cell matching dim: ",
        paste(dim(rna_counts2), collapse = " x "))

# For duplicate columns, keep the one with the highest total count
if (any(duplicated(colnames(rna_counts2)))) {
  message("检测到重复的 ArchR cellName 映射列名，将按列总和保留最大者。")
  
  cs <- Matrix::colSums(rna_counts2)
  o <- order(colnames(rna_counts2), -cs)
  
  rna_counts2 <- rna_counts2[, o, drop = FALSE]
  rna_counts2 <- rna_counts2[, !duplicated(colnames(rna_counts2)), drop = FALSE]
}

# With strictMatch = TRUE, seRNA column names must exactly match ArchR cellNames
archr_cells_use <- pro_and_sec_with_RIF_ArchR$cellNames

missing_in_rna <- setdiff(archr_cells_use, colnames(rna_counts2))
extra_in_rna <- setdiff(colnames(rna_counts2), archr_cells_use)

if (length(extra_in_rna) > 0) {
  message("RNA counts 中存在不属于当前 ArchR project 的细胞，将移除：", length(extra_in_rna))
  rna_counts2 <- rna_counts2[, setdiff(colnames(rna_counts2), extra_in_rna), drop = FALSE]
}

if (length(missing_in_rna) > 0) {
  stop(
    "有 ArchR 细胞没有对应的 RNA counts，strictMatch=TRUE 会失败。缺失细胞数：",
    length(missing_in_rna),
    "\n示例：",
    paste(head(missing_in_rna, 10), collapse = ", ")
  )
}

# Reorder to match ArchR cellNames
rna_counts2 <- rna_counts2[, archr_cells_use, drop = FALSE]

stopifnot(identical(colnames(rna_counts2), archr_cells_use))


# -----------------------------
# 4) Match Seurat gene symbols to ArchR gene annotations
# -----------------------------

genes_archr <- getGenes(pro_and_sec_with_RIF_ArchR)

# Allow gene symbol column names to vary across ArchR versions
gene_mcols <- as.data.frame(mcols(genes_archr))

if ("symbol" %in% colnames(gene_mcols)) {
  gene_symbol_col <- "symbol"
} else if ("name" %in% colnames(gene_mcols)) {
  gene_symbol_col <- "name"
} else {
  stop("getGenes() 返回的 GRanges 中找不到 symbol 或 name 列，请检查 gene annotation。")
}

archr_symbols <- as.character(gene_mcols[[gene_symbol_col]])
se_genes_raw <- rownames(rna_counts2)

idx <- match(se_genes_raw, archr_symbols)

# Fallback: Seurat gene names may use "." where ArchR uses "-"
need <- which(is.na(idx))
if (length(need) > 0) {
  idx2 <- match(gsub("\\.", "-", se_genes_raw[need]), archr_symbols)
  idx[need] <- idx2
}

keep_genes <- which(!is.na(idx))

if (length(keep_genes) == 0) {
  stop("Seurat 基因与 ArchR gene annotation 完全匹配不上，请检查 gene symbol 体系。")
}

message("Matched genes: ", length(keep_genes), " / ", length(se_genes_raw))


# -----------------------------
# 5) Resolve multiple Seurat genes mapping to one ArchR gene
#    This typically occurs when "." is replaced with "-"
#    Keep the entry with the highest total count
# -----------------------------

dup_archr_idx <- idx[keep_genes]
tab <- table(dup_archr_idx)
dup_vals <- as.integer(names(tab[tab > 1]))

if (length(dup_vals) > 0) {
  message("检测到多个 Seurat gene 映射到同一个 ArchR gene，数量：", length(dup_vals))
  message("将对每个重复 ArchR gene 保留 rowSums 最大的 Seurat gene。")
  
  rs <- Matrix::rowSums(rna_counts2[keep_genes, , drop = FALSE])
  keep_flag <- rep(TRUE, length(keep_genes))
  
  for (v in dup_vals) {
    pos <- which(dup_archr_idx == v)
    best <- pos[which.max(rs[pos])]
    rmv <- setdiff(pos, best)
    keep_flag[rmv] <- FALSE
  }
  
  keep_genes <- keep_genes[keep_flag]
}

rna_counts3 <- rna_counts2[keep_genes, , drop = FALSE]
genes_gr <- genes_archr[idx[keep_genes]]

# -----------------------------
# 6) Prepare ranges and names
#    Make row names unique, sync rowRanges names, and sort by genomic coordinates
# -----------------------------

gene_symbols_use <- archr_symbols[idx[keep_genes]]
gene_key <- make.unique(as.character(gene_symbols_use))

rownames(rna_counts3) <- gene_key
names(genes_gr) <- gene_key

# Write symbol/name fields so getGroupSE output can be matched to gene names
mcols(genes_gr)$symbol <- gene_key
mcols(genes_gr)$name <- gene_key

# Sort by genomic coordinates to keep rows aligned across write and read operations
ord <- order(
  as.character(seqnames(genes_gr)),
  start(genes_gr),
  end(genes_gr),
  gene_key
)

rna_counts3 <- rna_counts3[ord, , drop = FALSE]
genes_gr <- genes_gr[ord]

names(genes_gr) <- rownames(rna_counts3)

message("Final RNA matrix dim for ArchR GeneExpressionMatrix: ",
        paste(dim(rna_counts3), collapse = " x "))

stopifnot(identical(colnames(rna_counts3), pro_and_sec_with_RIF_ArchR$cellNames))
stopifnot(identical(rownames(rna_counts3), names(genes_gr)))


# -----------------------------
# 7) Build a SummarizedExperiment and write it to ArchR
# -----------------------------

seRNA <- SummarizedExperiment(
  assays = list(
    counts = Matrix::Matrix(rna_counts3, sparse = TRUE)
  ),
  rowRanges = genes_gr
)

rownames(seRNA) <- rownames(rna_counts3)
colnames(seRNA) <- colnames(rna_counts3)

stopifnot(identical(colnames(seRNA), pro_and_sec_with_RIF_ArchR$cellNames))
stopifnot(identical(rownames(seRNA), rownames(rna_counts3)))

pro_and_sec_with_RIF_ArchR <- addGeneExpressionMatrix(
  input       = pro_and_sec_with_RIF_ArchR,
  seRNA       = seRNA,
  strictMatch = TRUE,
  force       = TRUE,
  scaleTo     = 10000
)

# Check that the matrix was added
getAvailableMatrices(pro_and_sec_with_RIF_ArchR)






## Inspect GeneExpression data

library(ArchR)
library(SummarizedExperiment)
library(tibble)
library(tidyr)
library(ggplot2)

# Confirm the matrix name
getAvailableMatrices(pro_and_sec_with_RIF_ArchR)

# Genes to inspect
genes_use <- c("SOX2", "MKI67", "GATA1", "IGFBP1", 'XDH')

# Summarize GeneExpressionMatrix by stage
gse_stage <- getGroupSE(
  ArchRProj = pro_and_sec_with_RIF_ArchR,
  useMatrix = "GeneExpressionMatrix",
  groupBy = "stage",
  divideN = TRUE,
  scaleTo = 10000,
  verbose = FALSE
)

# Extract gene names
feature_names <- if ("name" %in% colnames(rowData(gse_stage))) {
  rowData(gse_stage)$name
} else if (!is.null(rownames(gse_stage))) {
  rownames(gse_stage)
} else {
  stop("找不到基因名，请先检查 rowData(gse_stage)")
}

# Match the genes of interest
idx <- match(genes_use, feature_names)
keep <- !is.na(idx)

mean_mat <- assay(gse_stage, assayNames(gse_stage)[1])[idx[keep], , drop = FALSE]
rownames(mean_mat) <- genes_use[keep]

# Reshape to long format
mean_df <- as.data.frame(as.matrix(mean_mat)) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "stage", values_to = "mean_expr")

# Plot stages in the order returned by ArchR
mean_df$stage <- factor(mean_df$stage, levels = colnames(gse_stage))

ggplot(mean_df, aes(stage, mean_expr, group = gene, color = gene)) +
  geom_point(size = 2) +
  geom_line() +
  theme_classic() +
  labs(x = "stage", y = "mean expression")


 