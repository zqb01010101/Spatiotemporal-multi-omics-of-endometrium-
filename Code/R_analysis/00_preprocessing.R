#### Create Seurat objects from ST Pipeline output


####Scirpt_01_get_seurat_from_STpineline####

#### Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(Matrix)
library(scales)
library(plotly)


## Read input files into a Seurat object
read_data_from_files <- function(data_name) {
  file_path <- get_file_path(data_name)
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



## Use read_data_from_files on the output from your ST Pipeline run.
## project_RNA and project_ATAC contain RNA and ATAC sample names, respectively.

## Store the Seurat objects in a list
file_data_list_RNA <- lapply(project_RNA, read_data_from_files)
names(file_data_list_RNA) <- project

file_data_list_ATAC <- lapply(project_ATAC, read_data_from_files)
names(file_data_list_ATAC) <- project




####Script_02_filtered_cell####

Min_nCounts=100
file_data_list_RNA <- lapply(file_data_list_RNA, function(seurat_obj) {
  subset(seurat_obj, subset = nCount_RNA > Min_nCounts)
})

file_data_list_ATAC <- lapply(file_data_list_ATAC, function(seurat_obj) {
  subset(seurat_obj, subset = nCount_RNA > Min_nCounts)
})

## Filter using the H&E images




####Script_03_get_merged_Seurat_and_cluster####
pro_and_sec_seurat <- merge(file_data_list_RNA[[1]], y = file_data_list_RNA[[-1]])
pro_and_sec_seurat_for_ArchR <- merge(file_data_list_ATAC[[1]], y = file_data_list_ATAC[[-1]])

