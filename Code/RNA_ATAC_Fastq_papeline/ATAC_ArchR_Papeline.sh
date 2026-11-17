

data_dir=$1
sample=$2

threads=10

cd ${data_dir}


cat > ${data_dir}/${sample}_ArchR.R << EOF
library(ArchR)
library(BSgenome)
addArchRThreads(threads = ${threads})
addArchRGenome("hg38")

args<-commandArgs()
ID_file = read.table( file = paste(args[6]) )
colnames(ID_file) = c('barcode' , 'A' , 'B')

sample='${sample}'
wd = '${data_dir}'
setwd(wd)

inputFiles = '../chromap_output/rmdup_chromap_${sample}.bam'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = '${sample}',
  bcTag = "CB" ,
  minFrags = 1 ,
  maxFrags = 9999999999999,
  offsetPlus = 4,
  offsetMinus = -5,
  addTileMat = F ,
  minTSS = 0 ,
  filterTSS = NULL ,
  #promoterRegion = c(-2000 , 2000) , ### cause error in creat arrow file
  addGeneScoreMat = F ,
  force = T ,
  logFile = 'ArchR_log.txt' ,
  threads = ${threads}
)

proj1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '${data_dir}' ,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)


library(ggplot2)
library(viridis)
library(ggpointdensity)

cell_df =  as.data.frame(proj1@cellColData)
total_fragment = sum(cell_df\$nFrags)
nCells_detected = nrow(cell_df)
median_frags = median(cell_df\$nFrags)
median_TSSe = median(cell_df\$TSSEnrichment)
mean_frags = mean(cell_df\$nFrags)
mean_TSSe = mean(cell_df\$TSSEnrichment)



sample='${sample}'
summary = as.data.frame(sample)
summary\$total_fragment = total_fragment 
summary\$nCells_detected = nCells_detected
summary\$median_frags = median_frags
summary\$median_TSSe = median_TSSe
summary\$mean_frags = mean_frags
summary\$mean_TSSe = mean_TSSe
f_name = paste( sample , 'summary.txt' , sep = '_')
write.table(t(summary) , file = f_name , quote = F , col.names = F)


FS_df <- plotFragmentSizes(ArchRProj = proj1,
                           returnDF = T      )


title = paste( 'Sample = ' , sample , '\n' ,
               'nFrags = ' , total_fragment , sep = '')
p_FS_df = ggplot() +
  geom_line(data = FS_df , aes(x= fragmentSize, y = fragmentPercent) , color= 'blue' ) +
  ggtitle(title) +
  theme_classic()
ggsave(  "ATAC_FragmentSize.pdf" , p_FS_df ,width = 10, height = 6 , units="in", scale=1)


TSS_df <- plotTSSEnrichment(ArchRProj = proj1,
                        returnDF = T )

p_TSS_df = ggplot() +
  geom_line(data =TSS_df , aes(x= x, y = normValue) , color= 'red' ) +
  geom_line(data =TSS_df , aes(x= x, y = smoothValue) , color= 'blue' )+
  geom_vline( xintercept = 0  , linetype = 'dashed') +
  geom_hline( yintercept = 1  , linetype = 'dashed') +
  xlab('TSS') + 
  ylab('TSS_enrichment') +
  ggtitle("${sample}") +
  theme_classic()
ggsave(  "ATAC_OverallTssEnrichment.pdf" , p_TSS_df ,width = 10, height = 6, units="in", scale=1)



title = paste( 'Sample = ' , sample ,  '\n' ,
               'nCells_detected = ' , nCells_detected  , '\n' ,
               'median_frags = ' , median_frags , '\n' ,
               'median_TSS_enrichment = ' , median_TSSe , 
               sep = '')
p_TSS_enrich_per_cell = ggplot( data = cell_df ,   aes(  y = TSSEnrichment  , x = log10( nFrags ) )    )+
  geom_pointdensity( size = 1 ) +
  scale_color_viridis() +
  geom_vline( xintercept = 2  , linetype = 'dashed') +  ### cell with fragment greater than 100
  geom_hline( yintercept = 5  , linetype = 'dashed') +  ### cell with TSS enrich greater than 5
  labs( color = "n_cells") +
  ggtitle(title) +
  theme_classic()
ggsave(  "ATAC_TSS_enrich_per_cell.pdf" , p_TSS_enrich_per_cell ,width = 10, height = 6, units="in", scale=1)


p_PromoterRatio = ggplot( data = cell_df ,   aes(  y = PromoterRatio  , x = Sample )   )+
  geom_violin( color = 'blue'  , width = 0.5)  +
  geom_boxplot( width = 0.2 ) +
  ggtitle("${sample}") +
  theme_classic()
ggsave(  "ATAC_PromoterRatio.pdf" , p_PromoterRatio  ,width = 3 , height = 5, units="in", scale=1)

p_BlacklistRatio = ggplot( data = cell_df ,   aes(  y =BlacklistRatio  , x = Sample )   )+
  geom_violin( color = 'blue'  , width = 0.5)  +
  geom_boxplot( width = 0.2 ) +
  ggtitle("${sample}") +
  theme_classic()
ggsave(  "ATAC_BlacklistRatio.pdf" , p_BlacklistRatio ,width = 3 , height = 5, units="in", scale=1)


cell_df\$barcode =  as.data.frame(str_split_fixed(row.names(cell_df)  ,'#' , 2))\$V2


cell_df = merge( cell_df , ID_file , by = 'barcode' , all.y = T )

p = ggplot(  data = cell_df   ,  aes(  y=A  , x=B , fill = nFrags )  )+
  xlim(0 , 97) +
  ylim(0 , 97) +
  scale_fill_gradientn( na.value = "#808080",  
                          colours = c("#1E90FF", "#FFFFFF","#f82239") )+
  theme_bw()+  
  geom_raster( ) +
  ggtitle("${sample}") +
  coord_fixed()
ggsave(  "ATAC_reads_per_pixel.pdf" , p ,width = 8 , height = 7.5 , units="in", scale=1)


p = ggplot(  data = cell_df   ,  aes(  y=A  , x=B , fill = log10(nFrags) )  )+
  xlim(0 , 97) +
  ylim(0 , 97) +
  scale_fill_gradientn( na.value = "#808080",  
                          colours = c("#1E90FF", "#FFFFFF","#f82239") )+
  theme_bw()+  
  geom_raster( ) +
  ggtitle("${sample}") +
  coord_fixed()
ggsave(  "ATAC_log_reads_per_pixel.pdf" , p ,width = 8 , height = 7.5 , units="in", scale=1)


p = ggplot(  data = cell_df   ,  aes(  y=A  , x=B , fill = TSSEnrichment   )  )+
  xlim(0 , 97) +
  ylim(0 , 97) +
  scale_fill_gradientn( na.value = "#808080",  
                        #  colours = c("#1E90FF", "#FFFFFF","#f82239") 
                        colours = c("#0000FF",  '#00FF00' , '#FFD700' , "#FF5800","#800000") ) +
  theme_bw()+  
  geom_raster( ) +
  ggtitle("${sample}") +
  coord_fixed()
ggsave(  "ATAC_spatial_TSSenrich.pdf" , p ,width = 8 , height = 7.5 , units="in", scale=1)


EOF

