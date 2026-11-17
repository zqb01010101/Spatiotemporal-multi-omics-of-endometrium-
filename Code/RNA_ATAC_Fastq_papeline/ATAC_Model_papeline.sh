#!/bin/sh
#SBATCH --output=/mnt/server1/data9/zouqb/Standard_oe/ATAC_papeline.o
#SBATCH --error=/mnt/server1/data9/zouqb/Standard_oe/ATAC_papeline.e 
#SBATCH --nodelist=liulab177
#SBATCH --mem=1G

#### Set up the environment
source /home/liulab/anaconda3/etc/profile.d/conda.sh
conda activate hjz_st_env
####


#### Set paths and core parameters for the library layout
Data_name=$1

stage=$(echo "$Data_name" | awk -F'_' '{print $1"_"$2}')
pre_name=$(echo "$Data_name" | awk -F'_' '{print $1"_"$2"_"$3}')
remaining=$(echo "$Data_name" | awk -F'_' '{for (i=5; i<=NF; i++) printf "%s%s", $i, (i<NF ? "_" : "")}')
Pure_data_name=${pre_name}'_'${remaining}    #sec_PO5_E19_240828_Add




Raw_data_dir=/mnt/server2/data2/mijie/p300s/mijie/Spatial_endometrium/01_Raw_data/ATAC+RNA/${stage}/${Data_name}/
Label=${Data_name}/
MACS2_gene_type=hs

#Whitelist=/mnt/server2/data5/Spatial_endometrium/spatial_RNA_ATAC/sec_PO7/sec_PO7_E8_ATAC_240906/barcode_whitelist.tsv
# Alternative ID file for the v6 barcode layout
# Alternative ID file for the legacy 96-channel chip in 10x loading order
ID_file=/mnt/server2/data3/zhangtzh/pnas/zhangtzh/02_Liver_space/core/combine_barcode.round2round1_index1_index2.v7_big_r1_25to97_r2_9to98.20241209.txt.IDFILE ## Barcode ID file for the v7 layout
# Alternative ID file for the legacy 96-channel chip with revised round-2 barcodes
Ref_dir=/mnt/server1/data5/zhanjh/reference/spatial/hg38/chromap_ref/




############

#### Derived settings
#Work_dir=/mnt/server1/data9/zouqb/Spatial_endometrium/${Pure_data_name}/ATAC/
Work_dir=/mnt/server3/data3/zouqb/Spatial_endometrium/${Pure_data_name}/ATAC/
mkdir -p ${Work_dir}
cd ${Work_dir}
cp ${ID_file}  ${Work_dir}/ID_file.txt
cut -f1 ${ID_file} > ${Work_dir}/barcode_whitelist.tsv
Whitelist=${Work_dir}/barcode_whitelist.tsv




Processed_data_dir=${Work_dir}/Processed_data
mkdir -p ${Processed_data_dir}
Summary_dir=${Work_dir}/Summary
mkdir -p ${Summary_dir}
summary=${Summary_dir}/${Data_name}_Summary_Report.txt

Ref_index=${Ref_dir}/index
Ref_FASTA=${Ref_dir}/hg38_chrall.fa
Ref_FASTA_index=${Ref_dir}/hg38_chrall.fa.fai

Data_name_1=${Data_name}_1
Data_name_2=${Data_name}_2
echo ${Label} > ${summary}
#########


######## Step-specific settings
#### umi_tools extract settings
umi_tools_extract_output_dir=${Processed_data_dir}/umi_tools_extract_output/
umi_tools_extract_output1=${umi_tools_extract_output_dir}/umi_tools_extract_${Data_name_1}.fq.gz
umi_tools_extract_output2=${umi_tools_extract_output_dir}/umi_tools_extract_${Data_name_2}.fq.gz
regular_expression="^(?P<umi_1>.)(?P<discard_1>ATCGGCGTACGACT){e<=1}.{8}(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=2}.{8}(?P<discard_3>CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG){e<=4}.*$"
filtered_out_output1=${umi_tools_extract_output_dir}/umi_tools_extract_filtered_out_output1_${Data_name_1}.fq.gz
filtered_out_output2=${umi_tools_extract_output_dir}/umi_tools_extract_filtered_out_output2_${Data_name_2}.fq.gz
umi_tools_extract_output2_gun=${umi_tools_extract_output_dir}/umi_tools_extract_${Data_name_2}.fq
seqkit_barcode_output=${umi_tools_extract_output_dir}/seqkit_barcode_.fq.gz
seqkit_genomic_output=${umi_tools_extract_output_dir}/asdf_.fq.gz


######

#### Trim Galore settings
trim_output_dir=${Processed_data_dir}/trim_output/
trim_basename=trim_${Data_name}
trim_output1=${trim_output_dir}/${trim_basename}_val_1.fq.gz
trim_output2=${trim_output_dir}/${trim_basename}_val_2.fq.gz
trimming_report=${trim_output_dir}/umi_tools_extract_${Data_name_2}.fq.gz_trimming_report.txt
####

#### Chromap settings
chromap_output_dir=${Processed_data_dir}/chromap_output/
chromap_output=${chromap_output_dir}/chromap_${Data_name}.sam
chromap_ref_index=${Ref_index}
chromap_ref_FASTA=${Ref_FASTA}
barcode=${seqkit_barcode_output}
chromap_summary=${chromap_output_dir}/chromap_summary.txt

sorted_chromap_output=${chromap_output_dir}/sorted_chromap_${Data_name}.bam
rmdup_chromap_output=${chromap_output_dir}/rmdup_chromap_${Data_name}.bam
rmdup_message=${chromap_output_dir}/${Data_name}_metrix.txt
####

#### MACS settings
MACS2_output_dir=${Processed_data_dir}/MACS2_output
MACS2_log=${MACS2_output_dir}/MACS2.log.txt

MACS2_chr_size=${MACS2_output_dir}/MACS2_chr_size.txt
sorted_MACS2_output=${MACS2_output_dir}/${Data_name}_treat_pileup_sorted.bdg
bTB_MACS2_output=${MACS2_output_dir}/${Data_name}_treat_pileup.sorted.bw
####

#### ArchR settings
ArchR_dir=${Processed_data_dir}/ArchR
ArchR_Rscript_dir=/mnt/server1/data9/zouqb/Spatial_ATAC_RNA_Papeline/ATAC_papeline/ATAC_ArchR_Papeline.sh

Cluster_dir=${Processed_data_dir}/Cluster/
Cluster_downstream_dir=/mnt/server1/data5/zhanjh/shell/spatial/downstream/ 
####

#### Signac settings
peak_matrix_dir=${Processed_data_dir}/Peak_matrix
sinto_fragments_output=${peak_matrix_dir}/sinto_fragments_${Data_name}.tsv.tmp
bedtools_output=${peak_matrix_dir}/sinto_fragments_${Data_name}.tsv
bgzip_output=${peak_matrix_dir}/sinto_fragments_${Data_name}.tsv.gz
####










cat > ${Work_dir}/${Data_name}_ATAC_spacial.sh <<EOF
#!/bin/sh
#SBATCH --job-name=${Data_name}  
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G 
#SBATCH --nodelist=liulab177 
#SBATCH --output=${Work_dir}/${Data_name}.o   # Standard output file
#SBATCH --error=${Work_dir}/${Data_name}.e    # Standard error file
#SBATCH --time=999:00:00      
#SBATCH --partition=control,compute



#### Set up the environment
source /home/liulab/anaconda3/etc/profile.d/conda.sh
conda activate hjz_st_env
#


## Check read quality with FastQC
#fastqc ${Raw_data_dir}/${Data_name_1}.fq.gz -q --threads 6 -o /mnt/server1/data9/zouqb/Standard_oe/
#fastqc ${Raw_data_dir}/${Data_name_2}.fq.gz -q --threads 6 -o /mnt/server1/data9/zouqb/Standard_oe/


#### Remove adapters with umi_tools
mkdir -p ${umi_tools_extract_output_dir}
umi_tools extract --extract-method=regex  --bc-pattern2="${regular_expression}" -I ${Raw_data_dir}/${Data_name_1}.fq.gz -S ${umi_tools_extract_output1} --read2-in=${Raw_data_dir}/${Data_name_2}.fq.gz --read2-out=${umi_tools_extract_output2} --filtered-out=${filtered_out_output1} --filtered-out2=${filtered_out_output2}  -L ${Work_dir}/umi_tools_extract_.log

echo -e '\n'check_barcode_pattern >> ${summary}
tail ${Work_dir}/umi_tools_extract_.log | grep Input >> ${summary}
tail ${Work_dir}/umi_tools_extract_.log | grep output >> ${summary}


# After umi_tools processing, the first 16 bases of read 2 retain the barcode.
# Use SeqKit to extract barcode positions from the umi_tools output, reusing its output directory.
gunzip  ${umi_tools_extract_output2}
seqkit --threads 6 subseq -r 1:16 ${umi_tools_extract_output2_gun} | gzip -c - > ${seqkit_barcode_output}
seqkit --threads 6 subseq -r 17:-1  ${umi_tools_extract_output2_gun} | gzip -c - >  ${seqkit_genomic_output}


mv ${seqkit_genomic_output} ${umi_tools_extract_output2}
seqkit_genomic_output=${umi_tools_extract_output2}







#### Trim adapter sequences
mkdir -p ${trim_output_dir}
trim_galore --cores  6  --paired ${umi_tools_extract_output1} ${umi_tools_extract_output2} --output_dir ${trim_output_dir} --basename ${trim_basename}
seqkit pair -1 ${trim_output1} -2 ${seqkit_barcode_output} 
mv ${trim_output_dir}/${trim_basename}_val_1.paired.fq.gz ${trim_output1}
mv ${trim_output_dir}/seqkit_barcode_.paired.fq.gz ${seqkit_barcode_output}

echo -e '\n'trim_adapter  ${trimming_report} >> ${summary}
grep 'Total reads processed:' ${trimming_report} >> ${summary}
grep 'Reads with adapters:' ${trimming_report} >> ${summary}
grep 'Reads written (passing filters):' ${trimming_report} >> ${summary}
grep 'shorter than the length cutoff' ${trimming_report} >> ${summary}






################
### Align reads with Chromap #######
################
mkdir -p ${chromap_output_dir}
chromap --num-threads 10 -x ${chromap_ref_index} -r ${chromap_ref_FASTA} -1 ${trim_output1} -2 ${trim_output2} -b ${barcode} --max-insert-size 2000 --MAPQ-threshold 30 --min-read-length 20 --barcode-whitelist ${Whitelist} --bc-error-threshold 1 --remove-pcr-duplicates-at-cell-level --SAM -o ${chromap_output} --summary ${chromap_summary}
echo -e '\n'mapping >> ${summary}
grep 'Number of reads:' ${Work_dir}/${Data_name}.e >> ${summary}
grep 'Number of output mappings (passed filters):'  ${Work_dir}/${Data_name}.e >> ${summary}
echo note:output_here_is_two_times_of_read_pairs_number >> ${summary}
# Deduplicate and index the Chromap output
## Sort
samtools sort -@ 10 ${chromap_output} > ${sorted_chromap_output}
## Remove duplicates
picard  MarkDuplicates  --REMOVE_DUPLICATES true -I ${sorted_chromap_output} -O ${rmdup_chromap_output} -M ${rmdup_message} --BARCODE_TAG CB
# Log progress
echo -e '\n'rmdup >> ${summary}
grep -E 'LIBRARY|Library' ${rmdup_message} | cut -f3 | sed ':a;N;s/\n/\t/g;ta' >> ${summary}
grep -E 'LIBRARY|Library' ${rmdup_message} | cut -f7 | sed ':a;N;s/\n/\t/g;ta' >> ${summary}
# Index the deduplicated file
samtools index ${rmdup_chromap_output}
# Summary statistics
samtools idxstats ${rmdup_chromap_output} 








# Call ATAC-seq peaks with MACS2

mkdir -p ${MACS2_output_dir}

## Adjust --nomodel, --nolambda, --keep-dup all, --SPMR, --shift -75, and -q 0.05 for the dataset. Set genome size (-g) above.
## MACS2 writes narrowPeak calls, BED summit positions, an R model script, and an XLS peak table to the output directory.
macs2 callpeak -n ${Data_name} -t ${rmdup_chromap_output} --outdir ${MACS2_output_dir} -g ${MACS2_gene_type} -B --nomodel --nolambda --keep-dup all --SPMR --shift -75 --extsize 150 -q 0.05 2>>${MACS2_log}

# Count peaks and save the count
peak_num=\$( cat ${MACS2_output_dir}/${Data_name}_peaks.narrowPeak | wc -l )
echo -ne '\n'peak_number_at_bulk_level'\t' >> ${summary}
echo \${peak_num} >> ${summary}

# Convert the ATAC-seq pileup signal to BigWig for genome browser visualization
cut -f1,2 ${Ref_FASTA_index} > ${MACS2_chr_size}
sort -k1,1 -k2,2n ${MACS2_output_dir}/${Data_name}_treat_pileup.bdg > ${sorted_MACS2_output}

bedGraphToBigWig ${sorted_MACS2_output} ${MACS2_chr_size}  ${bTB_MACS2_output}
rm ${sorted_MACS2_output} 
rm ${chromap_output}








## Switch environments
#export PATH=/mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/bin:$PATH
#source activate 
#conda deactivate
#conda activate /mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/envs/ATAC
conda deactivate
conda activate /mnt/server1/data9/zouqb/Env/anaconda/anaconda3/envs/zqb_base_env





#

# Run ArchR ATAC quality control
mkdir -p ${ArchR_dir}
bash ${ArchR_Rscript_dir} ${ArchR_dir} ${Data_name}  
#/opt/R/4.3.1/bin/Rscript /mnt/server2/data5/Spatial_endometrium/spatial_RNA_ATAC/sec_PO7/sec_PO7_E8_ATAC_240906/ArchR/sec_PO7_E8_ATAC_240906_ArchR.R   /mnt/server2/data5/Spatial_endometrium/core/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt
Rscript ${ArchR_dir}/${Data_name}_ArchR.R ${ID_file}

cat ${ArchR_dir}/${Data_name}_summary.txt >> ${summary} 
 
cp  ${Cluster_downstream_dir}/*  ${Cluster_dir}



### prepare fragment file for signac
### signac can call peak by invoke macs2 in R
conda deactivate
conda activate /mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/envs/ATAC
mkdir -p ${peak_matrix_dir}
sinto fragments -b ${rmdup_chromap_output} -t CB --shift_plus 4 --shift_minus -5  --min_mapq 30 --nproc 10 -f ${sinto_fragments_output} 

conda deactivate
conda activate /mnt/server1/data9/zouqb/Env/anaconda/anaconda3/envs/zqb_base_env
bedtools sort -i ${sinto_fragments_output} | grep -v chrM |  awk -vOFS='\t' -F '\t' '{ \$5=1 ; print \$0}' > ${bedtools_output} 

### Compress the processed fragment file
/mnt/server2/data1/Software/anaconda3/bin/bgzip ${bedtools_output} 
### Index the compressed file
/mnt/server2/data1/Software/anaconda3/bin/tabix -p bed ${bgzip_output} 

EOF


sbatch ${Work_dir}/${Data_name}_ATAC_spacial.sh
