#!/bin/sh
#SBATCH --output=/mnt/server1/data9/zouqb/Standard_oe/job.o   # Standard output file
#SBATCH --error=/mnt/server1/data9/zouqb/Standard_oe/job.e    # Standard error file
#SBATCH --nodelist=liulab177
#SBATCH --mem=1G

########################
#### Set the raw data, working directory, and reference paths for the library layout.
Data_name=$1    #sec_PO5_E19_RNA_240828_Add

stage=$(echo "$Data_name" | awk -F'_' '{print $1"_"$2}')
pre_name=$(echo "$Data_name" | awk -F'_' '{print $1"_"$2"_"$3}')
remaining=$(echo "$Data_name" | awk -F'_' '{for (i=5; i<=NF; i++) printf "%s%s", $i, (i<NF ? "_" : "")}')
Pure_name=${pre_name}'_'${remaining}    #sec_PO5_E19_240828_Add



#Work_dir=/mnt/server1/data9/zouqb/Spatial_endometrium/${Pure_name}/RNA/
Work_dir=/mnt/server3/data3/zouqb/Spatial_endometrium/${Pure_name}/RNA/
#Raw_data_dir=/mnt/server2/data2/mijie/p300s/mijie/Spatial_endometrium/01_Raw_data/ATAC+RNA/${stage}
Raw_data_dir=/mnt/server2/data2/mijie/p300s/mijie/Spatial_endometrium/01_Raw_data/RNA/


Reference_dir=/mnt/server1/data9/zouqb/Reference/human/hg38/
R_path=/mnt/server2/data1/Software/anaconda3/envs/Seurat4_env2/bin/Rscript
# Alternative ID file for the v6 barcode layout
# Alternative ID file for the legacy 96-channel chip in 10x loading order
# Alternative ID file for the v7 barcode layout
# Alternative ID file for the legacy 96-channel chip with revised round-2 barcodes
ID_file=/mnt/server1/data9/zouqb/ID_file/combine_barcode.round2round1_index1_index2.v6_big_25to97_9to98.txt
######################



## The generated script activates the environment






## Derived settings
Processed_data=${Work_dir}/Processed_data/

mkdir -p ${Work_dir}
mkdir -p ${Processed_data}

cp ${ID_file} ${Work_dir}/ID_file.txt



# umi_tools extract settings
umi_tools_output_dir=${Processed_data}/umi_tools/
Raw_data1=${Raw_data_dir}/${Data_name}/${Data_name}_1.fq.gz
Raw_data2=${Raw_data_dir}/${Data_name}/${Data_name}_2.fq.gz
umi_tools_output1=${umi_tools_output_dir}/${Data_name}_umi_tools_output_1.extract.fq.gz
umi_tools_output2=${umi_tools_output_dir}/${Data_name}_umi_tools_output_2.extract.fq.gz
umi_extracted_barcode=${umi_tools_output_dir}/${Data_name}_extracted_barcode




# Barcode check settings

check_barcode_script=/mnt/server1/data9/zouqb/Spatial_ATAC_RNA_Papeline/RNA_papeline/ZQB_RNA_pipeline_check_barcode.R
check_barcode_work_dir=${Processed_data}/check_barcode



# ST Pipeline settings
st_pipeline_work_dir=${Processed_data}/ST_pipeline/
st_pipeline_temp_dir=${st_pipeline_work_dir}/temp/
exp_name=${Data_name}_ST
ref_map=${Reference_dir}/STAR_index
ref_ann=${Reference_dir}/Homo_sapiens.GRCh38.109.chr.transcript_as_exon.gtf
ref_cont=${Reference_dir}/rDNA/STAR_index
ref_gene_model=${Reference_dir}/gene_model/hg38_RefSeq.bed
ppn=8         # Thread count
mem=120        # Memory allocation

# Cluster settings
Cluster_pipeline_dir=/mnt/server1/data9/zouqb/Spatial_ATAC_RNA_Papeline/RNA_papeline/Cluster/
Plot_ncount_pipeline=${Cluster_pipeline_dir}/01_Plot_nCount.R
Cluster_pipeline=${Cluster_pipeline_dir}/02_Cluster.R
Cluster_input_file=${st_pipeline_work_dir}/st_pipeline/${Data_name}_umi_tools_output_stdata.symbol.tsv


cat > ${Work_dir}/${Data_name}_spatial_RNA.sh << EOF
#!/bin/sh
#SBATCH --job-name=${Data_name} 
#SBATCH --cpus-per-task=8
#SBATCH --mem=140G 
#SBATCH --nodelist=liulab177
#SBATCH --output=${Work_dir}/${Data_name}.o   # Standard output file
#SBATCH --error=${Work_dir}/${Data_name}.e    # Standard error file
#SBATCH --time=999:00:00      
#SBATCH --partition=control,compute

# Activate the environment
source /home/liulab/anaconda3/etc/profile.d/conda.sh
conda activate zqb_env
#conda activate zqb_base_env




# Load the environment required for the current step
#export PATH=/mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/bin:$PATH ##5G
export PATH=/mnt/server2/data1/Software/anaconda3/bin:\$PATH  ##add seq


### Extract barcodes with umi_tools
mkdir -p ${umi_tools_output_dir}

umi_tools extract --extract-method=regex --bc-pattern2="^(?P<umi_1>.)(?P<discard_1>ATCGGCGTACGACT){s<=1}.{8}(?P<discard_2>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=1}.{8}(?P<discard_3>CCCATGATCGTCCGATGCAGTCGTGCCATGAGATGTGTATAAGAGACAG){s<=2,i<=1,d<=1}.{10}(?P<discard_4>.*)$" -I ${Raw_data1} -S ${umi_tools_output1} --read2-in=${Raw_data2} --read2-out=${umi_tools_output2} -L ${umi_tools_output_dir}/umi_tools_extract.log

zcat ${umi_tools_output2} | awk 'NR%4==2{print substr(\$1,1,16)}' > ${umi_extracted_barcode} 






#### Run the Python ST Pipeline with the configured parameters 
conda deactivate
conda activate /mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/envs/my_st_env

rm -rf ${st_pipeline_work_dir}
#rm -rf ${st_pipeline_temp_dir}

mkdir -p ${st_pipeline_work_dir}
#
#cd ${Processed_data}
 
#sh /mnt/server2/data5/Spatial_endometrium/core/7_st_pipeline.sh 8 80 ${exp_name} ${Processed_data} ${umi_tools_output2} ${umi_tools_output1} ${ref_map} ${ref_ann} ${ref_cont} ${ID_file} ${ref_gene_model}
sh /mnt/server1/data9/zouqb/Spatial_ATAC_RNA_Papeline/RNA_papeline/ZQB_RNA_pipeline_ST_pipeline.sh 8 80 ${exp_name} ${st_pipeline_work_dir} ${umi_tools_output2} ${umi_tools_output1} ${ref_map} ${ref_ann} ${ref_cont} ${ID_file} ${ref_gene_model}

conda deactivate
conda activate /mnt/server1/data9/zouqb/Env/anaconda/anaconda3/envs/zqb_base_env
#### Validate barcodes in R
mkdir -p ${check_barcode_work_dir}
#${R_path} ${check_barcode_script} -b ${ID_file} -p ${check_barcode_work_dir} -n ${umi_extracted_barcode}
Rscript ${check_barcode_script} -b ${ID_file} -p ${check_barcode_work_dir} -n ${umi_extracted_barcode}





EOF

echo -e  script_write_into:'\t\t\t'${Work_dir}/${Data_name}_spatial_RNA.sh
#date

sbatch ${Work_dir}/${Data_name}_spatial_RNA.sh

echo -e  job_finished_for:'\t\t\t'${Data_name}
#date

