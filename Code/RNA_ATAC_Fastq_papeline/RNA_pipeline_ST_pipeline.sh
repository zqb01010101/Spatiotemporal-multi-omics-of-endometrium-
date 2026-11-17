#!/bin/bash
if [ $# -lt 1 ]; then
    echo "error.. need args"
    exit 1
fi

ppn=$1 # 12
mem=$2 # 60 
name=$3 # E18_0910
path=$4 # /asnas/liujiang_group/yuhao/SingleCell_Space/12_E18_stain_20um_0910/scRNA
FW=$5 # /asnas/liujiang_group/yuhao/SingleCell_Space/12_E18_stain_20um_0910/scRNA/fastq/E18.5_stain_20um_RNA_2.extract.fq.gz
RV=$6 # /asnas/liujiang_group/yuhao/SingleCell_Space/12_E18_stain_20um_0910/scRNA/fastq/E18.5_stain_20um_RNA_1.extract.fq.gz
MAP=$7 # /asnas/liujiang_group/yuhao/Reference/mouse/mm10/UCSC/STAR_index
ANN=$8 # /asnas/liujiang_group/yuhao/Reference/mouse/mm10/UCSC/Mus_musculus.GRCm38.102.chr.gtf
CONT=$9 # /asnas/liujiang_group/yuhao/Reference/mouse/mm10/rDNA/STAR_index
ID_File=${10} # /asnas/liujiang_group/yuhao/SingleCell_Space/core/combine_barcode.round2round1_index1_index2.txt
gene_model=${11}
OUTPUT=st_pipeline

TMP=st_pipeline/temp

cd $path
mkdir -p ${TMP}

echo 321321
echo $ANN

st_pipeline_run.py \
  --output-folder $OUTPUT \
  --temp-folder $TMP \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --ids $ID_File \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $name \
  --htseq-no-ambiguous \
  --verbose \
  --threads $ppn \
  --log-file $OUTPUT/${name}_log.txt \
  --star-two-pass-mode \
  --no-clean-up \
  --contaminant-index $CONT \
  --disable-clipping \
  --min-length-qual-trimming 30 \
  --star-sort-mem-limit ${mem}000000000  \
  $FW $RV


echo 123123

/mnt/server2/data3/zhangtzh/pnas/Software/anaconda3/envs/my_st_env/bin/convertEnsemblToNames.py \
  --annotation $ANN \
  --output $OUTPUT/${name}_stdata.symbol.tsv \
  $OUTPUT/${name}_stdata.tsv

/mnt/server2/data1/Software_BIG/Format_trans/st_qa-new.py --input-data $OUTPUT/${name}_stdata.tsv

###Track
/mnt/server2/data1/Software/anaconda3/bin/samtools index $TMP/mapped.bam
/mnt/server2/data1/Software/anaconda3/envs/Deep_env/bin/bamCoverage -b $TMP/mapped.bam --normalizeUsing RPKM -o $TMP/mapped.fpkm.bw -p ${ppn} --binSize 10

### RSeQC
# http://rseqc.sourceforge.net/#read-distribution-py
/mnt/server2/data1/Software/anaconda3/envs/Deep_env/bin/read_distribution.py -i $TMP/mapped.bam -r $gene_model
 
# samtools
echo -e "\n## reads after remove rRNA ##"
samtools view -c $TMP/contaminated_clean.bam
echo -e "## reads after mapping ##"
samtools view -c $TMP/mapped.bam
