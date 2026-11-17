#!/bin/bash
#SBATCH --output=/mnt/server1/data9/zouqb/Standard_oe/merge.o
#SBATCH --error=/mnt/server1/data9/zouqb/Standard_oe/merge.e     
#SBATCH --nodelist=liulab177
#SBATCH --mem=1G

ppn=$1
mem=$2
Data=$3
Stage=$(echo "$Data" | awk -F'_' '{print $1"_"$2}')

work_1=/mnt/server2/data2/mijie/p300s/mijie/Spatial_endometrium/01_Raw_data/ATAC+RNA/${Stage}
work_2=/mnt/server2/data2/mijie/p300s/mijie/Spatial_endometrium/01_Raw_data/ATAC+RNA/
#work=/mnt/server2/data5/Spatial_endometrium

mkdir -p ${work_1}/${Data}_Add
#mkdir -p ${work}/${Data}_Add2

cd ${work_1}/${Data}_Add
#cd ${work}/${Data}_Add2
cat > ${work_1}/${Data}_Add/Merge_RawData_${Data}.sh << EOF
#!/bin/sh
#SBATCH --job-name=${Data}   
#SBATCH --cpus-per-task=$ppn     
#SBATCH --mem=${mem}G			
#SBATCH --nodelist=liulab177   
#SBATCH --output=${work_1}/${Data}_Add/${Data}_Add.o   
#SBATCH --error=${work_1}/${Data}_Add/${Data}_Add.e    
#SBATCH --time=99:10:00       
#SBATCH --partition=control,compute


##0
cd ${work_1}/${Data}_Add
#cd ${work}/${Data}_Add2
##1
mv ${work_1}/${Data}_add/${Data}_1.fq.gz ${work_1}/${Data}_add/${Data}_add_1.fq.gz
mv ${work_1}/${Data}_add/${Data}_2.fq.gz ${work_1}/${Data}_add/${Data}_add_2.fq.gz
cat ${work_1}/${Data}/${Data}_1.fq.gz ${work_1}/${Data}_add/${Data}_add_1.fq.gz > ${work_1}/${Data}_Add/${Data}_Add_1.fq.gz
cat ${work_1}/${Data}/${Data}_2.fq.gz ${work_1}/${Data}_add/${Data}_add_2.fq.gz > ${work_1}/${Data}_Add/${Data}_Add_2.fq.gz
#cat ${work}/${Data}_Add/${Data}_Add_1.fq.gz ${work}/${Data}_add2/${Data}_add2_1.fq.gz > ${work}/${Data}_Add2/${Data}_Add2_1.fq.gz
#cat ${work}/${Data}_Add/${Data}_Add_2.fq.gz ${work}/${Data}_add2/${Data}_add2_2.fq.gz > ${work}/${Data}_Add2/${Data}_Add2_2.fq.gz

##2
#rm -r ${work}/01_Raw_data/${Data}_add
#rm -r ${work}/${Data}_add2
EOF


sbatch Merge_RawData_${Data}.sh
