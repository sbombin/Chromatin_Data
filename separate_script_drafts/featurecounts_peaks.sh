#!/bin/bash

#SBATCH -p all
#SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 22
#SBATCH -J FeatureCounts
#SBATCH -o FeatureCounts.out
#SBATCH -e FeatureCounts.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sergei.bombin@emory.edu


module load samtools/1.14
module load miniconda3/2022.02
source activate /home/sbombin/miniconda3/envs/peaks

## bed files should be converted to SAF format to be used with featureCounts
## peaks bed files can be merged with HOMMER or bedtools merge before 

for file in $(ls *.bed | rev | cut -c 5- | rev | uniq)
do
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}'  ${file}.bed > ${file}.saf
done

for file in $(*.saf | cut -f1 -d "_" | uniq)
do
featureCounts -T 20 -p  -a ${file}*.saf -F SAF -o count_table_${file}.txt *shift.bam
done

conda deactivate
module unload samtools/1.14
module unload miniconda3/2022.02
