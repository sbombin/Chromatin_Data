#!/bin/bash

#SBATCH -p all
#SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 20
#SBATCH -J atac_2
#SBATCH -o ATAC_step2.out
#SBATCH -e ATAC_step2.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sergei.bombin@emory.edu

### Load modules and source
# source /sw/eicc/Modules/3.2.10/init/sh  # if default sh does not work 
# source /sw/eicc/Modules/3.2.10/init/zsh #  to use zsh
module load R/4.1.1
module load miniconda3/2022.02
source activate /home/sbombin/miniconda3/envs/peaks

## Path to Homer dependencies
PATH=$PATH:/home/sbombin/bin/.//bin/

# Reference Genome FASTA
ref_genome=/home/sbombin/ref_genomes/mm39_egfp.fa
# Reference Annotation
ref_annot==/home/sbombin/ref_genomes/mm39_egfp.gtf

## Create list of  control and treated samples. Needs to be changed for the  actual file names 
# Assuming metadata.csv file has samplenames in the first column 
control=$(grep "control" metadata.csv | cut -f1 -d ",")
treat=$(grep "treatment" metadata.csv | cut -f1 -d ",")

## Create directories
mkdir ../peaks_narrow

## Call Peaks with MACS2. In most cases option -c/--control should be added if there is a negative control
# Control group 
echo "${control}" | awk 'NF{print $0 "_final.bam"}' | paste -sd " " - | xargs -I {} macs2 callpeak -f BAMPE -t {}  -g mm -n Control -B --keep-dup all --extsize 200 --outdir ../peaks_narrow
# Treated group
echo "${treat}" | awk 'NF{print $0 "_final.bam"}' | paste -sd " " - | xargs -I {} macs2 callpeak -f BAMPE -t {}  -g mm -n Treat -B --keep-dup all --extsize 200 --outdir ../peaks_narrow

## Optional: filter out Encode black listed region. Not available for mm39 genome
#module load  bedtools/2.29.1
#for file in $(ls *.narrowPeak)
#do
#bedtools intersect -v -a ${file} -b /home/sbombin/references/hg38.blacklist.bed | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' |  grep -P 'chr[\dXY]+[ \t]'  | gzip -nc  > ${file}.filtr.gz 
#done
#module unload  bedtools/2.29.1

## Convert peak files into Homer format
for file in $(ls *.narrowPeak)
do
awk 'BEGIN{FS=OFS="\t"} {gsub(/\./, "+", $6)} 1' ${file} | cut -f1,2,3,4,5,6 > ${file}_hf.bed
done 

## Merge control and treat peaks
mergePeaks -d given -venn venn_merged_narrowPeak -matrix matrix_merged_narrowPeak  Control_peaks.narrowPeak_hf.bed Treat_peaks.narrowPeak_hf.bed > Merged_narrowPeak_hf.bed

cp Merged_narrowPeak_hf.bed ../tagdir
cd ../tagdir/

## Make reads count matrix with Homer
echo "${treat} ${control" | paste -sd " " - | xargs -I {} ${ref_genome} -gtf ${ref_annot} -raw -d {} > Count_merged_narrowPeak.txt

## Differential analysis with Homer and DESeq2
# Need to assign samples to Control or Tratment group after "-peaks" option (for now manually). Should be in the same order as in a count matrix
getDiffExpression.pl Count_merged_narrowPeak.txt -peaks Treat Treat Control Control -DESeq2 > DHMR_merged_narrowPeak.txt

conda deactivate
module unload miniconda3/2022.02
module unload R/4.1.1





