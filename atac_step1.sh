#!/bin/bash

#SBATCH -p all
#SBATCH -D .
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 20
#SBATCH -J atac_1
#SBATCH -o ATAC_step1.out
#SBATCH -e ATAC_step1.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sergei.bombin@emory.edu

### Load modules and source
# source /sw/eicc/Modules/3.2.10/init/sh  # if default sh does not work 
# source /sw/eicc/Modules/3.2.10/init/zsh #  to use zsh

module load java/1.8.0
module load samtools/1.14
module load picardtools/2.6.0
# Miniconda should be loaded as the last module
module load miniconda3/2022.02
source activate /home/sbombin/miniconda3/envs/bamqc
###

### Set up directories and variables

# Working directory
workdir = $(pwd)
# Configuration files
config_files=/home/sbombin/config_files/
# Binary
bin=/home/sbombin/bin
 
# BWA  Indexed Reference Genome
bwa_index=/home/sbombin/ref_genomes/INDX_BWA_mm39-egfp_ENSEMBL/INDX_BWA_mm39-egfp_ENSEMBL

# Reference Genome FASTA
ref_genome=/home/sbombin/ref_genomes/mm39_egfp.fa

# Reference Annotation
ref_annot==/home/sbombin/ref_genomes/mm39_egfp.gtf

## Create directories for output files
mkdir fastq_trimmed
mkdir -p fastq_trimmed/singletons
mkdir aligned
mkdir bam_filtered
mkdir bam_shifted
mkdir bam_final

for file in $(ls *.fastq | rev | cut -c 9- | rev | uniq)
do
## Trimm Illumina Adapters and low quality bases
java -Xmx20000m -jar $bin/trimmomatic-0.39.jar PE -phred33 -threads 20 ${file}_1.fastq ${file}_2.fastq fastq_trimmed/${file}_trim_1.fastq fastq_trimmed/singletons/${file}_trim_1_singletons.fastq fastq_trimmed/${file}_trim_2.fastq fastq_trimmed/singletons/${file}_trim_2_singletons.fastq ILLUMINACLIP:$config_files/full_adapters.fa:2:30:10:8:TRUE LEADING:10 TRAILING:10 MINLEN:20

## BWA MEM Alignment 
bin/bwa mem -t 20 -M $bwa_index fastq_trimmed/${file}_trim_1.fastq fastq_trimmed/${file}_trim_2.fastq > aligned_bwa/${file}.sam
cd aligned

## Convert SAM to BAM, sort, and index
samtools sort -@ 10 -o ${file}.bam -O bam ${file}.sam
samtools index -@ 10 ${file}.bam

## Mark duplicates
java -Xmx20000m -jar /sw/eicc/Pkgs/picardtools/2.6.0/picard.jar MarkDuplicates INPUT=${file}.bam OUTPUT=${file}_markdup.bam METRICS_FILE=${file}.markdup_metrics.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false CREATE_INDEX=true

## Filter BAM  files to remove low quality alignment, duplicates,and  mitochondrial reads. Sort and index after
# replace chrM with MT if  genome annotatio is in Ensembl format
samtools idxstats ${file}_markdup.bam | cut -f 1 | grep -v -P "^chrM$" | xargs samtools view -@ 20 -h -b -F 1804 -f 2 -q 20 ${file}_markdup.bam | samtools sort /dev/stdin -@ 10 -O bam -o ../bam_filtered/${file}_filtr.bam
cd ../bam_filtered/
samtools index -@ 10 ${file}_filtr.bam

## Shift alignment +4bp for positive strand and -5 for negative strand to achieve base-pair resolution of TF footprint
alignmentSieve --numberOfProcessors 20 --ATACshift --bam ${file}_filtr.bam -o ../bam_shifted/${file}_shift.bam 
cd ../bam_shifted/
# Should sort and index again
samtools sort -@ 10 -o ../bam_final/${file}_final.bam -O bam ${file}_shift.bam
cd ../bam_final/
samtools index -@ 10 ${file}_final.bam

## Convert BAM to HOMER format
makeTagDirectory ../tagdir/${file} ${file}_final.bam   
done

conda deactivate
module unload java/1.8.0
module unload samtools/1.14
module unload picardtools/2.6.0
module unload miniconda3/2022.02
