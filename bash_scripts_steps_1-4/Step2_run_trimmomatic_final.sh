#!/bin/bash

module load trimmomatic

ind_id=$1
fastq_1=$2
fastq_2=$3
out_dir=$4
MYADAPS='/apps/genomics/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa'

java -jar $TRIMMOMATIC PE ${fastq_1} ${fastq_2} \
    ${out_dir}${ind_id}'_f_paired.fq.gz' ${out_dir}${ind_id}'_f_unpaired.fq.gz' \
    ${out_dir}${ind_id}'_r_paired.fq.gz' ${out_dir}${ind_id}'_r_unpaired.fq.gz' \
    -threads 10 -phred33 ILLUMINACLIP:${MYADAPS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:50 AVGQUAL:20
