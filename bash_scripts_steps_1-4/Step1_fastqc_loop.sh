#!/bin/bash

#define paths for input fq/fastq file directory and output directory
path = /set/path/here/for/input/files.gz
output= /set/output/path/here/

#load required modules:
module load FastQC

#cd to correct directory for input files
cd "$path"

#Adjust the code of the for-loop to write the results to the newly created folder

#my file paths are used:

for file in *paired.fq.gz
do
fastqc -f fastq -o ${output} ${file}
done


for file in *paired.fastq.gz
do
fastqc -f fastq -o ${output} ${file}
done

 
 

