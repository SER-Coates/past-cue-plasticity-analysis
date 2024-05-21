#!/bin/bash
#script to loop over lists of files and map them to reference genome in STAR using while loop
#loop over 3 files in the same loop.
#paste them together.

#replace "path" with actual file paths required to input/output directories/files
#generate "STAR_output_alignments" directory in desired location for output files.

while IFS=$'\t' read -r f r name;
do
/path/STAR/source/STAR --genomeDir /path/STAR_Genome_index \
 --runThreadN 20 \
 --readFilesIn "$f" "$r" \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix  /path/STAR_output_alignments/"$name"_mapping/"$name"_ \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMstrandField intronMotif \
 --outSAMattributes NH HI AS nM XS \
 --quantMode GeneCounts ;
done < <(paste forward_reads.txt reverse_reads.txt STAR_prefixes.txt)
