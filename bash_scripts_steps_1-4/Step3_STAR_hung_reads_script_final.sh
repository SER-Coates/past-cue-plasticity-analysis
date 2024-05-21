#!/bin/bash
#Sarah Coates
#redone script for STAR transcript alignment to reference genome to map  2 read pairs that failed to finish in the previous script.

#path variable is the file path to the final directories used, replace with full file path if running analysis again.

/path/STAR/source/STAR --genomeDir /path/STAR_Genome_index \
 --runThreadN 20 \
 --readFilesIn /path/PP-RNA-1_C_1_1_f_paired.fq.gz /path/PP-RNA-1_C_1_1_r_paired.fq.gz \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix  /path/STAR_output_alignments/PP-RNA-1_C_1_mapping/PP-RNA-1_C_1_ \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMstrandField intronMotif \
 --outSAMattributes NH HI AS nM XS \
 --quantMode GeneCounts

/path/STAR/source/STAR --genomeDir /path/STAR_Genome_index \
 --runThreadN 20 \
 --readFilesIn /path/V300053233_L3_HK500SILghtEAAURAAPEI-592_1_f_paired.fq.gz /path/V300053233_L3_HK500SILghtEAAURAAPEI-592_1_r_paired.fq.gz \
 --readFilesCommand gunzip -c \
 --outFileNamePrefix  /path/STAR_output_alignments/V300053233_L3_HK500SILghtEAAURAAPEI-592_mapping/V300053233_L3_HK500SILghtEAAURAAPEI-592_ \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMunmapped Within \
 --outSAMstrandField intronMotif \
 --outSAMattributes NH HI AS nM XS \
 --quantMode GeneCounts
