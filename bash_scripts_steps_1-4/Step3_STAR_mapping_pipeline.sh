#Step 3 STAR mapping step pipeline

#The scripts and commands used are within each dotted line section

## N.B. replace "path" with actual file paths required to input/output directories/files
## generate "STAR_output_alignments" directory in desired location for output files.


## (A) script 1 - Genome index generation:
----------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash

#Sarah Coates May 2022
#star genome index generation script

/path/STAR/source/STAR --runThreadN 10 \
 --runMode genomeGenerate \
 --genomeDir /path/STAR_Genome_index \
 --genomeFastaFiles /path/STAR_input_genome_info/Su_DT_SLR_TGSGC_softmasked_rn.fasta \
 --sjdbGTFfile /path/STAR_input_genome_info/best.filt.annots.gtf2 \
 --sjdbOverhang 99
----------------------------------------------------------------------------------------------------------------------------------------------------


## (B) Commands to generate input files and output directories for mapping:
----------------------------------------------------------------------------------------------------------------------------------------------------
ls /path/ | grep '.*_f_.*' > forward_reads.txt
ls /path/ | grep '.*_r_.*' > reverse_reads.txt
cd /path/
ls /path/ | grep '.*_f_.*' | sed 's/_1_f_paired.fq.gz//' > STAR_prefixes.txt

while IFS=$'\t' read -r name; do mkdir "$name"_mapped ; done < /path/STAR_prefixes.txt
----------------------------------------------------------------------------------------------------------------------------------------------------


## (C) script 2 - genome mapping step
----------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#Sarah Coates May 2022
#STAR genome mapping script

#loop over lists of files and map them to reference genome using STAR
#loop over 3 files in the same loop by pasting them together.

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
----------------------------------------------------------------------------------------------------------------------------------------------------
#N.B. I had to edit one of the files to loop over with the command:
sed -i 's/ *$//' file
This removes trailing whitespace, which otherwise breaks the loop


## (D) script 3 - genome mapping for reads that did not complete due to an unknown error
----------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#Sarah Coates May 2022
#STAR genome mapping script hung reads

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
----------------------------------------------------------------------------------------------------------------------------------------------------