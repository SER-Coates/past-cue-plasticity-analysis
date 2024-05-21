#!/bin/bash

#Sarah Coates May 2022
#star genome index generation script

/path/STAR/source/STAR --runThreadN 10 \
 --runMode genomeGenerate \
 --genomeDir /path/STAR_Genome_index \
 --genomeFastaFiles /path/STAR_input_genome_info/Su_DT_SLR_TGSGC_softmasked_rn.fasta \
 --sjdbGTFfile /path/STAR_input_genome_info/best.filt.annots.gtf2 \
 --sjdbOverhang 99

