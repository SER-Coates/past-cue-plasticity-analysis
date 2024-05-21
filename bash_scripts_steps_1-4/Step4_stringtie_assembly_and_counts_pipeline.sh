########################
# Sarah E.R. Coates
# STEP 4  - reference-based transcript assembly and expression counts
# Analysis script and command pipeline for stringtie assembly of salt/zinc expression data
# and generation of expression counts for analysis

#replace "path" with required analysis paths to input/output files/directories

#############################################
# (A) transcriptome assembly steps
#############################################

## 1. generate input file lists, in the location containing sub-directories of the mapping files produced by STAR aligner in the previous step: 

### salt
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
ls | grep '^[BSGP]' | sed 's/_mapping//' > stringtie_salt_prefixes.txt
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

### zinc
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
ls | grep V3.* | sed 's/_mapping//' > stringtie_zinc_prefixes.txt 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

## 2. write scripts to assemble salt /zinc transcripts 

### salt transcripts
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#Script to generate transcript assembly via stringtie for salt data

while IFS='\t' read -r prefix;
do
/path/stringtie/stringtie -e -B -p 10 \
 -G /path/best.filt.annots.gff3 \
 -o /path/stringtie_salt_output_for_DE_analysis/"$prefix"/"$prefix".gtf /path/STAR_output_alignments/"$prefix"_mapping/"$prefix"_Aligned.sortedByCoord.out.bam ;
done < stringtie_salt_prefixes.txt
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


### zinc transcripts 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#!/bin/bash
#Script to generate transcript assembly via stringtie for zinc data

while IFS='\t' read -r prefix;
do
/path/stringtie/stringtie -e -B -p 10 \
 -G /path/best.filt.annots.gff3 \
 -o /path/stringtie_zinc_output_for_DE_analysis/"$prefix"/"$prefix".gtf /path/STAR_output_alignments/"$prefix"_mapping/"$prefix"_Aligned.sortedByCoord.out.bam ;
done < stringtie_zinc_prefixes.txt
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#############################################
# (B) Expression count matrix generation
#############################################

## Use python script included in the stringtie installation, prepDE.py

### 1. create input file list for the python script 
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
ls /path/stringtie_salt_output_for_DE_analysis/*/*.gtf > DEseq_prep_input_files_salt.txt 
paste stringtie_salt_prefixes.txt DEseq_prep_input_files_salt.txt > DESeq_prep_samples_list_salt.txt

ls /path/stringtie_zinc_output_for_DE_analysis/*/*.gtf > DEseq_prep_input_files_zinc.txt 
paste stringtie_zinc_prefixes.txt DEseq_prep_input_files_zinc.txt > DESeq_prep_samples_list_zinc.txt
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 2. run the python script for salt samples
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
python /path/stringtie/prepDE.py3 -i DEseq_prep_samples_list_salt.txt -

python /path/software/stringtie/prepDE.py3 -i DESeq_prep_samples_list_salt.txt 
	-g /path/stringtie_salt_output_for_DE_analysis/salt_gene_count_matrix.csv  
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 3. run the python script for zinc samples
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
python /path/stringtie/prepDE.py3 -i DEseq_prep_samples_list_zinc.txt -

python /path/software/stringtie/prepDE.py3 -i DESeq_prep_samples_list_zinc.txt 
	-g /path/stringtie_zinc_output_for_DE_analysis/zinc_gene_count_matrix.csv  
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------


