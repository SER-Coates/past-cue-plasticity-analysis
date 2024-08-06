# Past-cue plasticity analysis
Analysis to determine if past-cue plasticity facilitates adaptation.
Inputs are raw transcripts from an analysis of S. uniflora populations in response to zinc and salt treatments with hydroponics.
Analysis from pre-print: https://doi.org/10.1101/2024.05.06.592784 
# Steps with scripts to run
 #  Step 1: fastqc check data quality fastqc script name
  #script: Step1_fastqc_loop.sh
 #  Step 2: trimming to remove adapters with trimmomatic
#scripts: python script: Step2_run_trimomatic_master_final.py and run with bash script: Step2_run_trimmomatic_final.sh
# Step 3: mapping with STAR to s. uniflora reference - generate .BAM files
#Pipeline summary: Step_3_STAR_mapping_pipeline.sh  ##file with all steps includeing small commands
#Scripts:
         #a) Step3_STAR_index_script_final.sh
         #b) Step3_STAR_mapping_script_final.sh ##rerun
         #c) Step3_STAR_mapping_script_hung_reads.sh ##rerun certain reads if they fail in first script.
# Step 4: generate transcriptome annotation and read counts with stringtie using reference annotation
#scripts: Step4_stringtie_assembly_and_counts_pipeline.sh  ##file with all steps including small commands
# Step 5: import raw zinc/salt read counts into R and analyse data with DEseq2
#script: DEseq2_analysis2_step1_v6.R ##step1 R DEseq2 analysis script - generate dds objects, test differential expression comparisons (alpha = 0.05) and plot PCAs
#script: DEseq2_analysis2_step1_low_stringency.R ##step 1 as above but 0.1 alpha value for differential expression comparisons
# Step 6: analyse results of DEseq2 to address paper hypotheses and generate figures
#script: DEseq2_analysis2_step2_v6.R ##Using differentially expressed gene groups to test hypotheses and graph generation
#script: DEseq2_analysis2_step2_low_stringency.R ##step 2 as above but using inputdata from 0.1 alpha value script for differential expression comparisons
