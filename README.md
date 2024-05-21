# Past-cue plasticity analysis
Analysis to determine if past-cue plasticity facilitates adaptation.
Inputs are raw transcripts from an analysis of S. uniflora populations in response to zinc and salt treatments with hydroponics.
Analysis from pre-print: https://doi.org/10.1101/2024.05.06.592784 
# Steps with scripts to run
 #  Step 1: fastqc check data quality fastqc script name
  #script: 
 #  Step 2: trimming to remove adapters with trimmomatic
#scripts: 
# Step 3: mapping with STAR to s. uniflora reference - generate .BAM files
#scripts: 
# Step 4: generate transcriptome annotation and read counts with stringtie using reference annotation
#scripts:
# Step 5: import raw zinc/salt read counts into R and analyse data with DEseq2
#scripts:
# Step 6: analyse results of DEseq2 to address paper hypotheses and generate figures
#scripts:
