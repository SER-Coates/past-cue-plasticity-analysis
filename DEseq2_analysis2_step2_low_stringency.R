#########################################################################################-
#############                         SARAH COATES                            ###########-
#############           Analysis 2 of DEseq2 significantly DE loci            ###########-
#############                     Salt and zinc analysis       low stringency              ###########-
#############                                                                 ###########-
#########################################################################################-

######################################################-
#0. Set up workspace and load requirements ----
######################################################-

#clear workspace
remove(list=ls())

#set working directory if needed
#setwd()

##0.1 Install and load required packages ----

# install.packages(tidyverse)
# install.packages("data.table")
# install.packages("fplot")
# install.packages("lme4", lib = "C:/Program Files/R/R-4.3.0/library")
# installed.packages()[, c("Package", "LibPath")]
# install.packages("ggh4x")

#install colou related packages:
# install.packages(c("wesanderson", "ggsci", "viridis"))
# install.packages("rcartocolor")
# install.packages("pdftools")

#load packages - may not need all of these so can remove them, but unsure
library(DESeq2)
library(tidyverse)
library(data.table)
library(fplot)
library(lme4)
library(RColorBrewer)
library(wesanderson)
library(ggsci)
library(viridis)
library(rcartocolor)
library(pdftools)

#load ggh4x package for strip editing ability in facet wrap.
library(ggh4x)

packageVersion("DEseq2")

##0.2. Load functions and useful code ----

###a) splitting strings:

spltVector <- function(vector, pattern, position) {
  spltEach <- function(vector, pattern, position){
    unlist( strsplit(vector, pattern) )[position]
  }
  return( as.vector( mapply(spltEach, vector, pattern, position) ) )
}

#randomisation test for checking likelihood of overlapping sets of DE genes at random:
#randomisation_gene_set_test <- function(all_genes, subsample1, subsample2, overlapped_genes) { 
randomisation_gene_set_test <- function(all_genes, subsample1, subsample2, overlapped_genes) { 
  randmatch <- c() # an emtpy vector where we'll record the number of matching items during the random sampling
  for (i in 1:10000) { # 1:10000 is the number of randomisations to do, I'd probably change to 1:100000
    s1 <-sample(1:as.integer(all_genes), as.integer(subsample1), replace = F) # 1: all genes is a vector.
    s2 <-sample(1:as.integer(all_genes), as.integer(subsample2), replace = F) # as above
    randmatch <- c(randmatch, length(which(s1 %in% s2))) # append the number of matches between this round of sampling to the results vector
  }
  
  # analyse the results
  obsmatch <- as.integer(overlapped_genes) # how many matches (e.g. overexpressed genes) in your actual dataset
  
  # -- plot
  hist( randmatch, breaks = 0:(max(c(obsmatch,randmatch))+1), main=NULL, xlab="expected number" )
  abline(v = obsmatch, col="blue")
  
  #max and median of expected result
  medr <- median(randmatch)
  maxr <- max(randmatch)
  
  # -- calculate empirical p-value
  p_empir <- length(which(randmatch >= obsmatch)) / length(randmatch)
  return(list("p_empir" = p_empir, "medr" = medr, "max" = maxr))
} 

#Fisher test for two overlapping sets (of genes)
my_fisher_test = function(overlap, set1, set2, background) {
  test_m = matrix(c(overlap, set1-overlap, set2-overlap, background-(set1+set2-overlap)), nrow = 2)
  x = fisher.test(test_m, alternative = "greater")
  return(x)
}

#code on line to convert list to data frame and added some column names:
list2df_dt <- function(x, colA, colB) {
  tmp <- lapply(x, as.data.frame, stringsAsFactors = FALSE)
  tmp <- data.table::rbindlist(tmp, idcol = colA)
  colnames(tmp)[2] <-  colB
  tmp
}


#make notin operator from in operator:
`%notin%` <- Negate(`%in%`)

##0.3 set up my colour palette for later use in figures: ----

#load colours for graphs
col_coast = "#87CEFA" #lightblue for coast
col_mine = "orange"  #orange for mine
col_PPBD = "#046C9A" #dark blue for england
col_GRSA = "#F21A00" #bright_red for wales

#load colour blind friendly palette that I might use:
#https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
#original safe palette with 12 colours:
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#Maybe try this colour in the pallete instead of 117733 "#0c6735"?

#modified palette with 14 colours for bubble plots, 2 extra added from colorBlindGrey8 palette
safe_colorblind_palette1 <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#D55E00", "#6699CC", "#888888", "#661100", "#F0E442")

#6 default R colour palette colours: 
#hue_pal()(6)
# output: "#F8766D" "#B79F00" "#00BA38" "#00BFC4" "#619CFF" "#F564E3"

##0.4 load generally useful files (aka annotation file):

#read in the functional annotations table:
fun_annot <- setDT(read.csv("R_input_files/best.filt.annots.emapper_R_import.csv", sep = ",", header = T))
str(fun_annot)

##0.5 Important abbreviations

# DE = differentially expressed
# sa = Coast-W (South Aberystwyth)
# gr = Mine-W (Grogwynion mine)
# bd = Coast-E (Brean Down)
# pp = Mine-E (Priddy Pools)
# P = plasticity
# S = salt
# Z = zinc
# C = control
# C1 = control 1 (salt experiment control)
# C2 = control 2 (zinc experiment control)
# cst = coast
# min = mine
# A = Ancestral population (coast)
# D = Descendent population (mine)
# EC = evolutionary change
# sig = significant result from DE test
# no_DE = non-significant result from DE test

#############################################################-
#1.0 IMPORT RESULTS FROM DESEQ2 SCRIPT 1 with dds analysis  ----
###############################################################-

##1.1 Results from the combined data (DEalldds1) ----

##load counts file for combined experiment data (as out put from script step 1)
dds1_all_norm_counts <- setDT(read.csv("R_input_files/alldds1_norm_counts_18_05_23.csv", sep = ","))
rownames(dds1_all_norm_counts) <- dds1_all_norm_counts$X

#control1-2 comparison across experiments (as output from previous script step 1)
SA_C1_2_no_DE <- setDT(read.csv("R_input_files/LS_SA_C1_2_no_DE_17_07_24.csv"))
BD_C1_2_no_DE <- setDT(read.csv("R_input_files/LS_BD_C1_2_no_DE_17_07_24.csv"))
GR_C1_2_no_DE <- setDT(read.csv("R_input_files/LS_GR_C1_2_no_DE_17_07_24.csv"))
PP_C1_2_no_DE <- setDT(read.csv("R_input_files/LS_PP_C1_2_no_DE_17_07_24.csv"))

#SZ comparisons between both data
res_SZgrsa_no_DE <- setDT(read.csv("R_input_files/LS_res_SZ_grsa_no_DE_17_07_24.csv"))
res_SZppbd_no_DE <- setDT(read.csv("R_input_files/LS_res_SZ_ppbd_no_DE_17_07_24.csv"))
res_SZsa_no_DE <- setDT(read.csv("R_input_files/LS_res_sa_SZ_no_DE_17_07_24.csv"))
res_SZbd_no_DE <- setDT(read.csv("R_input_files/LS_res_bd_SZ_no_DE_17_07_24.csv"))
res_SZgr_no_DE <- setDT(read.csv("R_input_files/LS_res_gr_SZ_no_DE_17_07_24.csv"))
res_SZpp_no_DE <- setDT(read.csv("R_input_files/LS_res_pp_SZ_no_DE_17_07_24.csv"))
#significant DE:
res_SZgrsa_sig <- setDT(read.csv("R_input_files/LS_res_SZ_grsa_sig_17_07_24.csv"))
res_SZppbd_sig <- setDT(read.csv("R_input_files/LS_res_SZ_ppbd_sig_17_07_24.csv"))
res_SZsa_sig <- setDT(read.csv("R_input_files/LS_res_sa_SZ_sig_17_07_24.csv"))
res_SZbd_sig <- setDT(read.csv("R_input_files/LS_res_bd_SZ_sig_17_07_24.csv"))
res_SZgr_sig <- setDT(read.csv("R_input_files/LS_res_gr_SZ_sig_17_07_24.csv"))
res_SZpp_sig <- setDT(read.csv("R_input_files/LS_res_pp_SZ_sig_17_07_24.csv"))

##1.2 load salt experiment results ----

# significant results salt exp DE between population pairs within treatments:
resGR_SA_C_sig <- setDT(read.csv("R_input_files/LS_resGR_SA_C_sig_17_07_24.csv", sep = ","))
resPP_BD_C_sig <- setDT(read.csv("R_input_files/LS_resPP_BD_C_sig_17_07_24.csv", sep = ","))
resGR_SA_S_sig <- setDT(read.csv("R_input_files/LS_resGR_SA_S_sig_17_07_24.csv", sep = ","))
resPP_BD_S_sig <- setDT(read.csv("R_input_files/LS_resPP_BD_S_sig_17_07_24.csv", sep = ","))

#non-significant results salt exp between pop pairs - only import if needed

#significant results salt exp within population comparisons between treatments:
resSA_c_s_sig <- setDT(read.csv("R_input_files/LS_resSA_c_s_sig_17_07_24.csv", sep = ","))
resBD_c_s_sig <- setDT(read.csv("R_input_files/LS_resBD_c_s_sig_17_07_24.csv", sep = ","))
resGR_c_s_sig <- setDT(read.csv("R_input_files/LS_resGR_c_s_sig_17_07_24.csv", sep = ","))
resPP_c_s_sig <- setDT(read.csv("R_input_files/LS_resPP_c_s_sig_17_07_24.csv", sep = ","))

#no DE no response to salt.
resSA_c_s_noDE <- setDT(read.csv("R_input_files/LS_resSA_c_s_noDE_17_07_24.csv", sep = ","))
resBD_c_s_noDE <- setDT(read.csv("R_input_files/LS_resBD_c_s_noDE_17_07_24.csv", sep = ","))
resGR_c_s_noDE <- setDT(read.csv("R_input_files/LS_resGR_c_s_noDE_17_07_24.csv", sep = ","))
resPP_c_s_noDE <- setDT(read.csv("R_input_files/LS_resPP_c_s_noDE_17_07_24.csv", sep = ","))

##1.3 Add imports for zinc results ----

#dds1:
resZ_GR_SA_C_sig <- setDT(read.csv("R_input_files/LS_resZ_GR_SA_C_sig_17_07_24.csv", sep = ","))
resZ_PP_BD_C_sig <- setDT(read.csv("R_input_files/LS_resZ_PP_BD_C_sig_17_07_24.csv", sep = ","))
resZ_GR_SA_Z_sig <- setDT(read.csv("R_input_files/LS_resZ_GR_SA_Z_sig_17_07_24.csv", sep = ","))
resZ_PP_BD_Z_sig <- setDT(read.csv("R_input_files/LS_resZ_PP_BD_Z_sig_17_07_24.csv", sep = ","))

#non significant EC within zinc:
resZ_GR_SA_Z_noDE <- setDT(read.csv("R_input_files/LS_resZ_GR_SA_Z_noDE_17_07_24.csv", sep = ","))
resZ_PP_BD_Z_noDE <- setDT(read.csv("R_input_files/LS_resZ_PP_BD_Z_noDE_17_07_24.csv", sep = ","))

#dds2 significant results zinc exp within population comparisons between treatments:
resSA_c_z_sig <- setDT(read.csv("R_input_files/LS_resSA_c_z_sig_17_07_24.csv", sep = ","))
resBD_c_z_sig <- setDT(read.csv("R_input_files/LS_resBD_c_z_sig_17_07_24.csv", sep = ","))
resGR_c_z_sig <- setDT(read.csv("R_input_files/LS_resGR_c_z_sig_17_07_24.csv", sep = ","))
resPP_c_z_sig <- setDT(read.csv("R_input_files/LS_resPP_c_z_sig_17_07_24.csv", sep = ","))

#dds2 non-significant results zinc exp within population comparisons between treatments:
resSA_c_z_noDE <- setDT(read.csv("R_input_files/LS_resSA_c_z_noDE_17_07_24.csv", sep = ","))
resBD_c_z_noDE <- setDT(read.csv("R_input_files/LS_resBD_c_z_noDE_17_07_24.csv", sep = ","))
resGR_c_z_noDE <- setDT(read.csv("R_input_files/LS_resGR_c_z_noDE_17_07_24.csv", sep = ","))
resPP_c_z_noDE <- setDT(read.csv("R_input_files/LS_resPP_c_z_noDE_17_07_24.csv", sep = ","))

################################################################-
#2.0 FORMAT RESULTS ready for downstream analysis ----
################################################################-

##2.1 subset the datasets by columns of choice ----
colstokeep <- c("Names", "log2FoldChange", "padj")

###2.1.1 combined dataset (DEalldds1) ----

#no differential expression within control set of genes
CCsa_ns <- SA_C1_2_no_DE[, ..colstokeep]
CCbd_ns <- BD_C1_2_no_DE[, ..colstokeep]
CCgr_ns <- GR_C1_2_no_DE[, ..colstokeep]
CCpp_ns <- PP_C1_2_no_DE[, ..colstokeep]
#SZ diffs between pops and within pops (sig and ns)
SZgrsa_ns <- res_SZgrsa_no_DE[, ..colstokeep]
SZppbd_ns <- res_SZppbd_no_DE[, ..colstokeep]
SZsa_ns <- res_SZsa_no_DE[, ..colstokeep]
SZbd_ns <- res_SZbd_no_DE[, ..colstokeep]
SZgr_ns <- res_SZgr_no_DE[, ..colstokeep]
SZpp_ns <- res_SZpp_no_DE[, ..colstokeep]
#sig SZ diff
SZgrsa_sig <- res_SZgrsa_sig[, ..colstokeep]
SZppbd_sig <- res_SZppbd_sig[, ..colstokeep]
SZsa_sig <- res_SZsa_sig[, ..colstokeep]
SZbd_sig <- res_SZbd_sig[, ..colstokeep]
SZgr_sig <- res_SZgr_sig[, ..colstokeep]
SZpp_sig <- res_SZpp_sig[, ..colstokeep]

###2.1.2 salt experiment ----
#EC control 1
C1grsa_sig <- resGR_SA_C_sig[, ..colstokeep]
C1ppbd_sig <- resPP_BD_C_sig[, ..colstokeep]
#EC salt
Sgrsa_sig <- resGR_SA_S_sig[, ..colstokeep]
Sppbd_sig <- resPP_BD_S_sig[, ..colstokeep]
#control to salt DE (salt plasticity)
CSsa_sig <- resSA_c_s_sig[, ..colstokeep]
CSbd_sig <- resBD_c_s_sig[, ..colstokeep]
CSgr_sig <- resGR_c_s_sig[, ..colstokeep]
CSpp_sig <- resPP_c_s_sig[, ..colstokeep]
#no sig. salt plasticity
CSsa_ns <- resSA_c_s_noDE[, ..colstokeep]
CSbd_ns <- resBD_c_s_noDE[, ..colstokeep]
CSgr_ns <- resGR_c_s_noDE[, ..colstokeep]
CSpp_ns <- resPP_c_s_noDE[, ..colstokeep]

### 2.1.3 zinc experiment ----
#EC control 2
C2grsa_sig  <- resZ_GR_SA_C_sig[, ..colstokeep]
C2ppbd_sig <- resZ_PP_BD_C_sig[, ..colstokeep]
#EC zinc
Zgrsa_sig  <- resZ_GR_SA_Z_sig[, ..colstokeep]
Zppbd_sig  <- resZ_PP_BD_Z_sig[, ..colstokeep]
#no significant EC zinc
Zgrsa_ns <- resZ_GR_SA_Z_noDE[, ..colstokeep]
Zppbd_ns <- resZ_PP_BD_Z_noDE[, ..colstokeep]
#control zinc significant (zinc plasticity)
CZsa_sig <- resSA_c_z_sig[, ..colstokeep]
CZbd_sig <- resBD_c_z_sig[, ..colstokeep]
CZgr_sig <- resGR_c_z_sig[, ..colstokeep]
CZpp_sig <- resPP_c_z_sig[, ..colstokeep]
#no zinc plasticity
CZsa_ns <- resSA_c_z_noDE[, ..colstokeep]
CZbd_ns <- resBD_c_z_noDE[, ..colstokeep]
CZgr_ns <- resGR_c_z_noDE[, ..colstokeep]
CZpp_ns <- resPP_c_z_noDE[, ..colstokeep]

###2.1.4 check lengths: ----

####a) combined dataset

#no differential expression within control set of genes
dim(CCsa_ns)[1] 
dim(CCbd_ns)[1] 
dim(CCgr_ns)[1]
dim(CCpp_ns)[1]
#SZ diffs between pops and within pops (sig and ns)
dim(SZgrsa_ns)[1]
dim(SZppbd_ns)[1]
dim(SZsa_ns)[1]
dim(SZbd_ns)[1]
dim(SZgr_ns)[1]
dim(SZpp_ns)[1]
#sig SZ diff
dim(SZgrsa_sig)[1]
dim(SZppbd_sig)[1]
dim(SZsa_sig)[1]
dim(SZbd_sig)[1]
dim(SZgr_sig)[1]
dim(SZpp_sig)[1]

#### b) salt experiment 
#EC control 1
dim(C1grsa_sig)[1]
dim(C1ppbd_sig)[1]
#EC salt
dim(Sgrsa_sig)[1] 
dim(Sppbd_sig)[1]
#control to salt DE (salt plasticity)
dim(CSsa_sig)[1]
dim(CSbd_sig)[1]
dim(CSgr_sig)[1]
dim(CSpp_sig)[1]
#no sig. salt plasticity
dim(CSsa_ns)[1]
dim(CSbd_ns)[1]
dim(CSgr_ns)[1]
dim(CSpp_ns)[1]

#### c) zinc experiment
#EC control 2
dim(C2grsa_sig)[1]
dim(C2ppbd_sig)[1]
#EC zinc
dim(Zgrsa_sig)[1]
dim(Zppbd_sig )[1]
#no significant EC zinc
dim(Zgrsa_ns)[1]
dim(Zppbd_ns)[1]
#control zinc significant (zinc plasticity)
dim(CZsa_sig)[1]
dim(CZbd_sig)[1]
dim(CZgr_sig)[1]
dim(CZpp_sig)[1]
#no zinc plasticity
dim(CZsa_ns)[1]
dim(CZbd_ns)[1]
dim(CZgr_ns)[1]
dim(CZpp_ns)[1]


###########################################-
##2.3 rename columns in the dataset ----
###########################################-

###2.3.1 combined dataset ----

#no differential expression within control set of genes
setnames(CCsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CCsa_ns", "padj_CCsa_ns"))
setnames(CCbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CCbd_ns", "padj_CCbd_ns"))
setnames(CCgr_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CCgr_ns", "padj_CCgr_ns"))
setnames(CCpp_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CCpp_ns", "padj_CCpp_ns"))
#SZ diffs between pops and within pops (sig and ns)
setnames(SZgrsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZgrsa_ns", "padj_SZgrsa_ns"))
setnames(SZppbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZppbd_ns", "padj_SZppbd_ns"))
setnames(SZsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZsa_ns", "padj_SZsa_ns"))
setnames(SZbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZbd_ns", "padj_SZbd_ns"))
setnames(SZgr_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZgr_ns", "padj_SZgr_ns"))
setnames(SZpp_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZpp_ns", "padj_SZpp_ns"))
#sig SZ diff
setnames(SZgrsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZgrsa_sig", "padj_SZgrsa_sig"))
setnames(SZppbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZppbd_sig", "padj_SZppbd_sig"))
setnames(SZsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZsa_sig", "padj_SZsa_sig"))
setnames(SZbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZbd_sig", "padj_SZbd_sig"))
setnames(SZgr_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZgr_sig", "padj_SZgr_sig"))
setnames(SZpp_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_SZpp_sig", "padj_SZpp_sig"))

###2.3.2 salt experiment ----
#EC control 1
setnames(C1grsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_C1grsa_sig", "padj_C1grsa_sig"))
setnames(C1ppbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_C1ppbd_sig", "padj_C1ppbd_sig"))
#EC salt
setnames(Sgrsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Sgrsa_sig", "padj_Sgrsa_sig"))
setnames(Sppbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Sppbd_sig", "padj_Sppbd_sig"))
#control to salt DE (salt plasticity)
setnames(CSsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSsa_sig", "padj_CSsa_sig"))
setnames(CSbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSbd_sig", "padj_CSbd_sig"))
setnames(CSgr_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSgr_sig", "padj_CSgr_sig"))
setnames(CSpp_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSpp_sig", "padj_CSpp_sig"))
#no sig. salt plasticity
setnames(CSsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSsa_ns", "padj_CSsa_ns"))
setnames(CSbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSbd_ns", "padj_CSbd_ns"))
setnames(CSgr_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSgr_ns", "padj_CSgr_ns"))
setnames(CSpp_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CSpp_ns", "padj_CSpp_ns"))

###2.3.3 zinc experiment ----
#EC control 2
setnames(C2grsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_C2grsa_sig", "padj_C2grsa_sig"))
setnames(C2ppbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_C2ppbd_sig", "padj_C2ppbd_sig"))
#EC zinc
setnames(Zgrsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Zgrsa_sig", "padj_Zgrsa_sig"))
setnames(Zppbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Zppbd_sig", "padj_Zppbd_sig"))
#no significant EC zinc
setnames(Zgrsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Zgrsa_ns", "padj_Zgrsa_ns"))
setnames(Zppbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_Zppbd_ns", "padj_Zppbd_ns"))
#control zinc significant (zinc plasticity)
setnames(CZsa_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZsa_sig", "padj_CZsa_sig"))
setnames(CZbd_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZbd_sig", "padj_CZbd_sig"))
setnames(CZgr_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZgr_sig", "padj_CZgr_sig"))
setnames(CZpp_sig, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZpp_sig", "padj_CZpp_sig"))
#no zinc plasticity
setnames(CZsa_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZsa_ns", "padj_CZsa_ns"))
setnames(CZbd_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZbd_ns", "padj_CZbd_ns"))
setnames(CZgr_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZgr_ns", "padj_CZgr_ns"))
setnames(CZpp_ns, old = c("log2FoldChange", "padj"), new = c("log2FoldChange_CZpp_ns", "padj_CZpp_ns"))

###2.2.4 check lengths: ----

####combined dataset

#no differential expression within control set of genes
dim(CCsa_ns)[1] 
dim(CCbd_ns)[1] 
dim(CCgr_ns)[1]
dim(CCpp_ns)[1]
#SZ diffs between pops and within pops (sig and ns)
dim(SZgrsa_ns)[1]
dim(SZppbd_ns)[1]
dim(SZsa_ns)[1]
dim(SZbd_ns)[1]
dim(SZgr_ns)[1]
dim(SZpp_ns)[1]
#sig SZ diff
dim(SZgrsa_sig)[1]
dim(SZppbd_sig)[1]
dim(SZsa_sig)[1]
dim(SZbd_sig)[1]
dim(SZgr_sig)[1]
dim(SZpp_sig)[1]

####salt experiment (experiment 1) 
#EC in control 1 
dim(C1grsa_sig)[1]
dim(C1ppbd_sig)[1]
#EC salt
dim(Sgrsa_sig)[1] 
dim(Sppbd_sig)[1]
#control to salt DE (salt plasticity)
dim(CSsa_sig)[1]
dim(CSbd_sig)[1]
dim(CSgr_sig)[1]
dim(CSpp_sig)[1]
#no sig. salt plasticity
dim(CSsa_ns)[1]
dim(CSbd_ns)[1]
dim(CSgr_ns)[1]
dim(CSpp_ns)[1]

####zinc experiment (experiment 2)
#EC control 2
dim(C2grsa_sig)[1]
dim(C2ppbd_sig)[1]
#EC zinc
dim(Zgrsa_sig)[1]
dim(Zppbd_sig )[1]
#no significant EC zinc
dim(Zgrsa_ns)[1]
dim(Zppbd_ns)[1]
#control vs zinc significant (zinc plasticity)
dim(CZsa_sig)[1]
dim(CZbd_sig)[1]
dim(CZgr_sig)[1]
dim(CZpp_sig)[1]
#no sig. zinc plasticity
dim(CZsa_ns)[1]
dim(CZbd_ns)[1]
dim(CZgr_ns)[1]
dim(CZpp_ns)[1]

###################################################################################################################-
#3.0 CC GENES C1 VS C2 CONTROLS  no sig DE between experiments 1 (salt) and 2 (zinc) for all 4 populations ----
###################################################################################################################-

#no sig. DE between C1 and C2 within both coasts and both mines.
CC_cst <- merge(CCsa_ns, CCbd_ns, by = "Names")
CC_min <- merge(CCgr_ns, CCpp_ns, by = "Names")
#all 4 population pairs have no sig. DE in control:
CC_all <- merge(CC_cst, CC_min, by = "Names") #control1-control2 set of genes
CC_allNames <- CC_all$Names
CC_all_count <- length(CC_allNames) #23903 #20781 - smaller here as more were significant.

##3.1 filter all results data by this gene set ----

### 3.1.1 combined dataset ----

#SZ diffs between pops and within pops (sig and ns)
SZgrsa_ns1 <- SZgrsa_ns[which(SZgrsa_ns$Names %in% CC_allNames),]
SZppbd_ns1 <- SZppbd_ns[which(SZppbd_ns$Names %in% CC_allNames),]
SZsa_ns1 <- SZsa_ns[which(SZsa_ns$Names %in% CC_allNames),]
SZbd_ns1 <- SZbd_ns[which(SZbd_ns$Names %in% CC_allNames),]
SZgr_ns1 <- SZgr_ns[which(SZgr_ns$Names %in% CC_allNames),]
SZpp_ns1 <- SZpp_ns[which(SZpp_ns$Names %in% CC_allNames),]
#sig SZ diff
SZgrsa_sig1 <- SZgrsa_sig[which(SZgrsa_sig$Names %in% CC_allNames),]
SZppbd_sig1 <- SZppbd_sig[which(SZppbd_sig$Names %in% CC_allNames),]
SZsa_sig1  <- SZsa_sig[which(SZsa_sig$Names %in% CC_allNames),]
SZbd_sig1  <- SZbd_sig[which(SZbd_sig$Names %in% CC_allNames),] 
SZgr_sig1 <- SZgr_sig[which(SZgr_sig$Names %in% CC_allNames),]
SZpp_sig1 <- SZpp_sig[which(SZpp_sig$Names %in% CC_allNames),]

### 3.1.2 salt experiment dataset ----
#EC control 1
C1grsa_sig1 <- C1grsa_sig[which(C1grsa_sig$Names %in% CC_allNames),]
C1ppbd_sig1 <- C1ppbd_sig[which(C1ppbd_sig$Names %in% CC_allNames),]
#EC salt
Sgrsa_sig1 <- Sgrsa_sig[which(Sgrsa_sig$Names %in% CC_allNames),]
Sppbd_sig1 <- Sppbd_sig[which(Sppbd_sig$Names %in% CC_allNames),]
#control to salt DE (salt plasticity)
CSsa_sig1 <- CSsa_sig[which(CSsa_sig$Names %in% CC_allNames),]
CSbd_sig1 <- CSbd_sig[which(CSbd_sig$Names %in% CC_allNames),]
CSgr_sig1 <- CSgr_sig[which(CSgr_sig$Names %in% CC_allNames),]
CSpp_sig1 <- CSpp_sig[which(CSpp_sig$Names %in% CC_allNames),]
#no sig. salt plasticity
CSsa_ns1 <- CSsa_ns[which(CSsa_ns$Names %in% CC_allNames),]
CSbd_ns1 <- CSbd_ns[which(CSbd_ns$Names %in% CC_allNames),]
CSgr_ns1 <- CSgr_ns[which(CSgr_ns$Names %in% CC_allNames),]
CSpp_ns1 <- CSpp_ns[which(CSpp_ns$Names %in% CC_allNames),]

### 3.2.3 zinc experiment dataset ----
#EC control 2
C2grsa_sig1 <- C2grsa_sig[which(C2grsa_sig$Names %in% CC_allNames),]
C2ppbd_sig1 <- C2ppbd_sig[which(C2ppbd_sig$Names %in% CC_allNames),]
#EC zinc
Zgrsa_sig1 <- Zgrsa_sig[which(Zgrsa_sig$Names %in% CC_allNames),]
Zppbd_sig1 <- Zppbd_sig[which(Zppbd_sig$Names %in% CC_allNames),]
#no significant EC zinc
Zgrsa_ns1 <- Zgrsa_ns[which(Zgrsa_ns$Names %in% CC_allNames),]
Zppbd_ns1 <- Zppbd_ns[which(Zppbd_ns$Names %in% CC_allNames),]
#control zinc significant (zinc plasticity)
CZsa_sig1 <- CZsa_sig[which(CZsa_sig$Names %in% CC_allNames),]
CZbd_sig1 <- CZbd_sig[which(CZbd_sig$Names %in% CC_allNames),]
CZgr_sig1 <- CZgr_sig[which(CZgr_sig$Names %in% CC_allNames),]
CZpp_sig1 <- CZpp_sig[which(CZpp_sig$Names %in% CC_allNames),]
#no zinc plasticity
CZsa_ns1 <- CZsa_ns[which(CZsa_ns$Names %in% CC_allNames),]
CZbd_ns1 <- CZbd_ns[which(CZbd_ns$Names %in% CC_allNames),]
CZgr_ns1 <- CZgr_ns[which(CZgr_ns$Names %in% CC_allNames),]
CZpp_ns1 <- CZpp_ns[which(CZpp_ns$Names %in% CC_allNames),]

##3.2 Check lengths of filtered data: ----

###3.2.1 combined dataset ----

#SZ diffs between pops and within pops (sig and ns)
dim(SZgrsa_ns1)[1]
dim(SZgrsa_ns1)[1]
dim(SZppbd_ns1)[1]
dim(SZsa_ns1)[1]
dim(SZbd_ns1)[1]
dim(SZgr_ns1)[1]
dim(SZpp_ns1)[1]
#sig SZ diff
dim(SZgrsa_sig1)[1]
dim(SZppbd_sig1)[1]
dim(SZsa_sig1)[1]
dim(SZbd_sig1)[1]
dim(SZgr_sig1)[1]
dim(SZpp_sig1)[1]

###3.2.2 salt experiment dataset ----
#EC control 1
dim(C1grsa_sig1)[1]
dim(C1ppbd_sig1)[1]
#EC salt
dim(Sgrsa_sig1)[1] 
dim(Sppbd_sig1)[1]
#control to salt DE (salt plasticity)
dim(CSsa_sig1)[1]
dim(CSbd_sig1)[1]
dim(CSgr_sig1)[1]
dim(CSpp_sig1)[1]
#no sig. salt plasticity
dim(CSsa_ns1)[1]
dim(CSbd_ns1)[1]
dim(CSgr_ns1)[1]
dim(CSpp_ns1)[1]

###3.2.3 zinc experiment dataset ----
#EC control 2
dim(C2grsa_sig1)[1]
dim(C2ppbd_sig1)[1]
#EC zinc
dim(Zgrsa_sig1)[1]
dim(Zppbd_sig1)[1]
#no significant EC zinc
dim(Zgrsa_ns1)[1]
dim(Zppbd_ns1)[1]
#control zinc significant (zinc plasticity)
dim(CZsa_sig1)[1]
dim(CZbd_sig1)[1]
dim(CZgr_sig1)[1]
dim(CZpp_sig1)[1]
#no zinc plasticity
dim(CZsa_ns1)[1]
dim(CZbd_ns1)[1]
dim(CZgr_ns1)[1]
dim(CZpp_ns1)[1]

################################################################################-
#4.0 EC GENES Parallel evolutionary change / no evolutionary change ----
################################################################################-

##4.1 ECZ Evolutionary change in zinc response (EC zinc) ----

#find genes that evolve to zinc in both locations  
ECZ <- merge(Zgrsa_sig1, Zppbd_sig1, by = "Names") 
dim(ECZ)[1] # 9464
#multiply 2 LFCs - if product is > 0 then they have same direction of change
ECZ$sign_ECZ <- ECZ$log2FoldChange_Zgrsa_sig*ECZ$log2FoldChange_Zppbd_sig 
ECZ_same <- ECZ[ECZ$sign_ECZ > 0,]                                  
ECZ_same_count <- dim(ECZ_same)[1] #9121

ECZ_same_up <- ECZ_same[ECZ_same$log2FoldChange_Zgrsa_sig > 0,]
ECZ_same_down <- ECZ_same[ECZ_same$log2FoldChange_Zgrsa_sig < 0,]
dim(ECZ_same_down)

# #not using fisher test so code below commented out:
# #fisher test for data overlap:
# my_fisher_test(overlap = ECZ_same_count, set1 = dim(Zgrsa_sig1)[1],
#               set2 = dim(Zppbd_sig1)[1], background = CC_all_count)
# #p = 2.2x10-16 #significantly more than expected by chance?

##gain those genes with no EC in zinc in either set:
NO_ECZ <- merge(Zgrsa_ns1, Zppbd_ns1, by = "Names")
dim(NO_ECZ)[1] #6885 genes

#check sets are mutually exclusive:
which(ECZ$Names %in% NO_ECZ)

##4.2 ECC Evolutionary change in controls (EC control EC1/EC2) ----

#evolution in controls in salt experiment (exp1):
EC1 <- merge(C1grsa_sig1, C1ppbd_sig1, by = "Names")
dim(EC1)
EC1$sign_EC1 <- EC1$log2FoldChange_C1grsa_sig*EC1$log2FoldChange_C1ppbd_sig
EC1_same <- EC1[EC1$sign_EC1 > 0,]
dim(EC1_same)[1]

#evolution between controls in zinc experiment (exp2):
EC2 <- merge(C2grsa_sig1, C2ppbd_sig1, by = "Names")
EC2$sign_EC2 <- EC2$log2FoldChange_C2grsa_sig*EC2$log2FoldChange_C2ppbd_sig
EC2_same <- EC2[EC2$sign_EC2 > 0,]
dim(EC2_same)[1]

#evolutionary change in controls (ECC):
ECC <- merge(EC1_same, EC2_same, by = "Names")
dim(ECC)[1]
ECC$signECC <- ECC$log2FoldChange_C1grsa_sig*ECC$log2FoldChange_C2ppbd_sig
ECC_same <- ECC[ECC$signECC > 0,]
ECC_same_count <- dim(ECC_same)[1] #now 173

ECC_same_up <- ECC_same[ECC_same$log2FoldChange_C1grsa_sig > 0,]
ECC_same_down <- ECC_same[ECC_same$log2FoldChange_C1grsa_sig < 0,]
dim(ECC_same_down)

##4.3 evolutionary change in salt response if needed ----

#not needed.

##4.4 Randomisation tests: ----

# # Section commented out to prevent unnecessary runs:

# #randomisation test to see if both geog location comparisons overlap more than by chance:
# # n.b. max values are stochastic so change between repeated runs, but medians are consistent
#
#
# randomisation_gene_set_test(overlapped_genes = ECZ_same_count, subsample1 = dim(Zgrsa_sig1)[1],
#                             subsample2 = dim(Zppbd_sig1)[1], all_genes = CC_all_count)
# #
# 
# # #EC1 and EC2:
# randomisation_gene_set_test(overlapped_genes = dim(EC1_same)[1], subsample1 = dim(C1grsa_sig1)[1],
#                             subsample2 = dim(C1ppbd_sig1)[1], all_genes = CC_all_count)
# #medr 384, p empir = 0, overlap = 764
# 
# randomisation_gene_set_test(overlapped_genes = dim(EC2_same)[1], subsample1 = dim(C2grsa_sig1)[1],
#                             subsample2 = dim(C2ppbd_sig1)[1], all_genes = CC_all_count)
# #p = 0, 141 median, overlap = 278
# #randomisation test for ECC test:
# randomisation_gene_set_test(overlapped_genes = ECC_same_count, subsample1 = dim(EC1_same)[1],
#                             subsample2 = dim(EC2_same)[1], all_genes = CC_all_count)
# #p empir = 0, med = 9, overlap = 124

###########################################################################-
#5.0 COMMON PLASTICITY shared plastic genes common across populations ---- 
##########################################################################-

##5.1 Ancestral (Coastal) salt plasticity ----

ASP <- merge(CSsa_sig1, CSbd_sig1, by = "Names")
dim(ASP)[1] #959
ASP$signASP <- ASP$log2FoldChange_CSsa_sig*ASP$log2FoldChange_CSbd_sig
ASP_same <- ASP[ASP$signASP > 0,]
ASP_same_count <- dim(ASP_same)[1] #now 1085

#ASP upregulated and downregulated shared genes:

ASP_up <- ASP[ASP$log2FoldChange_CSbd_sig > 0,]
ASP_down <- ASP[ASP$log2FoldChange_CSbd_sig < 0,] 

dim(ASP_up)[1] #259
dim(ASP_down)[1] #700

##5.2 Descendant (Mine) salt plasticity ----

DSP <- merge(CSgr_sig1, CSpp_sig1, by = "Names")
dim(DSP)[1] #156
DSP$signDSP <- DSP$log2FoldChange_CSgr_sig*DSP$log2FoldChange_CSpp_sig
DSP_same <- DSP[DSP$signDSP > 0,]
DSP_same_count <- dim(DSP_same)[1] #155 #186

DSP_same_up <- DSP_same[DSP_same$log2FoldChange_CSgr_sig > 0,]
DSP_same_down <- DSP_same[DSP_same$log2FoldChange_CSgr_sig < 0,]

###5.2.1 Descendent (mine) genes with no significant salt plasticity: ----
no_DSP <-  merge(CSgr_ns1, CSpp_ns1, by = "Names")
dim(no_DSP)[1] #15477

# #unused code commented out:
# #filter for low LFCS only!  -0.1 to +0.1?
# no_DSPfilt1 <- no_DSP[no_DSP$log2FoldChange_CSgr_ns %between% c(-0.5, 0.5) ,]
# no_DSPfilt2 <- no_DSP[no_DSP$log2FoldChange_CSpp_ns %between% c(-0.5, 0.5) ,]


##5.3 Ancestral (Coastal) zinc plasticity ----

AZP <- merge(CZsa_sig1, CZbd_sig1, by = "Names")
dim(AZP)[1] #10998
AZP$signAZP <- AZP$log2FoldChange_CZsa_sig*AZP$log2FoldChange_CZbd_sig
AZP_same <- AZP[AZP$signAZP > 0,]
AZP_same_count <- dim(AZP_same)[1] #now 10890

AZP_up <- AZP[AZP$log2FoldChange_CZbd_sig > 0,]
AZP_down <- AZP[AZP$log2FoldChange_CZbd_sig < 0,] 

dim(AZP_up)[1] #4984 #5017 
dim(AZP_down)[1] #6014 #5993

###5.3.1 extra set comparing AZP and ASP directions:

ASP_AZP <- merge(ASP_same, AZP_same, by = "Names")
dim(ASP_AZP)[1] #577
ASP_AZP$signASP_AZP <- ASP_AZP$log2FoldChange_CZsa_sig*ASP_AZP$log2FoldChange_CSbd_sig

ASP_AZP_diff <- ASP_AZP[ASP_AZP$signASP_AZP < 0,]
dim(ASP_AZP_diff)[1] #411

ASP_AZP_diff_upsalt <- ASP_AZP_diff[ASP_AZP_diff$log2FoldChange_CSsa_sig > 0,] 
ASP_AZP_diff_downsalt <- ASP_AZP_diff[ASP_AZP_diff$log2FoldChange_CSsa_sig < 0,] 

#view(ASP_AZP_diff_upsalt)

##5.4 Descendant (Mine) zinc plasticity ----

DZP <- merge(CZgr_sig1, CZpp_sig1, by = "Names")
dim(DZP)[1] #143 #now 191
DZP$signDZP <- DZP$log2FoldChange_CZgr_sig*DZP$log2FoldChange_CZpp_sig
DZP_same <- DZP[DZP$signDZP > 0,]
DZP_same_count <- dim(DZP_same)[1] #190 now

DZP_same_up <- DZP_same[DZP_same$log2FoldChange_CZgr_sig > 0,]
DZP_same_down <- DZP_same[DZP_same$log2FoldChange_CZgr_sig < 0,]

###5.4.1 Derived zinc plasticity (DZP with ECZ) ----

EDZP <- merge(DZP_same, ECZ_same, by = "Names")
EDZP_count <- dim(EDZP)[1] #now 119

# #unused, DZP same as evolution direction:
#EDZP$signEDZP <- EDZP$log2FoldChange_CZgr_sig*EDZP$log2FoldChange_Zgrsa_sig
#EDZP_same <- EDZP[EDZP$signEDZP > 0,] #65 in the same direction in Zinc as evo. change.

###5.4.2 Genes not zinc plastic in mines (no DZP) ----

no_DZP <- merge(CZgr_ns1, CZpp_ns1, by = "Names")
dim(no_DZP)[1] #18394

no_DZPfilt1 <- no_DZP[no_DZP$log2FoldChange_CZgr_ns %between% c(-1, 1) ,]
no_DZPfilt2 <- no_DZP[no_DZP$log2FoldChange_CZpp_ns %between% c(-1, 1) ,]


##5.5 randomisation tests ASP, DSP, AZP, DZP ----

# # Section commented out to prevent unnecessary runs:
#
# #ASP
# randomisation_gene_set_test(overlapped_genes = ASP_same_count, subsample1 = dim(CSbd_sig1)[1],
#                             subsample2 = dim(CSsa_sig1)[1], all_genes = CC_all_count) #p_empir = 0 #median = 151
# 
# 
# #DSP
# randomisation_gene_set_test(overlapped_genes = DSP_same_count, subsample1 = dim(CSgr_sig1)[1],
#                             subsample2 = dim(CSpp_sig1)[1], all_genes = CC_all_count)
# 
# 
# #AZP
# randomisation_gene_set_test(overlapped_genes = AZP_same_count, subsample1 = dim(CZsa_sig1)[1],
#                             subsample2 = dim(CZbd_sig1)[1], all_genes = CC_all_count)
# 
# 
# #DZP
# randomisation_gene_set_test(overlapped_genes = DZP_same_count, subsample1 = dim(CZgr_sig1)[1],
#                             subsample2 = dim(CZpp_sig1)[1], all_genes = CC_all_count)
#
#
# #EDZP 
# randomisation_gene_set_test(overlapped_genes = EDZP_count, subsample1 = ECZ_same_count,
#                           subsample2 = DZP_same_count, all_genes = CC_all_count)


######################################################################################-
#6.0 Treatment-Treatment (SZ) comparisons across experiments ----
######################################################################################-

##6.1 DE between salt and zinc in Coasts ----
SZcst <- merge(SZsa_sig1, SZbd_sig1, by = "Names")
dim(SZcst)[1] #9955
SZcst$signSZcst <- SZcst$log2FoldChange_SZbd_sig*
  SZcst$log2FoldChange_SZsa_sig
SZcst_same <- SZcst[SZcst$signSZcst > 0] #9

##6.2 DE between salt and zinc in Mines ----
SZmin <- merge(SZgr_sig1, SZpp_sig1, by = "Names")
dim(SZmin)[1] #300
SZmin$signSZmin <- SZmin$log2FoldChange_SZgr_sig*
  SZmin$log2FoldChange_SZpp_sig
SZmin_same <- SZmin[SZmin$signSZmin > 0] #300


#####################################################################################################-
#7.0 Genes with Pre-adaptive plasticity ----
#####################################################################################################-

ASP_DZP <- merge(ASP_same, DZP_same, by = "Names")
ASP_DZP$signASPDZP <- ASP_DZP$log2FoldChange_CSbd_sig*ASP_DZP$log2FoldChange_CZgr_sig
ASP_DZP_same <- ASP_DZP[ASP_DZP$signASPDZP > 0,] #35

ASP_DZP_AZP <- merge(ASP_DZP_same, AZP_same, by = "Names") #20

#subset for AZP same direction as ASP and DZP:
ASP_DZP_AZP_same <- ASP_DZP_AZP[ASP_DZP_AZP$log2FoldChange_CZsa_sig*ASP_DZP_AZP$log2FoldChange_CSsa_sig > 0]

#subset those with same direction of DSP (strict pre-adaptive plasticity)
ASP_DZP_AZP_DSP <- merge(ASP_DZP_AZP_same, DSP_same, by = "Names")

dim(ASP_DZP_AZP_DSP)[1] #1 gene now in the group. 

#0 genes with pre-adaptive plasticity in strict sense.

##7.1 some additional control sets (not used in this version) ----

#SZ comparison in mines filtered for low LFC:

SZ_min_noDE <- merge(SZgr_ns1, SZpp_ns1, by = "Names")

dim(SZ_min_noDE)[1]

#view(SZ_min_noDE)

SZ_min_noDE_LFCfilt1 <- SZ_min_noDE[SZ_min_noDE$log2FoldChange_SZpp_ns %between% c(-0.1, 0.1),]
SZ_min_noDE_LFCfilt2 <- SZ_min_noDE_LFCfilt1[SZ_min_noDE_LFCfilt1$log2FoldChange_SZgr_ns %between% c(-0.1, 0.1),]
#920

#checking overlap between ancestral zinc and descendent zinc repsponse to compare numbers:

AZP_DZP <- merge(AZP_same, DZP_same, by = "Names") #90 genes
dim(AZP_DZP)[1]
AZP_DZP_same <- AZP_DZP[AZP_DZP$log2FoldChange_CZsa_sig*AZP_DZP$log2FoldChange_CZgr_sig > 0]
dim(AZP_DZP_same)[1] #57


####################################################-
#8.0 Genes showing Cue transfer pattern ----
####################################################-

#EC in zinc where genes are plastic - evolved descendent zinc plastic genes
ASP_EDZP <- merge(ASP_same, EDZP, by = "Names")
ASP_EDZP$signASPEDZP <- ASP_EDZP$log2FoldChange_CSbd_sig*ASP_EDZP$log2FoldChange_CZgr_sig
ASP_EDZP_same <- ASP_EDZP[ASP_EDZP$signASPEDZP > 0,] #32

#check the SZcst set to see which direction the change is in. It is positive, as I think would be expected.
ASP_EDZP_SZcst <- merge(ASP_EDZP_same, SZcst_same, by = "Names")
ASP_EDZP_SZcst$sign_check_szdir <- ASP_EDZP_SZcst$log2FoldChange_CSbd_sig*ASP_EDZP_SZcst$log2FoldChange_SZsa_sig
ASP_EDZP_SZcst_diff <- ASP_EDZP_SZcst[ASP_EDZP_SZcst$sign_check_szdir < 0,]
dim(ASP_EDZP_SZcst_diff)[1]

#see how many genes are opposing in the AZP gene set - filter this set to generate cue transfer final set!
ASP_EDZP_SZcst_AZP <- merge(ASP_EDZP_SZcst_diff, AZP_same) 
ASP_EDZP_SZcst_AZP$signASP_EDZP_SZcst_AZP <- ASP_EDZP_SZcst_AZP$log2FoldChange_CZsa_sig*ASP_EDZP_SZcst_AZP$log2FoldChange_CZpp_sig
ASP_EDZP_SZcst_AZP_same <- ASP_EDZP_SZcst_AZP[ASP_EDZP_SZcst_AZP$signASP_EDZP_SZcst_AZP > 0,] 
ASP_EDZP_SZcst_AZP_diff <- ASP_EDZP_SZcst_AZP[ASP_EDZP_SZcst_AZP$signASP_EDZP_SZcst_AZP < 0,]                               

#final filter on the cue_transfer_gene set
dim(ASP_EDZP_SZcst_AZP_same)[1] #0 AZP in the same direction as ASP and DZP
dim(ASP_EDZP_SZcst_AZP_diff)[1] #17 in the opposing direction

#final set of cue transfer genes:
ASP_EDZP_SZcst <- ASP_EDZP_SZcst[which(ASP_EDZP_SZcst$Names %notin% ASP_EDZP_SZcst_AZP_same$Names),]

############################################################-
##8.3 checking ancestral zinc response transfer to descendent response for comparison ----
###########################################################-

#AZP that has been modified during evolution:
AZP_EDZP <- merge(AZP_same, EDZP, by = "Names") #62 genes
dim(AZP_EDZP)[1]
AZP_EDZP_same <- AZP_EDZP[AZP_EDZP$log2FoldChange_CZsa_sig*AZP_EDZP$log2FoldChange_CZgr_sig > 0]
dim(AZP_EDZP_same)[1] #28


################################################################################################-
#9.0 Co-opted genes (parallel salt plasticity adopted into mine adaptive expression) ----
################################################################################################-

################################################################################################-
#9.0 Co-opted genes (genetically adopted) (parallel salt plasticity adopted into mine adaptive expression) ----
################################################################################################-

##9.1 Setting broader gene sets ----

#subset those with evolutionary change in controls but no zinc plasticity 
ECC_noDZP <- merge(ECC_same, no_DZP, by = "Names") #probably better to just subset - but might need LFCs downstream
dim(ECC_noDZP) #92 genes
dim(ECC)

#find those with ASP as well :)
ASP_ECC_noDZP <- merge(ASP_same, ECC_noDZP, by = "Names") #39 genes
ASP_ECC_noDZP_same <- ASP_ECC_noDZP[ASP_ECC_noDZP$log2FoldChange_C1grsa_sig*ASP_ECC_noDZP$log2FoldChange_CSsa_sig > 0]
dim(ASP_ECC_noDZP_same)[1] #38

#add SZcst comparison into set of genes:
ASP_ECC_noDZP_SZcst <- merge(ASP_ECC_noDZP_same, SZcst_same)
ASP_ECC_noDZP_SZcst$sign_check_szdir <- ASP_ECC_noDZP_SZcst$log2FoldChange_CSsa_sig*ASP_ECC_noDZP_SZcst$log2FoldChange_SZsa_sig
ASP_ECC_noDZP_SZcst_same <- ASP_ECC_noDZP_SZcst[ASP_ECC_noDZP_SZcst$sign_check_szdir > 0]
ASP_ECC_noDZP_SZcst_diff <- ASP_ECC_noDZP_SZcst[ASP_ECC_noDZP_SZcst$sign_check_szdir < 0]

dim(ASP_ECC_noDZP_SZcst_same)[1]
dim(ASP_ECC_noDZP_SZcst_diff)[1]

ASP_ECC_noDZP_SZcst_AZP <- merge(ASP_ECC_noDZP_SZcst_diff, AZP_same)
dim(ASP_ECC_noDZP_SZcst_AZP)[1]
ASP_ECC_noDZP_SZcst_AZP$sign_check_AZP <- ASP_ECC_noDZP_SZcst_AZP$log2FoldChange_CSsa_sig*ASP_ECC_noDZP_SZcst_AZP$log2FoldChange_CZbd_sig
ASP_ECC_noDZP_SZcst_AZP_same <- ASP_ECC_noDZP_SZcst_AZP[ASP_ECC_noDZP_SZcst_AZP$sign_check_AZP > 0] #1
ASP_ECC_noDZP_SZcst_AZP_diff <- ASP_ECC_noDZP_SZcst_AZP[ASP_ECC_noDZP_SZcst_AZP$sign_check_AZP < 0] #39

##9.2 Final filtered Co-option gene set: ----

#I have kept the set name of the last set used for the final filtering now to save renaming variables later.
ASP_ECC_noDZP_SZcst <- ASP_ECC_noDZP_SZcst_diff[which(ASP_ECC_noDZP_SZcst_diff$Names %notin% ASP_ECC_noDZP_SZcst_AZP_same$Names),]
upreg_coopt_no <- dim(ASP_ECC_noDZP_SZcst[ASP_ECC_noDZP_SZcst$log2FoldChange_CSsa_sig > 0])[1] #1 genes

dim(ASP_ECC_noDZP_SZcst)[1] #44

####9.3 assimilation for comparison to cooption ----

dim(ECC_noDZP)

AZP_ECC <- merge(AZP_same, ECC_same, by = "Names")
dim(AZP_ECC)
AZP_ECC_same <- AZP_ECC[AZP_ECC$log2FoldChange_CZsa_sig*AZP_ECC$log2FoldChange_C1grsa_sig > 0]

AZP_ECC_noDZP <- merge(AZP_same, ECC_noDZP, by = "Names") #47 genes - equivalent to 310 in Wood et al. I think
AZP_ECC_noDZP_same <- AZP_ECC_noDZP[AZP_ECC_noDZP$log2FoldChange_CZsa_sig*AZP_ECC_noDZP$log2FoldChange_C1grsa_sig > 0]
#14 
AZP_ECC_noDZP_ECZ <- merge(AZP_ECC_noDZP, ECZ, by = "Names")
dim(AZP_ECC_noDZP_ECZ) #43 with EC in zinc and evolved fix expression from PC

#count number of genes without ECZ change which is 8
length(which(AZP_ECC_noDZP$Names %notin% AZP_ECC_noDZP_ECZ$Names))


##################################################-
#10.0 AIM 1 results visualisation (graphs/tables) ----
##################################################-

##10.1 compare ASP vs DSP gene sets ----

#comparison of ancestral versus derived plasticity genes - LFC 
ASP_DSP_shared1 <- merge(CSsa_sig1, CSgr_sig1, by = "Names")
#241 are the same.
ASP_DSP_shared1_same <- ASP_DSP_shared1[ASP_DSP_shared1$log2FoldChange_CSsa_sig*ASP_DSP_shared1$log2FoldChange_CSsa_sig > 0,]
#all 241 in same direction.

ASP_DSP_shared2 <- merge(CSbd_sig1, CSpp_sig1, by = "Names")
#298
ASP_DSP_shared2_same <- ASP_DSP_shared2[ASP_DSP_shared2$log2FoldChange_CSbd_sig*ASP_DSP_shared2$log2FoldChange_CSpp_sig > 0,]

ASP_DSP_shared2[ASP_DSP_shared2$Names %notin% ASP_DSP_shared2_same$Names,]
#297
#1 gene DE in DSP in different direction to ancestral plasticity

ASP_DSP <- merge(ASP_same, DSP_same, by = "Names")
#paralel gene set common in ancestor and descendent = 132
ASP_DSP_same <- ASP_DSP[ASP_DSP$log2FoldChange_CSbd_sig*ASP_DSP$log2FoldChange_CSgr_sig > 0,] 
#132 out of 155 are the same plastic set of genes and change in same direction c vs s

#parallel loss of salt plasticity between both coasts and both mines
#salt plastic genes look for comparisons in those sets without plasticity in mines but plasticity in coasts
CS_diff_ASP_DSP <- setdiff(ASP$Names, DSP$Names) #827

CS_diff_cst1min1 <- setdiff(CSsa_sig1$Names, CSgr_sig1$Names) #1837
CS_diff_cst2min2  <- setdiff(CSbd_sig1$Names, CSpp_sig1$Names) #1378

length(intersect(CS_diff_cst1min1, CS_diff_ASP_DSP)) #764
length(intersect(CS_diff_cst2min2, CS_diff_ASP_DSP)) #682

#export table of gene names for this set:

##10.1 bar Graphs ----

###Data table for graphs:
# Maybe change how I am presenting this data and combine with parallel change in salt response set
pop_names <- c("Coast-W", "Coast-E", "Mine-W", "Mine-E", "Coast-B", "Mine-B")
Geography <- c("A", "B", "A")
ecotype <- c("Coast", "Coast", "Mine", "Mine", "Coast", "Mine")
DE_genes_salt <- c(dim(CSsa_sig1)[1], dim(CSbd_sig1)[1], dim(CSgr_sig1)[1], 
                   dim(CSpp_sig1)[1], dim(ASP_same)[1], dim(DSP_same)[1])
DE_genes_zinc <- c(dim(CZsa_sig1)[1], dim(CZbd_sig1)[1], dim(CZgr_sig1)[1], 
                   dim(CZpp_sig1)[1], dim(AZP_same)[1], dim(DZP_same)[1])
expr_table <- data.frame(pop_names, ecotype, DE_genes_salt, DE_genes_zinc)

# #both combined somehow:
# DE_genes_both <- c(dim(CSsa_sig1)[1], dim(CSbd_sig1)[1], dim(CSgr_sig1)[1], 
#                    dim(CSpp_sig1)[1], dim(ASP_same)[1], dim(DSP_same)[1], dim(CZsa_sig1)[1], dim(CZbd_sig1)[1], dim(CZgr_sig1)[1], 
#                    dim(CZpp_sig1)[1], dim(AZP_same)[1], dim(DZP_same)[1])
# 
# pop_names1 <- c("Coast1", "Coast2", "Mine1", "Mine2", "Coast1&2", "Mine1&2", "Coast1", "Coast2", "Mine1", "Mine2", "Coast1&2", "Mine1&2")
# Geography1 <- c("A", "B", "A", "A", "B", "A")
# ecotype1 <- c("Coast", "Coast", "Mine", "Mine", "Coast", "Mine", "Coast", "Coast", "Mine", "Mine", "Coast", "Mine")
# treatment <- c("Salt", "Salt", "Salt", "Salt", "Salt", "Salt", "Zinc", "Zinc", "Zinc", "Zinc", "Zinc", "Zinc")

#expr_table1 <- data.frame(pop_names1, ecotype1, treatment, DE_genes_both)

#both graph
# DE_genes_both_graph <- expr_table1 %>% mutate(pop_names=fct_relevel(pop_names1, c(c("Coast-W", "Coast-E", "Mine-W", "Mine-E", "Coast-W&E", "Mine-W&E")))) %>%
#   ggplot(., aes(x= pop_names1, y=DE_genes_both, fill = ecotype1))+
#   geom_text(aes(label=DE_genes_both), vjust=-1)+
#   geom_col()+
#   xlab("Population")+
#   ylab("DE gene number in salt response")+
#   scale_fill_manual(values = c(col_coast, col_mine)) +
#   facet_wrap(~ treatment, nrow= 1, scale = "fixed") +
#   #theme_classic()+
#   theme(panel.background = element_rect(colour = "grey40", fill = NA), panel.grid = element_blank(), 
#         strip.background = element_rect(colour = "grey40", fill = NA), strip.placement = "outside", strip.text = element_text(size = 18),
#         aspect.ratio = 1, axis.title = element_text(size = 15), axis.text = element_text(size = 14))
# DE_genes_both_graph


#salt graph same pattern as before with 0.05 alpha tests, just a few more genes in each group.
#set the order of factor levels to fix order first, then plot this order of info in bar chart
salt_DE_genes <- expr_table %>% mutate(pop_names=fct_relevel(pop_names, c("Coast-W", "Mine-W", "Coast-E", "Mine-E", "Coast-B", "Mine-B"))) %>%
  ggplot(., aes(x= pop_names, y=DE_genes_salt, fill = ecotype))+
  geom_col()+
  xlab("Population")+
  ylab("DE gene number in salt response")+
  scale_fill_manual(values = c(col_coast, col_mine)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
        axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 22), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.title = element_text(size = 22),
        aspect.ratio = 1)
#salt_DE_genes
# 
# 
# #zinc graph
# zinc_DE_genes <- expr_table %>% mutate(pop_names=fct_relevel(pop_names, c("Coast-W", "Mine-W", "Coast-E", "Mine-E", "Coast-B", "Mine-B"))) %>%
#   ggplot(., aes(x= pop_names, y=DE_genes_zinc, fill = ecotype))+
#   geom_col()+
#   xlab("Population")+
#   ylab("DE gene number in zinc response")+
#   scale_fill_manual(values = c(col_coast, col_mine)) +
#   theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
#         axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 15), axis.text = element_text(size = 13),
#         legend.text = element_text(size = 15), legend.title = element_text(size = 14),
#         aspect.ratio = 1.2)
# zinc_DE_genes


# #export to files


setFplot_page(page = "a4", margins = "normal", units = "tw", pt = 20, reset = FALSE, w2h = 1)
#
# pdf_fit(file = "salt_DE_genes_all_27_02_24_v10.pdf", pt =20, width = 1.15, w2h = 1)
# salt_DE_genes
# fit.off()

# #try just exporting with pdf function:
# pdf(file = "salt_DE_genes_all_27_02_24_v3.pdf")
# salt_DE_genes
# dev.off()

# 
# pdf(file = "zinc_DE_genes_all_17_11_23.pdf")
# zinc_DE_genes
# dev.off()

#original settings used:
# pdf(width = 14, file = "both_DE_genes_plot_15_11_23.pdf")
# DE_genes_both_graph
# dev.off()

#try a more precise text size export:


##10.2 results tables to use instead of bar charts ----

#make results tables for salt and zinc plasticity changes:
geographical_location <- c("A", "B", "shared")
Coastal_salt_plastic_genes<- c(dim(CSsa_sig1)[1], dim(CSbd_sig1)[1], dim(ASP_same)[1]) 
Mine_salt_plastic_genes <- c(dim(CSgr_sig1)[1], dim(CSpp_sig1)[1],  dim(DSP_same)[1])
Coastal_zinc_plastic_genes<- c(dim(CZsa_sig1)[1], dim(CZbd_sig1)[1], dim(AZP_same)[1]) 
Mine_zinc_plastic_genes <- c(dim(CZgr_sig1)[1], dim(CZpp_sig1)[1],  dim(DZP_same)[1])

#create salt table to see how salt plasticity changes during adaptation:
salt_table1 <- data.frame(geographical_location, Coastal_salt_plastic_genes, Mine_salt_plastic_genes) 

#create zinc table for comparison that can go in the supplementary info.
zinc_table1 <- data.frame(geographical_location, Coastal_zinc_plastic_genes, Mine_zinc_plastic_genes)

#export results tables:
#write.csv(salt_table1, "salt_plasticity change_22_05_23.csv")
#write.csv(zinc_table1, "zinc_plasticity change_22_05_23.csv")


################################-
#10.3 LFC plots ----
################################-

# #LFCS:
# #histograms/boxes summarising LFC distributions:
# boxplot(CSsa_sig1$log2FoldChange_CSsa_sig, ylim = c(-10,10), ylab = "log2 fold change s/c", xlab = "sa")
# boxplot(CSbd_sig1$log2FoldChange_CSbd_sig, ylim = c(-10,10), xlab = "bd")
# boxplot(CSgr_sig1$log2FoldChange_CSgr_sig, ylim = c(-10,10), xlab = "gr")
# boxplot(CSpp_sig1$log2FoldChange_CSpp_sig, ylim = c(-10,10), xlab = "pp")
# 
# boxplot(ASP$log2FoldChange_CSsa_sig, ylim = c(-10,10))
# boxplot(ASP$log2FoldChange_CSbd_sig, ylim = c(-10,10))
# boxplot(DSP$log2FoldChange_CSgr_sig, ylim = c(-10,10))
# boxplot(DSP$log2FoldChange_CSpp_sig, ylim = c(-10,10))
# 
# boxplot(ASP_DSP_same$log2FoldChange_CSsa_sig, ylim = c(-4,4))
# boxplot(ASP_DSP_same$log2FoldChange_CSgr_sig, ylim = c(-4,4))
# boxplot(ASP_DSP_same$log2FoldChange_CSbd_sig, ylim = c(-4,4))
# boxplot(ASP_DSP_same$log2FoldChange_CSpp_sig, ylim = c(-4,4))
# 
# #par(mfrow=c(1,4))
# par(mfrow=c(1,1))


#########################################################-
#11.0 MAIN AIM Visualising results ----
#########################################################-

##11.1  Make overall table of results for visuals using normalised counts ----

CC_norm_counts2 <- as.data.frame(dds1_all_norm_counts[which(rownames(dds1_all_norm_counts) %in% CC_allNames),])
colnames(CC_norm_counts2)
str(CC_norm_counts2)
dim(CC_norm_counts2)[1]

head(CC_norm_counts2)

#rename some columns to make them easier to transform.  
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.6_C_5"] <- "GRRNA6_C_5"
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.6_S_6"] <- "GRRNA6_S_6"
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.10_C_13"] <- "GRRNA10_C_13"
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.10_S_14"] <- "GRRNA10_S_14"
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.12_C_23"] <- "GRRNA12_C_23"
names(CC_norm_counts2)[names(CC_norm_counts2) == "GR.RNA.12_S_24"] <- "GRRNA12_S_24"
names(CC_norm_counts2)[names(CC_norm_counts2) == "PP.RNA.1_C_1"] <- "PPRNA1_C_1"
names(CC_norm_counts2)[names(CC_norm_counts2) == "PP.RNA.1_S_2"] <- "PPRNA1_S_2"

#write out for next step in new analysis script for Pst
#write.csv(CC_norm_counts2, "CC_norm_counts2_20_06_24.csv")

CC_norm_counts3 <- pivot_longer(CC_norm_counts2, 
                                cols = c(SA02_C_7,  SA04_C_15, SA07_C_11, SA6C, SA7C, SA8C,
                                         BD03_C_17,  BD05_C_9, BD07_C_19, BD1C, BD11C, BD12C,
                                         SA02_S_8, SA04_S_16, SA07_S_12, SA6Z, SA7Z, SA8Z,
                                         BD03_S_18, BD05_S_10, BD07_S_20, BD1Z, BD11Z, BD12Z,
                                         GRRNA6_C_5, GRRNA10_C_13, GRRNA12_C_23, GR10C, GR2C, GR8C, 
                                         PPRNA1_C_1, PP1_C_3, PP12_C_21, PP1C, PP4C, PP8C, 
                                         GRRNA6_S_6, GRRNA10_S_14, GRRNA12_S_24, GR10Z, GR2Z, GR8Z, 
                                         PPRNA1_S_2, PP1_S_4, PP12_S_22, PP1Z, PP4Z, PP8Z), 
                                names_to = c("individual", "treatment"), values_to = "norm_count",
                                names_pattern = "(^[A-Z]{2,5}[0-9]{0,2}.)([A-Z])")

#glitchy way of renaming the columns - kept breaking for unknown reasons
#CC_norm_counts3 <- rename(CC_norm_counts3, "gene_name" = "X")

#colnames(CC_norm_counts3)

colnames(CC_norm_counts3)[1] <- c("gene_name")

#colnames(CC_norm_counts3)

dim(CC_norm_counts3)[1]

head(CC_norm_counts3)

#write out for next step in new analysis script for Pst
#write.csv(CC_norm_counts3, "CC_norm_counts3_20_06_24.csv")


# # checking norm vs log norm count distribution
# hist(CC_norm_counts3$norm_count)
# hist(log(CC_norm_counts3$norm_count))
# # tried out mean with logs but ends up giving strange means for those with 0 values in the logged data.
# e.g. this gene goes weird and then I can't plot it. "evm.TU.s_12572.2883"
# # I am sticking to the original plan of arithmetic mean of the raw counts within treatments between ind/
# CC_norm_counts3$log_norm_count <- log(CC_norm_counts3$norm_count)
# 
# hist(CC_norm_counts3$log_norm_count)

#ecotype factor much quicker method than the loop:
ecotype <- rep(c("cst", "cst", "cst", "cst", "cst", "cst", 
  "cst", "cst", "cst", "cst", "cst", "cst",
  "cst", "cst", "cst", "cst", "cst", "cst",
  "cst", "cst", "cst", "cst", "cst", "cst",
  "min", "min", "min", "min", "min", "min",
  "min", "min", "min", "min", "min", "min",
  "min", "min", "min", "min", "min", "min",
  "min", "min", "min", "min", "min", "min"), 20781)
length(ecotype)

#geog factor made much quicker method than the loop:
geography <- rep(c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
                  "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
                  "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
                  "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"), 20781)

length(geography)


#add factor for geography, ecotype, combined eco-treat and ID
CC_norm_counts3$ecotype <- ecotype
CC_norm_counts3$ecotype_treatment <- paste(CC_norm_counts3$ecotype, 
                                           CC_norm_counts3$treatment, sep = "")
CC_norm_counts3$gene_ecotype_treatment <- paste(CC_norm_counts3$gene_name, 
                                                CC_norm_counts3$ecotype_treatment, sep = ",") #important!!
CC_norm_counts3$geog_pair <- geography

#most important factor column for next steps:
CC_norm_counts3$gene_eco_treat_geog <- paste(CC_norm_counts3$gene_ecotype_treatment, 
                                             CC_norm_counts3$geog_pair, sep = ",")

CC_norm_counts3$gene_name_geog_pair <- paste(CC_norm_counts3$gene_name, CC_norm_counts3$geog_pair)
CC_norm_counts3$eco_treat <- spltVector(CC_norm_counts3$gene_eco_treat_geog, ",", 2)


################################################################################-
##11.2 plot gene counts per gene per location (for cue transfer/cooption) ----
################################################################################-

#These plots are for supplementary info:

### create tables as subsets of the greater table for the plots:
#Cue transfer dataset:
Cue_transfer_norm_counts <- CC_norm_counts3[which(CC_norm_counts3$gene_name 
                                                 %in% ASP_EDZP_SZcst$Names),]

#cooption dataset:
Cooption_norm_counts <- CC_norm_counts3[which(CC_norm_counts3$gene_name 
                                              %in% ASP_ECC_noDZP_SZcst$Names),]

#view(Cooption_norm_counts)

### plot the data 
# #Cue transfer plots with all stuff needed
Cue_transfer_lineplots1 <- ggplot(Cue_transfer_norm_counts, aes(x = eco_treat, y = log(norm_count), color = geog_pair)) +
  geom_jitter(aes(color = geog_pair), position = position_dodge(0.2), cex = 1)+
  xlab("ecotype and treatment")+
  ylab("normalised gene counts")+
  scale_color_manual(values = c(col_GRSA, col_PPBD))+
  theme(panel.background = element_rect(colour = "grey40", fill = NA, linewidth = 1), panel.grid = element_blank(),
        strip.text = element_text(size = 10), strip.background = element_rect(colour = "grey40", fill = NA, linewidth = 1),
        aspect.ratio = 0.7, axis.line = element_blank(),
        axis.title = element_text(size = 10), axis.text = element_text(size = 8), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))+
  facet_wrap(~ gene_name, nrow= 8, scale = "free")
Cue_transfer_lineplots1

# pdf(file = "Cue_transfer_lineplots1_26_02_24.pdf", width = 14, height = 16)
# Cue_transfer_lineplots1
# dev.off()

# pdf_fit(file = "Cue_transfer_lineplots1_26_02_24_v2.pdf", w2h = 1, pt = 4)
# Cue_transfer_lineplots1
# fit.off()

# #ASP EDZP plots with all stuff needed
Cooption_lineplots1 <- ggplot(Cooption_norm_counts, aes(x = eco_treat, y = log(norm_count), color = geog_pair)) +
  geom_jitter(aes(color = geog_pair), position = position_dodge(0.2), cex = 1)+
  xlab("ecotype and treatment")+
  ylab("normalised gene counts")+
  scale_color_manual(values = c(col_GRSA, col_PPBD))+
  theme(panel.background = element_rect(colour = "grey40", fill = NA, linewidth = 1), panel.grid = element_blank(),
        strip.text = element_text(size = 10), strip.background = element_rect(colour = "grey40", fill = NA, linewidth = 1),
        aspect.ratio = 0.7, axis.line = element_blank(),
        axis.title = element_text(size = 10), axis.text = element_text(size = 8), 
        legend.background = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))+
  facet_wrap(~ gene_name, ncol = 6, scale = "free")
Cooption_lineplots1

# pdf(file = "Cooption_lineplots1_26_02_24.pdf", width = 14, height = 16)
# Cooption_lineplots1
# dev.off()


##############################################################################################-
##11.3 lineboxplots with means per ecotype and all genes in one graph ----
##############################################################################################-

### set up overall ecotype means table:

#edit this to be mean of the logged counts.

#version with non logged with arithmetic mean of raw data:
CC_norm_counts3_mean_ecotreat <- CC_norm_counts3 %>%
  group_by(gene_ecotype_treatment) %>%
  summarise_at(vars(norm_count), list(mean_norm_count = mean, SD_norm_count = sd, median_norm_count = median, max_norm_count = max, min_norm_count = min))

# CC_norm_counts3_mean_ecotreat <- CC_norm_counts3 %>%
#   group_by(gene_ecotype_treatment) %>%
#   summarise_at(vars(log_norm_count, norm_count), list(mean = mean, SD = sd, median = median, max = max, min = min))

#generate new columns from splitting the first column
CC_norm_counts3_mean_ecotreat$gene_name <- spltVector(CC_norm_counts3_mean_ecotreat$gene_ecotype_treatment,
                                                      ",", 1)
CC_norm_counts3_mean_ecotreat$eco_treat <- spltVector(CC_norm_counts3_mean_ecotreat$gene_ecotype_treatment, ",", 2)

CC_norm_counts3_mean_ecotreat$ecotype <- spltVector(CC_norm_counts3_mean_ecotreat$eco_treat, "[A-Z]", 1)

head(CC_norm_counts3_mean_ecotreat)

###11.3.1 tables to test Aim 1 ----

ASP_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_same$Names),]
remove_zinc <- c("cstZ", "minZ")
ASP_table <- ASP_table[which(ASP_table$eco_treat %notin% remove_zinc),] 
ASP_table$pattern <- if_else(ASP_table$gene_name %in% DSP_same$Names, "DSP_same", 
                                      if_else(ASP_table$gene_name %in%  no_DSP$Names, "noDSP", "other"))

ASP_DSP_table <- ASP_table[which(ASP_table$pattern == "DSP_same"),]
ASP_noDSP_table <- ASP_table[which(ASP_table$pattern == "noDSP"),]
ASP_other_patterns_table <- ASP_table[which(ASP_table$pattern == "other"),]

#comparing salt plasticity for mine genes and coast ones out of shared subset:
DSP_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% DSP_same$Names),]
DSP_table$pattern <- if_else(DSP_table$gene_name %in% ASP_DSP_same$Names, "ASP_DSP_same", "no_ASP_notsame")
DSP_table <- DSP_table[which(DSP_table$eco_treat %notin% remove_zinc),] 
DSP_noASP_table <- DSP_table[which(DSP_table$pattern == "no_ASP_notsame"),]


###11.3.2 tables for Cue transfer ----

#subset the broader set to generate the sets for the specific hypotheses of interest:
Cue_transfer_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_EDZP_SZcst$Names),]

#Cue_transfer_table$pattern <- if_else(Cue_transfer_table$gene_name %in% ASP_EDZP_SZcst_DSP$Names, "Harmonisation", 
#        if_else(Cue_transfer_table$gene_name %in% ASP_EDZP_SZcst_no_DSP$Names, "Cue-switching", "Other"))

Cue_transfer_table1 <- Cue_transfer_table
#Cue_transfer_table1$pattern <- factor(Cue_transfer_table1$pattern, levels = c("Harmonisation", "Cue-switching", "Other"))
Cue_transfer_table1$eco_treat <- factor(Cue_transfer_table1$eco_treat, levels = c("cstS", "cstC", "cstZ", "minS", "minC", "minZ"))

Cue_transfer_table2 <- Cue_transfer_table1
Cue_transfer_table2$eco_treat <- factor(Cue_transfer_table2$eco_treat, levels = c("cstC", "cstS", "cstZ", "minC", "minS", "minZ"))

# #smaller individual sets for individual graphs:
# ASP_EDZP_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_EDZP_same$Names),]
# groups_to_remove <- c("cstZ", "minS")
# ASP_EDZP_only_mean_ecotreat <- ASP_EDZP_mean_ecotreat[which(ASP_EDZP_mean_ecotreat$eco_treat %notin% groups_to_remove),] 
# #smaller gene sets:
# harmonisation_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_EDZP_SZcst_DSP$Names),]
# plast_switch_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_EDZP_noDSP$Names),]
# plast_switch2_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_EDZP_SZcst_no_DSP$Names),]

###11.3.3 tables for Cooption ----

cooption_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_ECC_noDZP_SZcst$Names),]
#cooption_table$pattern <- if_else(cooption_table$gene_name %in% ASP_ECC_noDZP_SZcst_noDSP$Names, "full", 
#                                      if_else(cooption_table$gene_name %in% ASP_ECC_noDZP_SZcst_DSP$Names, "partial", "other"))
#cooption_table$pattern <- factor(cooption_table$pattern, levels = c("full", "partial", "other"))

#make copy where factors are changed:
cooption_table1 <- cooption_table
#cooption_table1$pattern <- factor(cooption_table1$pattern, levels = c("full", "partial", "other"))
cooption_table1$eco_treat <- factor(cooption_table1$eco_treat, levels = c("cstS", "cstC", "cstZ", "minS", "minC", "minZ"))

cooption_table2 <- cooption_table1
cooption_table2$eco_treat <- factor(cooption_table2$eco_treat, levels = c("cstC", "cstS", "cstZ", "minC", "minS", "minZ"))


# #subpatterns:
# partial_cooption_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_ECC_noDZP_SZcst_DSP$Names),]
# partial_cooption_mean_ecotreat_up <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_ECC_noDZP_SZcst_DSP_up$Names),]
# partial_cooption_mean_ecotreat_down <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_ECC_noDZP_SZcst_DSP_down$Names),]
# #fully coopted/assimilated genes:
# full_cooption_mean_ecotreat <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_ECC_noDZP_SZcst_noDSP$Names),]

### 11.3.4 combined cue transfer and cooption table: ----

#The combined graph was ugly so code not used.
# colnames(Cue_transfer_table1)
# colnames(cooption_table1)
# 
# #add new column to become the factor to grid based on:
# setname_cue_transfer <- rep(c("Cue_transfer"), dim(Cue_transfer_table1)[1])
# Cue_transfer_table1$setname <- setname_cue_transfer
# 
# setname_cooption <- rep(c("Cooption"), dim(cooption_table1)[1])
# cooption_table1$setname <- setname_cooption
# 
# plast_facil_table <- rbind(Cue_transfer_table1, cooption_table1)
# 
# plast_facil_table
# 
# plast_facil_table$setname <- factor(plast_facil_table$setname, levels = (c("Cue_transfer", "Cooption")))
# 
# str(plast_facil_table)

###############################################-
###11.3.5 create table for the single strict pre-adaptive plastic gene
###############################################-

preadapt_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_DZP_AZP_DSP$Names),]
preadapt_table2 <- preadapt_table
preadapt_table2$eco_treat <- factor(preadapt_table2$eco_treat, levels = c("cstC", "cstS", "cstZ", "minC", "minS", "minZ"))

#looser preadaptive set:
Lpreadapt_table <- CC_norm_counts3_mean_ecotreat[which(CC_norm_counts3_mean_ecotreat$gene_name %in% ASP_DZP_AZP_same$Names),]
Lpreadapt_table2 <- Lpreadapt_table
Lpreadapt_table2$eco_treat <- factor(Lpreadapt_table2$eco_treat, levels = c("cstC", "cstS", "cstZ", "minC", "minS", "minZ"))

###############################################-
##11.4 Plot the line/box plots  ----
###############################################-

#prep:

# Only colour strips in x-direction
strip <- strip_themed(background_x = elem_list_rect(fill = c(col_coast, col_mine), color = c(F,F)))

#make labeller for the plots for each panel on facet wrap: 
#(name it as the variable name used for wrapping .labs and use the names of the key to be the current labels)
ecotype.labs <- c("Coast Ecotype", "Mine Ecotype")
names(ecotype.labs) <- c("cst", "min")


###11.4.1 plot Cue transfer ----


##alternate order of factors with control first:
Cue_transfer_simple_lineboxplot2 <- ggplot(Cue_transfer_table2, aes(x = eco_treat, y = log(mean_norm_count))) +
  #geom_errorbar(aes(ymin=log(max_norm_count), ymax=log(min_norm_count)), width=0.2, position = position_dodge(0.2)) +
  geom_boxplot(outlier.shape = NA, fill = "gray85", colour = "gray10") +
  geom_line(aes(group = gene_name), colour = "gray40") +
  geom_point(colour = "blue4", fill = "blue", alpha = 0.6, cex = 4, shape = 21)+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18))+
  xlab("Ecotype and treatment")+
  ylab("log of normalised gene counts")+
  scale_x_discrete(labels=c("cstC" = "Control", "cstS" = "Salt", "cstZ" = "Zinc",
                            "minC" = "Control", "minS" = "Salt", "minZ" = "Zinc"))+
  #facet_wrap(~ecotype, ncol = 2, scales="free_x", strip.position = "top")
  facet_wrap2(~ecotype, ncol = 2, scales = "free_x", strip.position = "top",
              labeller = labeller(ecotype = ecotype.labs), strip = strip)
Cue_transfer_simple_lineboxplot2

#colours that are colour-blind friendly from particular palette, run to list them:
#palette.colors(palette = "Okabe-Ito")
#could try again with the safe palette)


###11.4.2 plot Cooption ----

#factor order with control first:
cooption_mean_simple_lineboxplot2 <- ggplot(cooption_table2, aes(x = eco_treat, y = log(mean_norm_count))) +
  geom_boxplot(outlier.shape = NA, fill = "gray85") +
  geom_line(aes(group = gene_name), colour = "gray40") +
  geom_point(colour = "darkred", fill = "red", alpha = 0.6, cex = 4, shape = 21)+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18))+
  xlab("Ecotype and treatment")+
  ylab("log of normalised gene counts")+
  scale_x_discrete(labels=c("cstC" = "Control", "cstS" = "Salt",  "cstZ" = "Zinc",
                            "minC" = "Control", "minS" = "Salt",  "minZ" = "Zinc"))+
  #facet_wrap(~ecotype, ncol = 2, scales="free_x", strip.position = "top")
  facet_wrap2(~ecotype, ncol = 2, scales = "free_x", strip.position = "top",
              labeller = labeller(ecotype = ecotype.labs), strip = strip)
cooption_mean_simple_lineboxplot2


###11.4.3 plot the graph for the preadaptive_gene ----

#strict sense preadaptive plasticity (1 gene)
preadapt_simple_lineboxplot2 <- ggplot(preadapt_table2, aes(x = eco_treat, y = log(mean_norm_count))) +
  #geom_errorbar(aes(ymin=log(max_norm_count), ymax=log(min_norm_count)), width=0.2, position = position_dodge(0.2)) +
  geom_boxplot(outlier.shape = NA, fill = "gray85", colour = "gray10") +
  geom_line(aes(group = gene_name), colour = "gray40") +
  geom_point(colour = "blue4", fill = "blue", alpha = 0.6, cex = 4, shape = 21)+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18))+
  xlab("Ecotype and treatment")+
  ylab("log of normalised gene counts")+
  scale_x_discrete(labels=c("cstC" = "Control", "cstS" = "Salt", "cstZ" = "Zinc",
                            "minC" = "Control", "minS" = "Salt", "minZ" = "Zinc"))+
  #facet_wrap(~ecotype, ncol = 2, scales="free_x", strip.position = "top")
  facet_wrap2(~ecotype, ncol = 2, scales = "free_x", strip.position = "top",
              labeller = labeller(ecotype = ecotype.labs), strip = strip)
preadapt_simple_lineboxplot2

#relaxed sense preadaptive plasticity (allows for random plasticity value in salt in mines)
Lpreadapt_simple_lineboxplot2 <- ggplot(Lpreadapt_table2, aes(x = eco_treat, y = log(mean_norm_count))) +
  #geom_errorbar(aes(ymin=log(max_norm_count), ymax=log(min_norm_count)), width=0.2, position = position_dodge(0.2)) +
  geom_boxplot(outlier.shape = NA, fill = "gray85", colour = "gray10") +
  geom_line(aes(group = gene_name), colour = "gray40") +
  geom_point(colour = "blue4", fill = "blue", alpha = 0.6, cex = 4, shape = 21)+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18))+
  xlab("Ecotype and treatment")+
  ylab("log of normalised gene counts")+
  scale_x_discrete(labels=c("cstC" = "Control", "cstS" = "Salt", "cstZ" = "Zinc",
                            "minC" = "Control", "minS" = "Salt", "minZ" = "Zinc"))+
  #facet_wrap(~ecotype, ncol = 2, scales="free_x", strip.position = "top")
  facet_wrap2(~ecotype, ncol = 2, scales = "free_x", strip.position = "top",
              labeller = labeller(ecotype = ecotype.labs), strip = strip)
Lpreadapt_simple_lineboxplot2



###11.4.3 export plots to pdfs ----

# 
# # #export with control over dimensions:
# # 
setFplot_page(page = "a4", margins = "normal", units = "tw",pt = 20, w2h = 1.8, reset = FALSE)
# # 
# 
# pdf_fit(file = "cue_transfer_lineboxplot1_16_11_23.pdf", pt =16, width = 1.8, w2h = 0.8)
# cue_transfer_mean_lineboxplot1
# dev.off()
# 
# pdf_fit(file = "cooption_lineboxplot1_16_11_23.pdf", pt =16, width = 1.8, w2h = 0.8)
# cooption_mean_lineboxplot1
# dev.off()
# 

# #cue transfer simpler graph for presentations: 3rd year review
# pdf_fit(file = "Cue_transfer_simple_lineboxplot1_30_01_24.pdf", pt = 26, width = 1.8, w2h = 0.8)
# Cue_transfer_simple_lineboxplot1
# dev.off()

pdf_fit(file = "Cue_transfer_simple_lineboxplot2_09_04_24.pdf", pt = 26, width = 1.8, w2h = 0.8)
Cue_transfer_simple_lineboxplot2
dev.off()
# 
# #cooption simpler graph for presentations: 3rd year review
# pdf_fit(file = "cooption_mean_simple_lineboxplot_30_01_24.pdf", pt = 26, width = 1.8, w2h = 0.8)
# cooption_mean_simple_lineboxplot
# dev.off()
# 
pdf_fit(file = "cooption_mean_simple_lineboxplot2_09_04_24.pdf", pt = 26, width = 1.8, w2h = 0.8)
cooption_mean_simple_lineboxplot2
dev.off()
#
# pdf_fit(file = "both_simple_lineboxplot_30_01_23.pdf", pt = 26, width = 1.8, w2h = 0.8)
# both_simple_lineboxplot1
# fit.off()

pdf_fit(file = "preadapt_simple_lineboxplot2_30_07_24.pdf", pt = 26, width = 1.8, w2h = 0.8)
preadapt_simple_lineboxplot2
dev.off()

################################################################################-
#12.0 Functional analyses - annotations and GO ontology enrichment ----
################################################################################-

##12.1 sort out annotations for the total no. genes in the CC_all set used across the experiment ----

dim(CC_all)
#generate names for the set of genes that have model instead of TU as a title.
CC_all_model_names <- gsub(".TU.", ".model.", c(CC_all$Names))

#create annotation set that consists of only those within the 23,093 genes consistently expressed across experiments
CC_all_annot <- fun_annot[which(fun_annot$mRNA %in% CC_all_model_names),]

#only 19134 genes with annotations out of 23093
dim(CC_all_annot)[1]

#number of genes not in fun_annot or not matching a gene TU name out of the CC_all set
dim(CC_all)[1] - dim(CC_all_annot)[1] #3959
#names of genes not in the annotation file as far as I know.
no_annotation_gene_names <- CC_all_model_names[CC_all_model_names %notin% CC_all_annot$mRNA]
length(no_annotation_gene_names) #3959 genes - will use as a check against full file.

##########################################################-
##12.2 Generate annotation files and gene name lists ----
##########################################################-

#generate names for the set of genes that have model instead of TU as a title.

#overlapping ancestral and descendent salt plasticity:
ASP_DSP_same_model_names <- gsub(".TU.", ".model.", c(ASP_DSP_same$Names))
#fixed evolved changes:
ECC_model_names <- gsub(".TU.", ".model.", c(ECC_noDZP$Names))
#evolve changes to zinc
ECZ_model_names <- gsub(".TU.", ".model.", c(ECZ_same$Names))
#ASP
ASP_model_names <- gsub(".TU.", ".model.", c(ASP_same$Names))
#DSP
DSP_model_names <- gsub(".TU.", ".model.", c(DSP_same$Names))
#AZP
AZP_model_names <- gsub(".TU.", ".model.", c(AZP_same$Names))
#DZP
DZP_model_names <- gsub(".TU.", ".model.", c(DZP_same$Names))
#EDZP
EDZP_model_names <- gsub(".TU.", ".model.", c(EDZP$Names))

#plasticity to salt in coasts having new zinc cue total set
cue_transfer_model_names <- gsub(".TU.", ".model.", c(ASP_EDZP_SZcst$Names))
#harmonisation gene set
harmonisation_model_names <- gsub(".TU.", ".model.", c(ASP_EDZP_SZcst_DSP$Names))
#plasticity switching 5 gene set:
Pl_switch_model_names <- gsub(".TU.", ".model.", c(ASP_EDZP_SZcst_no_DSP$Names))

#evolved change to zinc and plastic to salt in coasts total set 
cooption_model_names <- gsub(".TU.", ".model.", c(ASP_ECC_noDZP_SZcst$Names))
#evolved subsets full cooption - no salt plasticity
cooption_noDSP_model_names <- gsub(".TU.", ".model.", c(ASP_ECC_noDZP_SZcst_noDSP$Names))
#partial cooption - salt plasticity remains (not sure this group is informative)
cooption_DSP_model_names <- gsub(".TU.", ".model.", c(ASP_ECC_noDZP_SZcst_DSP$Names))

#annotation set that consists of subset of the 2 genegroups of interest subsetted from CC_all_annot (19,304)

ASP_DSP_same_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% ASP_DSP_same_model_names),]
#ECC with no DZP (Evo change in controls that is constitutive)
ECC_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% ECC_model_names),]
#dim(ECC_annots)[1] #73
#ECZ all evolutionary changes in zinc expression levels
ECZ_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% ECZ_model_names),]
#dim(ECZ_annots)[1] #8136
#ASP
ASP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% ASP_model_names),]
#dim(ASP_annots)[1] #876
view(ASP_annots)
#DSP
DSP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% DSP_model_names),]
#dim(DSP_annots)[1] #144
#AZP
AZP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% AZP_model_names),]
#dim(AZP_annots)[1] #9640
#DZP
DZP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% DZP_model_names),]
#dim(DZP_annots)[1] #125
#EDZP
EDZP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% EDZP_model_names),]
#dim(EDZP_annots)[1] #81

#COOPTION ASP_ECC_noDZP_SZcst
cooption_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% cooption_model_names),]
dim(cooption_annots)[1] #36
#full cooption
cooption_noDSP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% cooption_noDSP_model_names),]
dim(cooption_noDSP_annots)[1] #22
#partial cooption
cooption_DSP_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% cooption_DSP_model_names),]
dim(cooption_DSP_annots)[1] #10

#CUE TRANSFER (past plasticity to salt also present in response to zinc)
cue_transfer_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% cue_transfer_model_names),]
dim(cue_transfer_annots)[1] #26
#harmonisation
harmonisation_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% harmonisation_model_names),]
dim(harmonisation_annots)[1] #13
#plast switching
Pl_switch_annots <- CC_all_annot[which(CC_all_annot$mRNA %in% Pl_switch_model_names),]
dim(Pl_switch_annots)[1] #5

##export full files:
# write.csv(ASP_DSP_same_annots, "ASP_DSP_same_annots_26_05_23.csv")
# write.csv(ECC_annots, "ECC_annots_23_05_23.csv")
# write.csv(ECZ_annots, "ECZ_annots_23_05_23.csv")
#write.csv(ASP_annots, "ASP_annots_23_05_23.csv")
#write.csv(DSP_annots, "DSP_annots_23_05_23.csv")
# write.csv(AZP_annots, "AZP_annots_23_05_23.csv")
# write.csv(DZP_annots, "DZP_annots_23_05_23.csv")
# write.csv(EDZP_annots, "EDZP_annots_23_05_23.csv")

## patterns for hypotheses:
write.csv(cooption_annots, "cooption__annots_19_01_24.csv")
write.csv(cue_transfer_annots, "Cue_transfer_annots_19_01_24.csv")
## 
# #extra ones I added later:
# write.csv(cooption_noDSP_annots, "cooption_noDSP_annots_23_06_23.csv")
# write.csv(cooption_DSP_annots, "cooption_DSP_annots_23_06_23.csv")
# write.csv(harmonisation_annots, "harmonisation_annots_23_06_23.csv")
# write.csv(Pl_switch_annots, "Pl_switch_annots_23_10_23.csv")

ECC_annots_Names <- ECC_annots$mRNA
ECZ_annots_Names <- ECZ_annots$mRNA
ASP_annots_Names <- ASP_annots$mRNA
DSP_annots_Names <- DSP_annots$mRNA
AZP_annots_Names <- AZP_annots$mRNA
DZP_annots_Names <- DZP_annots$mRNA
EDZP_annots_Names <- EDZP_annots$mRNA
ASP_DSP_same_annots_Names <- ASP_DSP_same_annots$mRNA
cooption_annots_Names <- cooption_annots$mRNA
cooption_noDSP_annots_Names <- cooption_noDSP_annots$mRNA
cooption_DSP_annots_Names <- cooption_DSP_annots$mRNA
cue_transfer_annots_Names <- cue_transfer_annots$mRNA
harmonisation_annots_Names <- harmonisation_annots$mRNA
Pl_switch_annots_Names <- Pl_switch_annots$mRNA

#export files that are just names (maybe don't need this):
# write.csv(ASP_DSP_same_annots$mRNA, "ASP_DSP_same_annots_Names_26_05_23.csv")
# write.csv(ECC_annots_Names, "ECC_annots_Names_23_05_23")
# write.csv(ECZ_annots_Names, "ECZ_annots_Names_23_05_23")
# write.csv(ASP_annots_Names, "ASP_annots_Names_23_05_23")
# write.csv(DSP_annots_Names, "DSP_annots_Names_23_05_23")
# write.csv(AZP_annots_Names, "AZP_annots_Names_23_05_23")
# write.csv(DZP_annots_Names, "DZP_annots_Names_23_05_23")
# write.csv(EDZP_annots_Names, "EDZP_annots_Names_23_05_23")

# write.csv(cooption_annots_Names, "cooption_annots_Names_19_01_24")
# write.csv(cue_transfer_annots_Names, "cue_transfer_annots_Names_19_01_24")


#######################################-
##12.3 TOPGO analysis ----
#######################################-

###12.3.1 format fun_annot into GeneIdtoGO format and export to then be able to have it import in the correct format ---- 

dim(fun_annot)
dim(CC_all_annot)

Geneidtogo <- CC_all_annot %>% dplyr::select(., mRNA, GOs)
dim(Geneidtogo)[1]
str(Geneidtogo)

readr::write_tsv(Geneidtogo, file = "geneIDtoGO_23_05_23.map", col_names = F)

###12.3.2 Load topGO and geneIDtoGO map file ----

#only load topgo here as has select function which intereferes with the dplyr select in previous lines of code.
library(topGO) 
packageVersion("topGO")

#set up geneIDtoGO file!
geneID2GO <-readMappings(file = "geneIDtoGO_23_05_23.map")
length(geneID2GO)
geneNames <- names(geneID2GO)
head(geneNames)
length(geneNames)

###12.3.3 Define genelist tables ----

#genelists are those sets of genes you find interesting.

# geneList_ECC <- factor(as.integer(geneNames %in% ECC_annots_Names))
# names(geneList_ECC) <- geneNames
# str(geneList_ECC)
# 
# geneList_ECZ <- factor(as.integer(geneNames %in% ECZ_annots_Names))
# names(geneList_ECZ) <- geneNames
# str(geneList_ECZ)

geneList_ASP <- factor(as.integer(geneNames %in% ASP_annots_Names))
names(geneList_ASP) <- geneNames
str(geneList_ASP)

geneList_DSP <- factor(as.integer(geneNames %in% DSP_annots_Names))
names(geneList_DSP) <- geneNames
str(geneList_DSP)

# geneList_ASP_DSP <- factor(as.integer(geneNames %in% ASP_DSP_same_annots_Names))
# names(geneList_ASP_DSP) <- geneNames
# str(geneList_ASP_DSP)

# geneList_AZP <- factor(as.integer(geneNames %in% AZP_annots_Names))
# names(geneList_AZP) <- geneNames
# str(geneList_AZP)
# 
# geneList_DZP <- factor(as.integer(geneNames %in% DZP_annots_Names))
# names(geneList_DZP) <- geneNames
# str(geneList_DZP)
# 
# geneList_EDZP <- factor(as.integer(geneNames %in% EDZP_annots_Names))
# names(geneList_EDZP) <- geneNames
# str(geneList_EDZP)

geneList_cue_transfer <- factor(as.integer(geneNames %in% cue_transfer_annots$mRNA))
names(geneList_cue_transfer) <- geneNames
str(geneList_cue_transfer)

# geneList_harmonisation <- factor(as.integer(geneNames %in% harmonisation_annots_Names))
# names(geneList_harmonisation) <- geneNames
# 
# geneList_Pl_switch <- factor(as.integer(geneNames %in% Pl_switch_annots_Names))
# names(geneList_Pl_switch) <- geneNames

geneList_cooption <- factor(as.integer(geneNames %in% cooption_annots_Names))
names(geneList_cooption) <- geneNames
str(geneList_cooption)

# geneList_cooption_noDSP <- factor(as.integer(geneNames %in% cooption_noDSP_annots_Names))
# names(geneList_cooption_noDSP) <- geneNames
# str(geneList_cooption_noDSP)
# 
# geneList_cooption_DSP <- factor(as.integer(geneNames %in% cooption_DSP_annots_Names))
# names(geneList_cooption_DSP) <- geneNames
# str(geneList_cooption_DSP)


###12.3.4 topGO data set up a new one for each topGO ontology analysis for each of my main gene sets of interest ----

#onotology options are BP = biological process, MF  = molecular function, CC = cellular component 
# GOdata_ECC <- new("topGOdata", ontology = "BP", allGenes = geneList_ECC, 
#                       annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_ECC)
# 
# GOdata_ECZ <- new("topGOdata", ontology = "BP", allGenes = geneList_ECZ, 
#                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_ECZ)

GOdata_ASP <- new("topGOdata", ontology = "BP", allGenes = geneList_ASP, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
numGenes(GOdata_ASP)

GOdata_DSP <- new("topGOdata", ontology = "BP", allGenes = geneList_DSP, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
numGenes(GOdata_DSP)

# GOdata_ASP_DSP <- new("topGOdata", ontology = "BP", allGenes = geneList_ASP_DSP, 
#                       annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_ASP_DSP)
# 
# GOdata_AZP <- new("topGOdata", ontology = "BP", allGenes = geneList_AZP, 
#                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_AZP)
# 
# GOdata_DZP <- new("topGOdata", ontology = "BP", allGenes = geneList_DZP, 
#                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_DZP)
# 
# GOdata_EDZP <- new("topGOdata", ontology = "BP", allGenes = geneList_EDZP, 
#                    annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_EDZP)

GOdata_cue_transfer <- new("topGOdata", ontology = "BP", allGenes = geneList_cue_transfer, 
                             annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
numGenes(GOdata_cue_transfer)

# GOdata_harmonisation <- new("topGOdata", ontology = "BP", allGenes = geneList_harmonisation, 
#                                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_harmonisation)

GOdata_cooption <- new("topGOdata", ontology = "BP", allGenes = geneList_cooption, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
numGenes(GOdata_cooption)

# GOdata_cooption_noDSP <- new("topGOdata", ontology = "BP", allGenes = geneList_cooption_noDSP, 
#                             annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_cooption_noDSP)
# 
# GOdata_cooption_DSP <- new("topGOdata", ontology = "BP", allGenes = geneList_cooption_DSP, 
#                             annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
# numGenes(GOdata_cooption_DSP)


#run the topGO ontology test:
# fisher_weight_ECC <- runTest(GOdata_ECC, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_ECC))
# geneData(fisher_weight_ECC)
# 
# fisher_weight_ECZ <- runTest(GOdata_ECZ, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_ECZ))
# geneData(fisher_weight_ECZ)

fisher_weight_ASP <- runTest(GOdata_ASP, algorithm = "weight", statistic = "fisher")
head(score(fisher_weight_ASP))
geneData(fisher_weight_ASP)

fisher_weight_DSP <- runTest(GOdata_DSP, algorithm = "weight", statistic = "fisher")
head(score(fisher_weight_DSP))
geneData(fisher_weight_DSP)

# fisher_weight_ASP_DSP <- runTest(GOdata_ASP_DSP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_ASP_DSP))
# geneData(fisher_weight_ASP_DSP)
# 
# fisher_weight_AZP <- runTest(GOdata_AZP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_AZP))
# geneData(fisher_weight_AZP)
# 
# fisher_weight_DZP <- runTest(GOdata_DZP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_DZP))
# geneData(fisher_weight_DZP)

# fisher_weight_EDZP <- runTest(GOdata_EDZP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_EDZP))
# geneData(fisher_weight_EDZP)

fisher_weight_cue_transfer <- runTest(GOdata_cue_transfer, algorithm = "weight", statistic = "fisher")
head(score(fisher_weight_cue_transfer))
geneData(fisher_weight_cue_transfer)

# fisher_weight_harmonisation <- runTest(GOdata_harmonisation, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_harmonisation))
# geneData(fisher_weight_harmonisation)

fisher_weight_cooption <- runTest(GOdata_cooption, algorithm = "weight", statistic = "fisher")
head(score(fisher_weight_cooption))
geneData(fisher_weight_cooption)

# fisher_weight_cooption_noDSP <- runTest(GOdata_cooption_noDSP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_cooption_noDSP))
# geneData(fisher_weight_cooption_noDSP)
# 
# fisher_weight_cooption_DSP <- runTest(GOdata_cooption_DSP, algorithm = "weight", statistic = "fisher")
# head(score(fisher_weight_cooption_DSP))
# geneData(fisher_weight_cooption_DSP)


###12.3.5 TOPGO look at results - generate results tables and export significant results ----

#gene tables
# res_fisher_weight_ECC <- GenTable(GOdata_ECC, weight_fisher_P = fisher_weight_ECC, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_ECC)
# 
# res_fisher_weight_ECZ <- GenTable(GOdata_ECZ, weight_fisher_P = fisher_weight_ECZ, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_ECZ)

res_fisher_weight_ASP <- GenTable(GOdata_ASP, weight_fisher_P = fisher_weight_ASP, topNodes = 50, numChar = 1000)
#view(res_fisher_weight_ASP)

res_fisher_weight_DSP <- GenTable(GOdata_DSP, weight_fisher_P = fisher_weight_DSP, topNodes = 50, numChar = 1000)
#view(res_fisher_weight_DSP)

# res_fisher_weight_ASP_DSP <- GenTable(GOdata_ASP_DSP, weight_fisher_P = fisher_weight_ASP_DSP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_ASP)
# 
# res_fisher_weight_AZP <- GenTable(GOdata_AZP, weight_fisher_P = fisher_weight_AZP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_AZP)
# 
# res_fisher_weight_DZP <- GenTable(GOdata_DZP, weight_fisher_P = fisher_weight_DZP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_DZP)

# res_fisher_weight_EDZP <- GenTable(GOdata_EDZP, weight_fisher_P = fisher_weight_EDZP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_EDZP)

res_fisher_weight_cue_transfer <- GenTable(GOdata_cue_transfer, weight_fisher_P = fisher_weight_cue_transfer, topNodes = 50, numChar = 1000)
#view(res_fisher_weight_cue_transfer)

# res_fisher_weight_harmonisation <- GenTable(GOdata_harmonisation, weight_fisher_P = fisher_weight_harmonisation, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_harmonisation)

res_fisher_weight_cooption <- GenTable(GOdata_cooption, weight_fisher_P = fisher_weight_cooption, topNodes = 50, numChar = 1000)
#view(res_fisher_cooption)

# res_fisher_weight_cooption_noDSP <- GenTable(GOdata_cooption_noDSP, weight_fisher_P = fisher_weight_cooption_noDSP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_cooption_noDSP)
# 
# res_fisher_weight_cooption_DSP <- GenTable(GOdata_cooption_DSP, weight_fisher_P = fisher_weight_cooption_DSP, topNodes = 50, numChar = 1000)
# #view(res_fisher_weight_cooption_DSP)


##significant GO terms for the gene sets:
# sign_fisher_weight_ECC <- subset(res_fisher_weight_ECC, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_ECZ <- subset(res_fisher_weight_ECZ, as.numeric(weight_fisher_P) < 0.05)
sign_fisher_weight_ASP <- subset(res_fisher_weight_ASP, as.numeric(weight_fisher_P) < 0.05)
sign_fisher_weight_DSP <- subset(res_fisher_weight_DSP, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_ASP_DSP <- subset(res_fisher_weight_ASP_DSP, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_AZP <- subset(res_fisher_weight_AZP, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_DZP <- subset(res_fisher_weight_DZP, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_EDZP <- subset(res_fisher_weight_EDZP, as.numeric(weight_fisher_P) < 0.05)

sign_fisher_weight_cue_transfer <- subset(res_fisher_weight_cue_transfer, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_harmonisation <- subset(res_fisher_weight_harmonisation, as.numeric(weight_fisher_P) < 0.05)

sign_fisher_weight_cooption <- subset(res_fisher_weight_cooption, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_cooption_noDSP <- subset(res_fisher_weight_cooption_noDSP, as.numeric(weight_fisher_P) < 0.05)
# sign_fisher_weight_cooption_DSP <- subset(res_fisher_weight_cooption_DSP, as.numeric(weight_fisher_P) < 0.05)


##write files for the significant fisher results so I have the tables saved as files:
# write.csv(sign_fisher_weight_ECC, file = "topGO_sign_fisher_weight_ECC_23_05_23.csv")
# write.csv(sign_fisher_weight_ECZ, file = "topGO_sign_fisher_weight_ECZ_23_05_23.csv")
# write.csv(sign_fisher_weight_ASP, file = "topGO_sign_fisher_weight_ASP_23_05_23.csv")
# write.csv(sign_fisher_weight_DSP, file = "topGO_sign_fisher_weight_DSP_23_05_23.csv")
# write.csv(sign_fisher_weight_AZP, file = "topGO_sign_fisher_weight_AZP_23_05_23.csv")
# write.csv(sign_fisher_weight_DZP, file = "topGO_sign_fisher_weight_DZP_23_05_23.csv")
# write.csv(sign_fisher_weight_EDZP, file = "topGO_sign_fisher_weight_EDZP_23_05_23.csv")
##extra single enrichments added later on:
#write.csv(sign_fisher_weight_ASP_DSP, "topGO_sign_fisher_weight_ASP_DSP_23_06_23.csv")
#
#write.csv(sign_fisher_weight_cue_transfer, file = "topGO_sign_fisher_weight_cue_transfer_19_01_24.csv")
#write.csv(sign_fisher_weight_harmonisation, file = "topGO_sign_fisher_weight_harmonisation_19_01_24.csv")
#
#write.csv(sign_fisher_weight_cooption, file = "topGO_sign_fisher_weight_cooption_19_01_24.csv")
#write.csv(sign_fisher_weight_cooption_noDSP, "topGO_sign_fisher_weight_cooption_noDSP_19_01_24.csv")
#write.csv(sign_fisher_weight_cooption_DSP, "topGO_sign_fisher_weight_cooption_DSP_19_01_24.csv")


#look for overlaps between enriched functions for gene groups of interest: 
ASP_DSP_GO_overlap <- sign_fisher_weight_DSP[sign_fisher_weight_DSP$GO.ID %in% sign_fisher_weight_ASP$GO.ID,]
ASP_DSP_GO_overlap_terms <- dplyr::select(ASP_DSP_GO_overlap, GO.ID, Term)

# #overlaps write:
# write.csv(ASP_DSP_GO_overlap_terms, file = "ASP_DSP_shared_GO_terms_26_05_23.csv")

############################################################################################-
#13.0 Go results visualisation for aim 1 comparing salt response in mine and coast ----
###########################################################################################-

###plot GO graphs using data structure similar to the package GOplot but made from topGO results

# #install.packages("GOplot")
# 
# library(GOplot)
# 
# data(EC)
# head(EC$david)
# str(EC$david)
# str(EC$genelist)
# test_circ <- circle_dat(terms = EC$david, EC$genelist)
# str(test_circ)

##13.1  try to make LFC tables for genes to see if this is useful at all ----

# ASP_same$ASP_logfc <- rowMeans(ASP_same[,c(2,4)])
# colnames(ASP_same)
# keepcols <- c("Names", "ASP_logfc", "padj_CSbd_sig")
# hist(ASP_same$ASP_logfc)
# boxplot(ASP_same$ASP_logfc)
# 
# DSP_same$DSP_logfc <- rowMeans(DSP_same[,c(2,4)])
# hist(DSP_same$DSP_logfc)
# boxplot(DSP_same$DSP_logfc)
# 
# ASP_DSP_same$ASP_lfc <- rowMeans(ASP_DSP_same[,c(2,4)])
# ASP_DSP_same$DSP_lfc <- rowMeans(ASP_DSP_same[,c(7,9)])
# 
# boxplot(ASP_DSP_same$ASP_lfc, ASP_DSP_same$DSP_lfc)
# boxplot(ASP_same$ASP_logfc, DSP_same$DSP_logfc)


##13.2 GO analysis create ASP and DSP tables for bubble plotting ----

### 13.2.1 create tables of mean ASP and DSP logFC for sign. DE genes ----
#ASP

keepcols <- c("Names", "ASP_logfc", "padj_CSbd_sig")
ASP_same$ASP_logfc <- rowMeans(ASP_same[, c(2,4)])
ASP_goplot <- ASP_same[, ..keepcols]
ASP_goplot$Genes <- ASP_model_names
head(ASP_goplot)
str(ASP_goplot)

ASP_goplot <- ASP_goplot[which(ASP_goplot$Genes %in% ASP_annots$mRNA),]

#DSP
keepcols1 <- c("Names", "DSP_logfc", "padj_CSgr_sig")
DSP_same$DSP_logfc <- rowMeans(DSP_same[, c(2,4)])
DSP_goplot <- DSP_same[, ..keepcols1]
DSP_goplot$Genes <- DSP_model_names
DSP_goplot <- DSP_goplot[which(DSP_goplot$Genes %in% DSP_annots$mRNA),]

###13.2.2 get GO gene names mapped to topGO terms ----

#get go gene names from the significant set in a list:
ASP_sign_GOID <- c(sign_fisher_weight_ASP$GO.ID)
sign_ASP_GOs <- genesInTerm(GOdata_ASP, ASP_sign_GOID)
str(sign_ASP_GOs)

DSP_sign_GOID <- c(sign_fisher_weight_DSP$GO.ID)
sign_DSP_GOs <- genesInTerm(GOdata_DSP, DSP_sign_GOID)
str(sign_DSP_GOs)

#create table of GOs and gene names using function list to data frame
gene_to_GO_ASP <- list2df_dt(sign_ASP_GOs, "ID", "Genes")
gene_to_GO_DSP <- list2df_dt(sign_DSP_GOs, "ID", "Genes")

# #check details of the tables I have made:
# head(ASP_goplot)
# dim(ASP_goplot)
# head(gene_to_GO_ASP)
# dim(gene_to_GO_ASP)
# 
# #checking content of the table:
# ASP_goplot[which(ASP_goplot$Genes %in% gene_to_GO_ASP$Genes),]
# gene_to_GO_ASP[which(gene_to_GO_ASP$Genes %in% ASP_goplot$Genes),]

###13.2.3 combine gene+GOID tables with topGO p value results ----

#merge the lfc results plus gene names table with gene+GOIDs

ASP_bubble_table <- merge(gene_to_GO_ASP, ASP_goplot, by = "Genes")
DSP_bubble_table <- merge(gene_to_GO_DSP, DSP_goplot, by = "Genes")

#create final combined table of LFCs, gene names and topGO fisher table
sign_fisher_weight_ASP <- rename(sign_fisher_weight_ASP, "ID" = "GO.ID")
sign_fisher_weight_DSP <- rename(sign_fisher_weight_DSP, "ID" = "GO.ID")

ASP_bubble_table_final <- merge(ASP_bubble_table, sign_fisher_weight_ASP, by = "ID")
DSP_bubble_table_final <- merge(DSP_bubble_table, sign_fisher_weight_DSP, by = "ID")

#could add column for direction, but not longer needed.
# ASP_bubble_table_final$lfc_dir <- ifelse(ASP_bubble_table_final$ASP_logfc > 0, "+", "-")
# view(ASP_bubble_table_final)

# ASP_bubble_table_mean <- ASP_bubble_table_final %>% group_by(ID) %>%
#   summarise_at(vars(ASP_logfc), list(GO_mean_logfc = mean))

### practice code for developing the zscore for the dataset:

# # formula:
# #z = # up - # down / count^-1
# 
# fold <- c(1,-1,1,-2,-4,-3)
# ID <- c("a", "a", "a", "a", "b", "b")
# gene <- c("lemon", "orangefruit", "banana", "blueberry", "blueberry1", "blueberry2")
# 
# zscore_test <- data.frame(cbind(fold, ID, gene))
# 
# # ASP_bubble_plot_table1 <- ASP_bubble_table %>% group_by(ID) %>%
# # summarise_at(vars(ASP_logfc), list(zscore = zscorefun)
# 
# #manually make zscore for one dataset
# (length(which(zscore_test$fold > 0)) - length(which(zscore_test$fold < 0))) / length(zscore_test$fold)
# 
# #length(which(.x > 0)) - length(which(.x < 0))) / length(.x)
# #test the thing as a lambda function inside summarise:
# zscore_test %>%
#   group_by(ID) %>%
#   summarise(across(c(fold), list(zscore = ~ (length(which(.x > 0)) - length(which(.x < 0)))/ length(.x))))


###13.2.6 Calculate zscore and add to new tables ready for plotting: ----

# # formula:
# #z = # up - # down / count^-1

ASP_bubble_table_zscore <- ASP_bubble_table_final %>% 
  group_by(ID) %>%
  summarise(across(c(ASP_logfc), list(zscore = ~ (length(which(.x > 0)) - length(which(.x < 0)))/ length(.x))))
ASP_bubble_table_zscore  <- rename(ASP_bubble_table_zscore, "zscore" = "ASP_logfc_zscore")

DSP_bubble_table_zscore <- DSP_bubble_table_final %>% 
  group_by(ID) %>%
  summarise(across(c(DSP_logfc), list(zscore = ~ (length(which(.x > 0)) - length(which(.x < 0)))/ length(.x))))
#rename zscore to just be zscore:
DSP_bubble_table_zscore  <- rename(DSP_bubble_table_zscore, "zscore" = "DSP_logfc_zscore")

##merge the tables with zscore with the p value tables from topGO: 

ASP_bubble_table_zscore_plot <- merge(sign_fisher_weight_ASP, ASP_bubble_table_zscore, by = "ID")
#view(ASP_bubble_table_zscore_plot)

DSP_bubble_table_zscore_plot <- merge(sign_fisher_weight_DSP, DSP_bubble_table_zscore, by = "ID")

### 13.2.7 add to tables the columns of broader function categories I defined: ----

#import tables:
ASP_broader_func <- data.frame(read.csv("ASP_GOID_broader_function.csv"))
DSP_broader_func <- data.frame(read.csv("DSP_GOID_broader_function.csv"))

ASP_broader_func1 <- data.frame(read.csv("ASP_GOs_7_functions.csv"))
DSP_broader_func1 <- data.frame(read.csv("DSP_GOs_7_functions.csv"))

#merge tables together via the ID column:
ASP_bubble_table_zscore_plot1 <- merge(ASP_bubble_table_zscore_plot, ASP_broader_func1, by = "ID")
DSP_bubble_table_zscore_plot1 <- merge(DSP_bubble_table_zscore_plot, DSP_broader_func1, by = "ID")

str(ASP_bubble_table_zscore_plot1)

ASP_bubble_table_zscore_plot1$ID <- factor(ASP_bubble_table_zscore_plot1$ID)
range(ASP_bubble_table_zscore_plot1$Significant)
ASP_bubble_table_zscore_plot1$weight_fisher_P <- as.numeric(ASP_bubble_table_zscore_plot1$weight_fisher_P)
ASP_bubble_table_zscore_plot1$broader_function <- factor(ASP_bubble_table_zscore_plot1$broader_function)
ASP_rep <- rep(c("ASP"), 50)
ASP_bubble_table_zscore_plot1$gene_set <- ASP_rep

DSP_bubble_table_zscore_plot1$ID <- factor(DSP_bubble_table_zscore_plot1$ID)
range(DSP_bubble_table_zscore_plot1$Significant)
DSP_bubble_table_zscore_plot1$weight_fisher_P <- as.numeric(DSP_bubble_table_zscore_plot1$weight_fisher_P)
DSP_bubble_table_zscore_plot1$broader_function <- factor(DSP_bubble_table_zscore_plot1$broader_function)
DSP_rep <- rep(c("DSP"), 34)
DSP_bubble_table_zscore_plot1$gene_set <- DSP_rep


##make combined table:
#view(ASP_bubble_table_zscore_plot1)
#view(DSP_bubble_table_zscore_plot1)
dim(ASP_bubble_table_zscore_plot1)
dim(DSP_bubble_table_zscore_plot1)
names(ASP_bubble_table_zscore_plot1) == names(DSP_bubble_table_zscore_plot1)
 
salt_bubble_table <- rbind(ASP_bubble_table_zscore_plot1, DSP_bubble_table_zscore_plot1)

#view(salt_bubble_table)

#create labelling selection from the dataset:
# subset(salt_bubble_table, ID %in% c("GO:0006970", "GO:0010200"))

shared_GOs <- intersect(ASP_bubble_table_zscore_plot1$ID, DSP_bubble_table_zscore_plot1$ID)

length(shared_GOs)

#create tables for export of shared and unique Go terms:
salt_bubble_table$shared_GO <- ifelse(salt_bubble_table$ID %in% shared_GOs, "yes", "no")
#view(salt_bubble_table)
dim(salt_bubble_table)

shared_bubble_table <- salt_bubble_table[salt_bubble_table$shared_GO == "yes",]
dim(shared_bubble_table)
shared_bubble_table_ASP <- shared_bubble_table[shared_bubble_table$gene_set == "ASP" ,]
shared_bubble_table_DSP <- shared_bubble_table[shared_bubble_table$gene_set == "DSP" ,]
shared_bubble_table_wide <- merge(shared_bubble_table_ASP, shared_bubble_table_DSP, by = "ID", suffixes = c("_ASP", "_DSP"))
#view(shared_bubble_table_wide)

shared_bubble_table_wide <- arrange(shared_bubble_table_wide, weight_fisher_P_ASP)

unique_bubble_table <- salt_bubble_table[salt_bubble_table$shared_GO == "no",]
dim(unique_bubble_table)
ASP_only_bubble_table <- unique_bubble_table[unique_bubble_table$gene_set == "ASP",]
dim(ASP_only_bubble_table)
ASP_only_bubble_table <- arrange(ASP_only_bubble_table, weight_fisher_P)
DSP_only_bubble_table <- unique_bubble_table[unique_bubble_table$gene_set == "DSP",] 
dim(DSP_only_bubble_table)  
DSP_only_bubble_table <- arrange(DSP_only_bubble_table, weight_fisher_P)


#write out the table in case it is useful :)
#write.csv(salt_bubble_table, "salt_bubble_table_28_02_24.csv")

# write.csv(shared_bubble_table_wide, "shared_bubble_table_wide_28_02_24_v1.csv")
# write.csv(ASP_only_bubble_table, "ASP_only_bubble_table_28_02_24_v1.csv")
# write.csv(DSP_only_bubble_table, "DSP_only_bubble_table_28_02_24_v1.csv")

#############-
##13.3 Plot the bubble plots ----
#############-

###13.3.2 stack the 2 sets together to make 1 graph together:

# Only colour strips in x-direction
strip <- strip_themed(background_x = elem_list_rect(fill = c(col_coast, col_mine)))
#create label for strip set to override current names:
gene_set.labs <- c("Coast salt DE genes", "Mine salt DE genes")
names(gene_set.labs) <- c("ASP", "DSP")

#in the scaling by size added a guide of none to get no key, but I wanted a key in the end:
#guide = "none"

salt_bubble_plot1 <- ggplot(salt_bubble_table, aes(x= zscore, 
                                                              y = -log10(weight_fisher_P), 
                                                              colour = broader_function, 
                                                              fill = broader_function,
                                                   label = shared_GO)) +
  geom_point(aes(size = Significant), shape = 21, alpha=0.6)+
  scale_size(range = c(range(salt_bubble_table$Significant*0.5)), name = "Gene number", breaks = c(5, 25, 50, 75))+
  #geom_text(data= subset(salt_bubble_table, ID %in% shared_GOs), aes(label = ID),
  #         show.legend = F, hjust = 0, vjust = 0.5, angle = 30, nudge_x = 0, check_overlap = T)+
  geom_point(data= subset(salt_bubble_table, ID %in% shared_GOs), shape = 23, size = 3, stroke = 1.5, show.legend = F, colour = "white")+
  xlab("zscore")+
  ylab("-log10 P value")+
  scale_colour_manual("Function", values = safe_colorblind_palette1)+
  scale_fill_manual("Function", values = safe_colorblind_palette1)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  coord_cartesian(ylim = c(1,9), xlim = c(-1.0, 1.2), expand = T)+
  theme(panel.background = element_rect(colour = "grey40", fill = NA), panel.grid = element_blank(),
        strip.background = element_rect(fill = c(col_coast, col_mine)), strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18), 
        legend.background = element_blank(), legend.key = element_blank(), legend.title = element_text(size = 22), legend.text = element_text(size = 18))+
  facet_wrap2(~gene_set, ncol = 2, scales = "fixed", strip.position = "top", strip = strip, 
              labeller = labeller(gene_set = gene_set.labs))
salt_bubble_plot1

# #export the plot:
# setFplot_page(page = "a4", margins = "normal", units = "tw",pt = 20, w2h = 1.8, reset = FALSE)
# 
# pdf_fit(file = "salt_bubble_plot1_IDlabels_24_01_24.pdf", pt =16, width = 2.8, w2h = 1.8)
# salt_bubble_plot1
# fit.off()

# #"tw" is scale to fraction of test size. or pw is page width
setFplot_page(page = "a4", margins = "normal", units = "tw", pt = 20, reset = FALSE, w2h = 2)
#
# version 16 is my preferred version so far: (with text size 20 and 18)
# pdf_fit(file = "salt_bubble_plot1_27_02_24_v16.pdf", pt =20, width = 2.8, w2h = 1.8)
# salt_bubble_plot1
# fit.off()

# pdf_fit(file = "salt_bubble_plot1_27_02_24_v17.pdf", pt =20, width = 3.4, w2h = 1)
# salt_bubble_plot1
# fit.off()

#try with pdf function. not nice text leading, although this might not be an issue:
# pdf(file = "salt_bubble_plot1_27_02_24_v18.pdf",width = 20, height = 15)
# salt_bubble_plot1
# dev.off()

#this one exported a file with text size ggplot edit to 22 and 18
# pdf_fit(file = "salt_bubble_plot1_27_02_24_v19.pdf", pt =20, width = 2.8, w2h = 1)
# salt_bubble_plot1
# fit.off()
# 
# pdf_fit(file = "salt_bubble_plot1_27_02_24_v20.pdf", pt =20, width = 3, w2h = 1)
# salt_bubble_plot1
# fit.off()

pdf_fit(file = "salt_bubble_plot1_27_02_24_v21.pdf", pt =20, width = 3, w2h = 1)
salt_bubble_plot1
fit.off()

pdf_fit(file = "salt_bubble_plot1_27_02_24_v23.pdf", pt =20, width = 3, w2h = 1)
salt_bubble_plot1
fit.off()


# # colour options for use in ggplot:
# scale_color_viridis(discrete = TRUE, option = "D")+ #contrast was a bit rubbish, default r colour palette better

#2nd version of plot with no labelling for the shared points:
salt_bubble_plot2 <- ggplot(salt_bubble_table, aes(x= zscore, 
                                                   y = -log10(weight_fisher_P), 
                                                   colour = broader_function, 
                                                   fill = broader_function,
                                                   label = shared_GO)) +
  geom_point(aes(size = Significant), shape = 21, alpha=0.6)+
  scale_size(range = c(range(salt_bubble_table$Significant*0.5)), name = "Gene number", breaks = c(5, 25, 50, 75))+
  #geom_text(data= subset(salt_bubble_table, ID %in% shared_GOs), 
  #        show.legend = F, hjust = 0.5, vjust = 0.5, angle = 30, nudge_x = 0.1, check_overlap = T)+
  #geom_point(data= subset(salt_bubble_table, ID %in% shared_GOs), shape = 23, size = 3, stroke = 0.8, show.legend = F, colour = "white")+
  xlab("zscore")+
  ylab("-log10 P value")+
  scale_colour_manual("Function", values = safe_colorblind_palette1)+
  scale_fill_manual("Function", values = safe_colorblind_palette1)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  coord_cartesian(ylim = c(1,9), xlim = c(-1.0, 1.2), expand = T)+
  theme(panel.background = element_rect(colour ="grey40", fill = NA), panel.grid = element_blank(),
        strip.background = element_rect(fill = c(col_coast, col_mine)), strip.text = element_text(size = 22),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18), 
        legend.background = element_blank(), legend.key = element_blank(), legend.title = element_text(size = 22), legend.text = element_text(size = 18))+
  facet_wrap2(~gene_set, ncol = 2, scales = "fixed", strip.position = "top", strip = strip, 
              labeller = labeller(gene_set = gene_set.labs))
#salt_bubble_plot2


# pdf_fit(file = "salt_bubble_plot2_27_02_24.pdf", pt =20, width = 3.0, w2h = 1.8)
# salt_bubble_plot2
# fit.off()


#3rd version of plot with outline to determine shared vs non shared. I'm not a fan.
salt_bubble_plot3 <- ggplot(salt_bubble_table, aes(x= zscore, 
                                                   y = -log10(weight_fisher_P), 
                                                   colour = broader_function, 
                                                   fill = broader_function,
                                                   label = shared_GO)) +
  geom_point(aes(size = Significant, colour = shared_GO), shape = 21, alpha=0.6, stroke = 1.5)+
  scale_size(range = c(range(salt_bubble_table$Significant*0.5)), name = "Gene number", breaks = c(5, 25, 50, 75))+
  #geom_text(data= subset(salt_bubble_table, ID %in% shared_GOs), 
  #        show.legend = F, hjust = 0.5, vjust = 0.5, angle = 30, nudge_x = 0.1, check_overlap = T)+
  #geom_point(data= subset(salt_bubble_table, ID %in% shared_GOs), shape = 23, size = 3, stroke = 0.8, show.legend = F, colour = "white")+
  xlab("zscore")+
  ylab("-log10 P value")+
  scale_colour_manual(values = c("white", "grey20"))+
  scale_fill_manual("Function", values = safe_colorblind_palette1)+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  coord_cartesian(ylim = c(1,9), xlim = c(-1.0, 1.2), expand = T)+
  theme(panel.background = element_rect(colour = "grey40", fill = NA), panel.grid = element_blank(),
        strip.background = element_rect(fill = c(col_coast, col_mine)), strip.text = element_text(size = 18),
        aspect.ratio = 1.2, axis.title = element_text(size = 22), axis.text = element_text(size = 18), 
        legend.background = element_blank(), legend.key = element_blank(), legend.title = element_text(size = 22), legend.text = element_text(size = 18))+
  facet_wrap2(~gene_set, ncol = 2, scales = "fixed", strip.position = "top", strip = strip, 
              labeller = labeller(gene_set = gene_set.labs))
salt_bubble_plot3


pdf_fit(file = "salt_bubble_plot3_27_02_24_v2.pdf", pt =20, width = 3.0, w2h = 1)
salt_bubble_plot3
fit.off()




############################### END #################################################


# # old code I used while working on making tables/merging them
# # view(bubble_table)
# # 
# test_DF <- left_join(gene_to_GO_ASP, ASP_goplot, by = "Genes")
# view(test_DF)
# # 
# test_DF <- na.omit(test_DF)
# # 
#  view(test_DF)
# # 
# setnames(sign_fisher_weight_ASP, old = c("GO.ID"), new = c("ID"))
# bubble_plot_table <- left_join(test_DF, sign_fisher_weight_ASP, by = "ID")
# bubble_plot_table2 <- merge(bubble_table, sign_fisher_weight_ASP, by = "ID")
# 
# 
# merge()
# 
# view(bubble_plot_table)

# which(bubble_plot_table$ID == "GO:0010200")
# length(which(bubble_table$ID == "GO:0010200"))
# 
# #calculate the z score:
# 
# #z = # up - # down / count^-1






