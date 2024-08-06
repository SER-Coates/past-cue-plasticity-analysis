#########################################################################################-
#############                         SARAH COATES                            ###########-
#############           Analysis 3 Pst analysis of expression data              ###########-
#############                     Salt and zinc analysis       v6               ###########-
#############                                                                 ###########-
#########################################################################################-


################################################-
#0. Set up workspace and load requirements ----
################################################-

#clear anything previous:
#remove(list=ls())
#set working directory if needed
#setwd()

#0.1 Install and load required packages ####
#install.packages("Pstat")
#install.packages("ggpubr")

library(Pstat)
library(tidyverse)
library(data.table)
library(ggpubr)
library(pdftools)
library(fplot)


##0.2. Load functions and useful code ----

###a) splitting strings:

spltVector <- function(vector, pattern, position) {
  spltEach <- function(vector, pattern, position){
    unlist( strsplit(vector, pattern) )[position]
  }
  return( as.vector( mapply(spltEach, vector, pattern, position) ) )
}


#colours for plotting
#from colour resources website:
#https://derekogle.com/NCGraphing/resources/colors

seaishgreen <- "#009E73"
dark_hydrangea <- "#CC79A7"
sea_steel_blue <- "#0072B2"
colourblind7 <-  "#44AA99"
"#DDCC77"


#############################################################-
#1.0 IMPORT FILES REQUIRED  ----
###############################################################-

##1.1 Results from the combined data (DEalldds1) normcounts ----

##load counts file for combined experiment data (as out put from script step 1)
norm_counts <- setDT(read.csv("CC_norm_counts3_20_06_24.csv", sep = ","))
dim(norm_counts) # check dimensions

#gene sets to use for Pst calculation
cooption_gene_set <- setDT(read.csv("cooption_gene_set_02_08_24.csv", sep = ","))
cue_transfer_table <- setDT(read.csv("cue_transfer_02_08_24.csv"))
DZP_table <- setDT(read.csv("DZP_same_143_18_07_24.csv"))
evolvedzinc_table <- setDT(read.csv("EDZP_91_18_07_24.csv"))
evolvedcontrols_table <- setDT(read.csv("ECC_same_124_18_07_24.csv"))

## load LFC values for all c vs z comparisons:
resSA_c_z_10_07_24 <- setDT(read.csv("resSA_c_z_10_07_24.csv", sep = ","))
resGR_c_z_10_07_24 <-  setDT(read.csv("resGR_c_z_10_07_24.csv", sep = ","))
resBD_c_z_10_07_24 <- setDT(read.csv("resBD_c_z_10_07_24.csv", sep = ","))
resPP_c_z_10_07_24 <- setDT(read.csv("resPP_c_z_10_07_24.csv", sep = ","))

##1.2 Import Fst results from RAD data analysis (Papadopulos et al., 2021)

Fst_SAGR <- setDT(read.csv("R_input_files/GR_AS_RAD_5kb_fst", sep = "\t"))
Fst_BDPP <- setDT(read.csv("R_input_files/PP_BD_RAD_5kb_fst", sep = "\t"))

#####################################################-
#2.0 Cooption Pst analysis ----
#####################################################-

##2.1 Format table for the control Pst pairwise comparisons analysis ----

norm_counts1 <- norm_counts[,2:5]

#norm_counts1$gene_individual <- paste(norm_counts1$gene_name)

control_counts <- norm_counts1[which(norm_counts1$treatment == "C"),]

head(control_counts)

control_counts

#make the data table wider to have each gene as a column with expression values for each individual:
control_Pst_data <- pivot_wider(control_counts, names_from = c(gene_name), 
                                values_from = norm_count)

print(control_Pst_data[,1], n = 24)

control_Pst_data <- control_Pst_data[,-2] 


#make column for population:

Pop <- c(paste(spltVector(control_Pst_data$individual, "", 1), spltVector(control_Pst_data$individual, "", 2), sep = ""))

control_Pst_data$Pop <- Pop

dim(control_Pst_data)

#put the Pop column in first position:
control_Pst_data1 <- control_Pst_data[,c(23095, 1:23094)]

control_Pst_data1 <- control_Pst_data1[,-2]  

 
#################################-
##2.2 Conduct the analysis ----
#################################-

#run Pst test on the data. Takes way too long. Gonna subset.

#subset data to test out running the null distribution of Pst so it doesn't take too long to run 
subsetting_vec <- sample.int(23094, 2000, replace = F)
subsetting_vec1 <- sample.int(23094, 10, replace = F)

control_subset <- control_Pst_data1[, c(1, subsetting_vec)] 

control_test_subset <- control_Pst_data1[, c(1, subsetting_vec1)]

#get coopted gene subset to do calculations on in table:
coopted_genes <- cooption_gene_set$x

cooption_Pst_data <- control_Pst_data1[, c(1, which(colnames(control_Pst_data1) %in% coopted_genes))] 

cue_transfer_genes <- cue_transfer_table$Names

#cue_transfer_Pst_data <- control_Pst_data1[, c(1, which(colnames(control_Pst_data1) %in% cue_transfer_genes))]

grep

#Pstat::Pst(control_test_subset, csh = 1, Pw = c("SA", "GR"))

#SAGR_control_Pst <-  Pstat::Pst(control_subset, csh = 1, Pw = c("SA", "GR"))

#BDPP_control_Pst <- Pstat::Pst(control_subset, csh = 1, Pw = c("BD", "PP"))

#str(BDPP_control_Pst)

SAGR_cooption_Pst  <- Pstat::Pst(cooption_Pst_data, csh = 1, boot = 5, Pw = c("SA", "GR"))
BDPP_cooption_Pst <- Pstat::Pst(cooption_Pst_data, csh = 1, boot = 5, Pw = c("BD", "PP"))

min(BDPP_cooption_Pst$Pst_Values) 

#####################################-
##2.3 visualise the results ----
#####################################-

# Create a dual-variable histogram with transparency
# minx <- min(SAGR_control_Pst$Pst_Values, SAGR_cooption_Pst$Pst_Values)
# maxx <- max(SAGR_control_Pst$Pst_Values, SAGR_cooption_Pst$Pst_Values)
# 
# hist(SAGR_control_Pst$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(SAGR_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# hist(BDPP_control_Pst$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(BDPP_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))


#use plots of Fst for comparison:
range_coopt_SAGR <- range(SAGR_cooption_Pst$Pst_Values)
range_coopt_BDPP <- range(BDPP_cooption_Pst$Pst_Values)

# hist(Fst_SAGR$MEAN_FST, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(SAGR_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# hist(Fst_SAGR$MEAN_FST, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_coopt_SAGR[1], col="red")
# abline (v =  range_coopt_SAGR[2], col="red")
# 
# hist(Fst_BDPP$MEAN_FST, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_coopt_BDPP[1], col="red")
# abline (v =  range_coopt_BDPP[2], col="red")
# 
# dfstSAGR <- density(Fst_SAGR$MEAN_FST)
# dpstSAGR <- density(SAGR_cooption_Pst$Pst_Values)
# plot(dfstSAGR, xlim = range(dfstSAGR$x, dpstSAGR$x), ylim = range(dfstSAGR$y, dpstSAGR$y), 
#      type = "l", col = "blue")
# lines(dpstSAGR, col = "red")
# 
# dfstBDPP <- density(Fst_BDPP$MEAN_FST)
# dpstBDPP <- density(BDPP_cooption_Pst$Pst_Values)
# plot(dfstBDPP, xlim = range(dfstBDPP$x, dpstBDPP$x),ylim = range(dfstBDPP$y, dpstBDPP$y), 
#      type = "l", col = "blue")
# lines(dpstBDPP, col = "red")

#hist(BDPP_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
#legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))

##########################################-
## 2.4 statistical summary of results
#########################################-

#compare Fst quantile above 95% to the values of the Pst - how many fall into this category?

SA_GR_top_95 <- quantile(x = Fst_SAGR$MEAN_FST, probs = c(0.95, 0.98, 0.99, 0.995))
BD_PP_top_95 <- quantile(x = Fst_BDPP$MEAN_FST, probs = c(0.95, 0.98, 0.99, 0.995))

SA_GR_top_95[1]  

#number above top 95% of all SAGR Fst values - all 33 :)
dim(SAGR_cooption_Pst[which(SAGR_cooption_Pst$Pst_Values > SA_GR_top_95[1]),])

#number above top 98% of all Fst values - all 33 :)
dim(SAGR_cooption_Pst[which(SAGR_cooption_Pst$Pst_Values > SA_GR_top_95[2]),])

#number above top 99.5% 21/33 
dim(SAGR_cooption_Pst[which(SAGR_cooption_Pst$Pst_Values > SA_GR_top_95[4]),])

#number above top 95% of all BD PP Fst values - all 33 also above 95% for BDPP  :)
dim(BDPP_cooption_Pst[which(BDPP_cooption_Pst$Pst_Values > BD_PP_top_95[1]),])

#number above top 98% of all BD PP Fst values - all 33 also above 98% for BDPP  :)
dim(BDPP_cooption_Pst[which(BDPP_cooption_Pst$Pst_Values > BD_PP_top_95[2]),])

#Conclusion - coopted gene Pst values all fall beyond the 98% quantile for Fst
#In both SA vs GR and BD vs PP

#####################################################-
#3.0 Cue transfer Pst analysis ----
#####################################################-

############################################################################-
##3.1 Format table for the control Pst pairwise comparisons analysis ----
############################################################################-

head(norm_counts1)
tail(norm_counts1)

CZ_counts <- norm_counts1[which(norm_counts1$treatment %in% c("C","Z")),]
CZ_counts_zincexp <- CZ_counts[which(CZ_counts$individual %in% c("SA6", "SA7", "SA8", 
                                                                                "BD1","BD11", "BD12", 
                                                                                "GR10", "GR2",  "GR8" , 
                                                                                "PP1" , "PP4",  "PP8")),]

CZ_counts_wide <- pivot_wider(CZ_counts_zincexp, names_from = c(treatment), 
                                values_from = norm_count)

CZ_counts_wide$zinc_expr_diff <- CZ_counts_wide$Z - CZ_counts_wide$C

CZ_counts_wide$zinc_relative_change <- (CZ_counts_wide$Z - CZ_counts_wide$C) / CZ_counts_wide$C

CZ_counts_wide$zinc_zdivc <- CZ_counts_wide$Z/CZ_counts_wide$C

CZ_counts_wide$zinc_log2_zdivc <- log2(CZ_counts_wide$Z/CZ_counts_wide$C)


#view(CZ_counts_wide)

CZ_delta_expr <- CZ_counts_wide[,c(1,2,5)]
CZ_zdivc_expr <- CZ_counts_wide[,c(1,2,7)]

CZ_Pst_data <- pivot_wider(CZ_delta_expr, names_from = c(gene_name), 
  values_from = zinc_expr_diff)

CZ_Pst_data1 <- pivot_wider(CZ_zdivc_expr, names_from = c(gene_name), 
                           values_from = zinc_zdivc)

dim(CZ_Pst_data1)

PopCZ <- c(paste(spltVector(CZ_Pst_data$individual, "", 1), spltVector(CZ_Pst_data$individual, "", 2), sep = ""))

CZ_Pst_data$Pop <- PopCZ
CZ_Pst_data1$Pop <- PopCZ

dim(CZ_Pst_data)

#put the Pop column in first position:
CZ_Pst_data <- CZ_Pst_data[,c(23095, 1:23094)]
CZ_Pst_data1 <- CZ_Pst_data1[,c(23095, 1:23094)]

CZ_Pst_data <- CZ_Pst_data[,-2]  
CZ_Pst_data1 <- CZ_Pst_data1[,-2]  

CZ_Pst_data

cue_transfer_Pst_data <- CZ_Pst_data[, c(1, which(colnames(CZ_Pst_data) %in% cue_transfer_genes))]
cue_transfer_Pst_data1 <- CZ_Pst_data1[, c(1, which(colnames(CZ_Pst_data1) %in% cue_transfer_genes))]

###3.1.1 make ZZ comparison ----
 
Z_counts <- norm_counts1[which(norm_counts1$treatment == "Z"),]

dim(Z_counts)

Z_Pst_data <- pivot_wider(Z_counts, names_from = c(gene_name), 
                              values_from = norm_count)
dim(Z_Pst_data)

Z_Pst_data <- Z_Pst_data[,-2]

PopZ <- c(paste(spltVector(Z_Pst_data$individual, "", 1), spltVector(Z_Pst_data$individual, "", 2), sep = ""))

Z_Pst_data$Pop <- PopZ

Z_Pst_data[,23095]

#put the Pop column in first position:
Z_Pst_data <- Z_Pst_data[,c(23095, 1:23094)]

Z_Pst_data <- Z_Pst_data[,-2]  

Z_Pst_data

cue_transfer_Z_Pst_data <- Z_Pst_data[, c(1, which(colnames(Z_Pst_data) %in% cue_transfer_genes))]

#Z_Pst_data_subset <- Z_Pst_data[, c(1, subsetting_vec)] 

####################################################################-
### 3.1.2 Try log fold change for Zinc vs Control comparisons ----
####################################################################-

#need LFCs for the comparisons between CZ in pp, bd, sa and gr
#probably use shrunken ones if I can.
#only one per population as is a summary statistic, need to calculate a manual version
#from norm counts table.
#Tried just doing z/c as a fraction to plot it.

# colnames(cue_transfer_table)
# 
# colstokeep1 <- c("X", "log2FoldChange") 
# 
# resSA_c_z_10_07_24 <- resSA_c_z_10_07_24[, ..colstokeep1]
# resGR_c_z_10_07_24 <- resGR_c_z_10_07_24[, ..colstokeep1]
# resBD_c_z_10_07_24 <- resBD_c_z_10_07_24[, ..colstokeep1]
# resPP_c_z_10_07_24 <- resPP_c_z_10_07_24[, ..colstokeep1]
# 
# setnames(resSA_c_z_10_07_24, old = c("X", "log2FoldChange"), new = c("gene_name","LFC_CZ_SA"))
# setnames(resGR_c_z_10_07_24, old = c("X", "log2FoldChange"), new = c("gene_name","LFC_CZ_GR"))
# setnames(resBD_c_z_10_07_24, old = c("X", "log2FoldChange"), new = c("gene_name","LFC_CZ_BD"))
# setnames(resPP_c_z_10_07_24, old = c("X", "log2FoldChange"), new = c("gene_name","LFC_CZ_PP"))
# 
# dim(resSA_c_z_10_07_24)
# 
# #make combined LFC sets:
# LFC_SAGR <- merge(resSA_c_z_10_07_24, resGR_c_z_10_07_24, by = "gene_name")
# LFC_BDPP <- merge(resBD_c_z_10_07_24, resPP_c_z_10_07_24, by = "gene_name")
# 
# #try calculating the log fold change rudimentary one for each gene:
# CZ_counts_zincexp

########################################################################-
##3.2 Conduct the analysis ----
########################################################################-

# CZ_Pst_control_subset <- CZ_Pst_data[, c(1, subsetting_vec)] 
# CZ_Pst_control_subset1 <- CZ_Pst_data1[, c(1, subsetting_vec)] 
# 
# SAGR_CZ_Pst <-  Pstat::Pst(CZ_Pst_control_subset, csh = 1, Pw = c("SA", "GR"))
# BDPP_CZ_Pst <- Pstat::Pst(CZ_Pst_control_subset, csh = 1, Pw = c("BD", "PP"))
# 
# SAGR_CZ_Pst1 <-  Pstat::Pst(CZ_Pst_control_subset1, csh = 1, Pw = c("SA", "GR"))
# BDPP_CZ_Pst1 <- Pstat::Pst(CZ_Pst_control_subset1, csh = 1, Pw = c("BD", "PP"))

SAGR_cue_transfer_Pst <-  Pstat::Pst(cue_transfer_Pst_data, csh = 1, Pw = c("SA", "GR"))
BDPP_cue_transfer_Pst <- Pstat::Pst(cue_transfer_Pst_data, csh = 1, Pw = c("BD", "PP"))

SAGR_cue_transfer_Pst1 <-  Pstat::Pst(cue_transfer_Pst_data1, csh = 1, Pw = c("SA", "GR"))
BDPP_cue_transfer_Pst1 <- Pstat::Pst(cue_transfer_Pst_data1, csh = 1, Pw = c("BD", "PP"))

### 3.2.1 try ZZ comparison instead ----

#SAGR_Z_Pst <-  Pstat::Pst(Z_Pst_data_subset, csh = 1, Pw = c("SA", "GR"))

#BDPP_Z_Pst <- Pstat::Pst(Z_Pst_data_subset, csh = 1, Pw = c("BD", "PP"))

SAGR_cue_transfer_Z_Pst <-  Pstat::Pst(cue_transfer_Z_Pst_data, csh = 1, Pw = c("SA", "GR"))

BDPP_cue_transfer_Z_Pst <- Pstat::Pst(cue_transfer_Z_Pst_data, csh = 1, Pw = c("BD", "PP"))

#CZ might be an issue with relative changes being small, so difference in plasticity small.
#if log fold changes it would be scaled to the gene size, so a diff in plasticity of 1 maybe big in the end.

#use log fold change values for C vs Z 

#####################################-
##3.3 visualise the results ----
#####################################-

# #get ranges for plots:
# range_cuet_SAGR <- range(SAGR_cue_transfer_Pst1$Pst_Values)
# range_cuet_BDPP <-range(BDPP_cue_transfer_Pst1$Pst_Values)
# 
# hist(SAGR_CZ_Pst1$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(SAGR_cue_transfer_Pst1$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topright", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# hist(BDPP_CZ_Pst1$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(BDPP_cue_transfer_Pst1$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topright", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# 
# hist(SAGR_CZ_Pst1$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_cuet_SAGR[1], col="red")
# abline (v =  range_cuet_SAGR[2], col="red")
# #abline (v =  SA_GR_top_95[2], col="blue")
# 
# hist(BDPP_CZ_Pst1$Pst_Values, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_cuet_BDPP[1], col="red")
# abline (v =  range_cuet_BDPP[2], col="red")
# #abline (v =  BD_PP_top_95[2], col="blue")
# 
# hist(Fst_SAGR$MEAN_FST, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_cuet_SAGR[1], col="red")
# abline (v =  range_cuet_SAGR[2], col="red")
# #abline (v =  SA_GR_top_95[2], col="blue")
# 
# hist(Fst_BDPP$MEAN_FST, xlab="Value",
#      ylab="Frequency", 
#      col=rgb(0, 0, 1, alpha=0.5))
# abline (v =  range_cuet_BDPP[1], col="red")
# abline (v =  range_cuet_BDPP[2], col="red")
# #abline (v =  BD_PP_top_95[2], col="blue")
# 
# 
# #z/c self-calculated fold change (I think I can call it that) plotted in density plot
# czdfstSAGR <- density(Fst_SAGR$MEAN_FST)
# czdpstSAGR <- density(SAGR_cue_transfer_Pst1$Pst_Values)
# plot(czdfstSAGR, xlim = range(czdfstSAGR$x, czdpstSAGR$x), ylim = range(czdfstSAGR$y, czdpstSAGR$y), 
#      type = "l", col = "blue")
# lines(czdpstSAGR, col = "red")
# legend("center", legend=c("SAGR Fst", "SAGR_cue_transfer_genes"), fill=c("blue", "red"))
# 
# czdfstBDPP <- density(Fst_BDPP$MEAN_FST)
# czdpstBDPP <- density(BDPP_cue_transfer_Pst1$Pst_Values)
# plot(czdfstBDPP, xlim = range(czdfstBDPP$x, czdpstBDPP$x),ylim = range(czdfstBDPP$y, czdpstBDPP$y), 
#      type = "l", col = "blue")
# lines(czdpstBDPP, col = "red")
# 
# 
# ### 3.3.1 plot Z dataset for comparison ----
# 
# #Pst vs Pst fopr Z
# hist(SAGR_Z_Pst$Pst_Values, xlab="Value",
#      ylab="Frequency",
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(SAGR_cue_transfer_Z_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topleft", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# hist(BDPP_Z_Pst$Pst_Values, xlab="Value",
#      ylab="Frequency",
#      col=rgb(0, 0, 1, alpha=0.5))
# hist(BDPP_cue_transfer_Z_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
# legend("topleft", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))
# 
# #FST vs Pst for Z
# zdfstSAGR <- density(Fst_SAGR$MEAN_FST)
# zdpstSAGR <- density(SAGR_cue_transfer_Z_Pst$Pst_Values)
# plot(zdfstSAGR, xlim = range(zdfstSAGR$x, zdpstSAGR$x), ylim = range(zdfstSAGR$y, zdpstSAGR$y), 
#      type = "l", col = "blue")
# lines(zdpstSAGR, col = "red")
# 
# zdfstBDPP <- density(Fst_BDPP$MEAN_FST)
# zdpstBDPP <- density(BDPP_cue_transfer_Z_Pst$Pst_Values)
# plot(dfstBDPP, xlim = range(zdfstBDPP$x, zdpstBDPP$x),ylim = range(zdfstBDPP$y, zdpstBDPP$y), 
#      type = "l", col = "blue")
# lines(zdpstBDPP, col = "red")


############################################################-
#4.0 paper figures for both cue transfer and cooption Pst-Fst results ----
############################################################-

##4.1 create combined Fst-Pst tables:

Fst_SAGR
Fst_BDPP

#edit Fst tables:

# geography <- rep(c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
#                    "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
#                    "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B",
#                    "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"), 23093)

Fst_SAGR$data_type <- rep(c("Fst"), 7691)
Fst_BDPP$data_type <- rep(c("Fst"), 7769)

SAGR_cooption_Pst$data_type <- rep(c("Pst"), 30)
BDPP_cooption_Pst$data_type <- rep(c("Pst"), 30)

SAGR_cue_transfer_Pst1$data_type <- rep(c("Pst"), 28)
BDPP_cue_transfer_Pst1$data_type <- rep(c("Pst"), 28)

#filter and stick data together:

#resSA_c_z_10_07_24 <- resSA_c_z_10_07_24[, ..colstokeep1]

col_filter_Fst <- c("MEAN_FST", "data_type")
col_filter_Pst <- c("Pst_Values", "data_type")

Fst_SAGR_2 <-  Fst_SAGR[, ..col_filter_Fst]
Fst_BDPP_2 <-  Fst_BDPP[, ..col_filter_Fst]

setnames(Fst_SAGR_2, old = c("MEAN_FST"), new = c("Value"))
setnames(Fst_BDPP_2, old = c("MEAN_FST"), new = c("Value"))

SAGR_cooption_Pst_2 <- setDT(SAGR_cooption_Pst)[, ..col_filter_Pst]
BDPP_cooption_Pst_2 <- setDT(BDPP_cooption_Pst)[, ..col_filter_Pst]

SAGR_cue_transfer_Pst1_2 <- setDT(SAGR_cue_transfer_Pst1)[, ..col_filter_Pst]
BDPP_cue_transfer_Pst1_2 <- setDT(BDPP_cue_transfer_Pst1)[, ..col_filter_Pst]

setnames(SAGR_cooption_Pst_2, old = c("Pst_Values"), new = c("Value"))
setnames(BDPP_cooption_Pst_2, old = c("Pst_Values"), new = c("Value"))
setnames(SAGR_cue_transfer_Pst1_2, old = c("Pst_Values"), new = c("Value"))
setnames(BDPP_cue_transfer_Pst1_2, old = c("Pst_Values"), new = c("Value"))

SAGR_cooption_Pst_2
Fst_SAGR_2

SAGR_Fst_Pst_cooption <- rbind(Fst_SAGR_2, SAGR_cooption_Pst_2)
SAGR_Fst_Pst_cue_transfer <- rbind(Fst_SAGR_2, SAGR_cue_transfer_Pst1_2)

BDPP_Fst_Pst_cooption <- rbind(Fst_BDPP_2, BDPP_cooption_Pst_2)
BDPP_Fst_Pst_cue_transfer <- rbind(Fst_BDPP_2, BDPP_cue_transfer_Pst1_2)

dim(SAGR_Fst_Pst_cue_transfer) #7719
dim(BDPP_Fst_Pst_cue_transfer) #7797
dim(SAGR_Fst_Pst_cooption) #7724
dim(BDPP_Fst_Pst_cooption) #7802

#create data-frames for cue transfer and cooption combined:

#add factor I need:
SAGR_Fst_Pst_cue_transfer$location <-  rep(c("Wales"), 7719)
BDPP_Fst_Pst_cue_transfer$location <- rep(c("England"), 7797)

SAGR_Fst_Pst_cooption$location <- rep(c("Wales"), 7721)
BDPP_Fst_Pst_cooption$location <- rep(c("England"), 7799)

#stick the tables together to make 2 tables:
cue_transfer_fst_pst <- rbind(SAGR_Fst_Pst_cue_transfer, BDPP_Fst_Pst_cue_transfer)

cooption_fst_pst <- rbind(SAGR_Fst_Pst_cooption, BDPP_Fst_Pst_cooption)

  
###########################################################-
##4.2 create quantiles maybe to plot onto the graphs too:
###########################################################-

SA_GR_top_95 <- quantile(x = Fst_SAGR$MEAN_FST, probs = c(0.90, 0.95, 0.975, 0.98, 0.99, 0.995))
BD_PP_top_95 <- quantile(x = Fst_BDPP$MEAN_FST, probs = c(0.90, 0.95, 0.975, 0.98, 0.99, 0.995))

#number of genes in cooption that are above 97.5% quantile (top 5% outliers 1 tailed) 
dim(SAGR_cooption_Pst[which(SAGR_cooption_Pst$Pst_Values > SA_GR_top_95[3]),]) #all 33
dim(BDPP_cooption_Pst[which(BDPP_cooption_Pst$Pst_Values > BD_PP_top_95[3]),]) #all 33

#number of genes in cue-transfer that are above 97.5% quantile (top 5% outliers 1 tailed) 
dim(SAGR_cue_transfer_Pst1[which(SAGR_cue_transfer_Pst1$Pst_Values > SA_GR_top_95[3]),]) 
#all 28
dim(BDPP_cue_transfer_Pst1[which(BDPP_cue_transfer_Pst1$Pst_Values > BD_PP_top_95[3]),])
#25/28
#check those in top 10% of outliers: (above 95% quantile)
dim(BDPP_cue_transfer_Pst1[which(BDPP_cue_transfer_Pst1$Pst_Values > BD_PP_top_95[2]),])
#27/28, not bad.
#check those in top 20% of outliers: (above 90% quantile)
dim(BDPP_cue_transfer_Pst1[which(BDPP_cue_transfer_Pst1$Pst_Values > BD_PP_top_95[1]),])
#28

#number of genes in cue-transfer that are above 97.5% quantile (top 5% outliers 1 tailed) 
dim(SAGR_cooption_Pst_2[which(SAGR_cooption_Pst_2$Value > SA_GR_top_95[3]),]) 
#all 33
dim(BDPP_cooption_Pst_2[which(BDPP_cooption_Pst_2$Value > BD_PP_top_95[3]),])
#all 33

###################################-
##4.3 make nice plots ----
###################################-


###4.3.1 make individual plots for each comparison of the 4 comparisons: ----
# Density plot for SAGR cooption
# SAGR_cooption_plot <- ggplot(SAGR_Fst_Pst_cooption, aes(x=Value, 
#                                                         color=data_type, 
#                                                         fill = data_type)) +
#   geom_density(size = 0.6, alpha = 0.4)+
#   geom_vline(xintercept = SA_GR_top_95[3], color = "grey30", size = 0.6)+
#   scale_colour_manual(values = c(seaishgreen, "darkred"))+
#   scale_fill_manual(values = c(seaishgreen, "red"))+
#   theme_classic()+
#   theme(panel.background = element_blank(), 
#         strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
#         aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 18))
# SAGR_cooption_plot
# 
# # Density plot for SAGR cue-transfer 
# SAGR_cue_transfer_plot <- ggplot(SAGR_Fst_Pst_cue_transfer, aes(x=Value, 
#                                                         color=data_type, 
#                                                         fill = data_type)) +
#   geom_density(size = 0.6, alpha = 0.4)+
#   geom_vline(xintercept = SA_GR_top_95[3], color = "grey30", size = 0.6)+
#   scale_colour_manual(values = c(seaishgreen, "blue4"))+
#   scale_fill_manual(values = c(seaishgreen, "blue"))+
#   theme_classic()+
#   theme(panel.background = element_blank(), 
#         strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
#         aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 18))
# SAGR_cue_transfer_plot
# 
# 
# # Density plot for BDPP cooption
# BDPP_cooption_plot <- ggplot(BDPP_Fst_Pst_cooption, aes(x=Value, 
#                                                         color=data_type, 
#                                                         fill = data_type)) +
#   geom_density(size = 0.6, alpha = 0.4)+
#   geom_vline(xintercept = BD_PP_top_95[3], color = "grey30", size = 0.6)+
#   scale_colour_manual(values = c(seaishgreen, "darkred"))+
#   scale_fill_manual(values = c(seaishgreen, "red"))+
#   theme_classic()+
#   theme(panel.background = element_blank(), 
#         strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
#         aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 18))
# BDPP_cooption_plot
# 
# # Density plot for BDPP cue-transfer 
# BDPP_cue_transfer_plot <- ggplot(BDPP_Fst_Pst_cue_transfer, aes(x=Value, 
#                                                                 color=data_type, 
#                                                                 fill = data_type)) +
#   geom_density(size = 0.6, alpha = 0.4)+
#   geom_vline(xintercept = BD_PP_top_95[3], color = "grey30", size = 0.6)+
#   scale_colour_manual(values = c(seaishgreen, "blue4"))+
#   scale_fill_manual(values = c(seaishgreen, "blue"))+
#   theme_classic()+
#   theme(panel.background = element_blank(), 
#         strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 22),
#         aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 18))
# BDPP_cue_transfer_plot
# 
# all4_plots <- ggarrange(SAGR_cue_transfer_plot, 
#                         BDPP_cue_transfer_plot, 
#                         SAGR_cooption_plot, 
#                         BDPP_cooption_plot, 
#                         align = "hv")

#double plots with facet wrap:

# Only colour strips in x-direction
# strip <- strip_themed(background_x = elem_list_rect(fill = c(col_coast, col_mine), color = c(F,F)))
# 
# #make labeller for the plots for each panel on facet wrap: 
# #(name it as the variable name used for wrapping .labs and use the names of the key to be the current labels)
# ecotype.labs <- c("Coast Ecotype", "Mine Ecotype")
# names(ecotype.labs) <- c("cst", "min")

###4.3.2 Make 4-way combined plot with 2-panel plots each for cue transfer and genetic adoption  ----

location_cutoffs <- data.frame(location = c("Wales", "England"), Z = c(SA_GR_top_95[3], BD_PP_top_95[3]))

cue_transfer_pstfst_plot <- ggplot(cue_transfer_fst_pst, aes(x=Value, 
                                                     color=data_type, 
                                                     fill = data_type)) +
  geom_density(size = 0.6, alpha = 0.5)+
  scale_colour_manual(values = c(seaishgreen, "blue4"))+
  scale_fill_manual(values = c(seaishgreen, "blue"))+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 20),
        aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 18),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))+
  facet_wrap(~location, ncol=2)+
  geom_vline(data = location_cutoffs, aes(xintercept = Z))   
#cue_transfer_pstfst_plot 


cooption_pstfst_plot <- ggplot(cooption_fst_pst, aes(x=Value, 
                                                        color=data_type, 
                                                        fill = data_type)) +
  geom_density(size = 0.6, alpha = 0.5)+
  scale_colour_manual(values = c(seaishgreen, "darkred"))+
  scale_fill_manual(values = c(seaishgreen, "red"))+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 20),
        aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))+
 facet_wrap(~location, ncol=2)+
 geom_vline(data = location_cutoffs, aes(xintercept = Z))   
#cooption_pstfst_plot 

cooption_cue_transfer_combined <- ggarrange(ncol = 1, cue_transfer_pstfst_plot, cooption_pstfst_plot, align = "hv")
cooption_cue_transfer_combined

pdf(width = 14, file = "coopt_cue_transfer_PstFst_05_08_24.pdf")
cooption_cue_transfer_combined
dev.off()


# facet_wrap2(~ecotype, ncol = 2, scales = "free_x", strip.position = "top",
#           labeller = labeller(ecotype = ecotype.labs), strip = strip)

#############################################################################-
#5.0 Pst analysis for other gene sets ----
#############################################################################-

#adaptive genes more broadly:

#DZP descendent zinc plasticity - mine zinc response shared
#EDZP evolved descendent zinc plasticity, DZP with evolutionary change in zinc

#ECC evolved control constitutive evolution.

###########################-
##5.1 table formatting
###########################-

#use CZ CZ_Pst_data1 as basis and control_Pst_data1

DZP_table <- setDT(read.csv("DZP_same_143_18_07_24.csv"))
evolvedzinc_table <- setDT(read.csv("EDZP_91_18_07_24.csv"))
evolvedcontrols_table <- setDT(read.csv("ECC_same_124_18_07_24.csv"))

#subset the CZ and control Pst tables for the gene sets of interest:
#plasticity that may be adaptive/ hasevolved zinc values
DZP_Pst_data <- CZ_Pst_data1[, c(1, which(colnames(CZ_Pst_data1) %in% DZP_table$Names))]
EDZP_Pst_data <- CZ_Pst_data1[, c(1, which(colnames(CZ_Pst_data1) %in% evolvedzinc_table$Names))]

#constitutive evolved genes:
ECC_Pst_data <- control_Pst_data1[, c(1, which(colnames(control_Pst_data1) %in% evolvedcontrols_table$Names))]

#########################################################-
##5.2 conduct Pst analysis on plasticity and evolved genes ----
#########################################################-

col_filter_Pst <- c("Pst_Values", "data_type")

SAGR_DZP_Pst_out <-  Pstat::Pst(DZP_Pst_data, csh = 1, Pw = c("SA", "GR"))
BDPP_DZP_Pst_out <- Pstat::Pst(DZP_Pst_data, csh = 1, Pw = c("BD", "PP"))

SAGR_EDZP_Pst_out <-  Pstat::Pst(EDZP_Pst_data, csh = 1, Pw = c("SA", "GR"))
BDPP_EDZP_Pst_out <- Pstat::Pst(EDZP_Pst_data, csh = 1, Pw = c("BD", "PP"))

SAGR_ECC_Pst_out <- Pstat::Pst(ECC_Pst_data, csh = 1, Pw = c("SA", "GR"))
BDPP_ECC_Pst_out <- Pstat::Pst(ECC_Pst_data, csh = 1, Pw = c("BD", "PP"))

Fst_SAGR_2
Fst_BDPP_2

#####################################################################-
##5.3 process the output data to combine with FST for graphs ----
#####################################################################-

#create new columns in the data:
SAGR_DZP_Pst_out$data_type <- rep(c("Pst_DZP"), 143)
BDPP_DZP_Pst_out$data_type <- rep(c("Pst_DZP"), 143)
SAGR_EDZP_Pst_out$data_type <- rep(c("Pst_EDZP"), 91)
BDPP_EDZP_Pst_out$data_type <- rep(c("Pst_EDZP"), 91)
SAGR_ECC_Pst_out$data_type <- rep(c("Pst_ECC"), 124)
BDPP_ECC_Pst_out$data_type <- rep(c("Pst_ECC"), 124)

#check data:
SAGR_DZP_Pst_out
BDPP_DZP_Pst_out
SAGR_EDZP_Pst_out
BDPP_EDZP_Pst_out
SAGR_ECC_Pst_out
BDPP_ECC_Pst_out

#filter column with gene names out into new DT:
SAGR_DZP_Pst_out_2 <- setDT(SAGR_DZP_Pst_out)[, ..col_filter_Pst]
BDPP_DZP_Pst_out_2 <- setDT(BDPP_DZP_Pst_out)[, ..col_filter_Pst]
SAGR_EDZP_Pst_out_2 <- setDT(SAGR_EDZP_Pst_out)[, ..col_filter_Pst]
BDPP_EDZP_Pst_out_2 <- setDT(BDPP_EDZP_Pst_out)[, ..col_filter_Pst]
SAGR_ECC_Pst_out_2 <- setDT(SAGR_ECC_Pst_out)[, ..col_filter_Pst]
BDPP_ECC_Pst_out_2 <- setDT(BDPP_ECC_Pst_out)[, ..col_filter_Pst]

setnames(SAGR_DZP_Pst_out_2, old = "Pst_Values", new = "Value")
setnames(BDPP_DZP_Pst_out_2 , old = "Pst_Values", new = "Value")
setnames(SAGR_EDZP_Pst_out_2, old = "Pst_Values", new = "Value")
setnames(BDPP_EDZP_Pst_out_2, old = "Pst_Values", new = "Value")
setnames(SAGR_ECC_Pst_out_2, old = "Pst_Values", new = "Value")
setnames(BDPP_ECC_Pst_out_2, old = "Pst_Values", new = "Value")

SAGR_Fst_Pst_DZP <- rbind(Fst_SAGR_2, SAGR_DZP_Pst_out_2)
BDPP_Fst_Pst_DZP <- rbind(Fst_BDPP_2, BDPP_DZP_Pst_out_2)
SAGR_Fst_Pst_EDZP <- rbind(Fst_SAGR_2, SAGR_EDZP_Pst_out_2)
BDPP_Fst_Pst_EDZP <- rbind(Fst_BDPP_2, BDPP_EDZP_Pst_out_2)
SAGR_Fst_Pst_ECC <- rbind(Fst_SAGR_2, SAGR_ECC_Pst_out_2)
BDPP_Fst_Pst_ECC <- rbind(Fst_BDPP_2, BDPP_ECC_Pst_out_2)


dim(SAGR_Fst_Pst_DZP)[1] #7834
dim(BDPP_Fst_Pst_DZP)[1] #7912
dim(SAGR_Fst_Pst_EDZP)[1] #7782  
dim(BDPP_Fst_Pst_EDZP)[1]  #7860
dim(SAGR_Fst_Pst_ECC)[1]  #7815
dim(BDPP_Fst_Pst_ECC)[1] #7893


#add factor for location to each PstFST table:
SAGR_Fst_Pst_DZP$location <-  rep(c("Wales"), dim(SAGR_Fst_Pst_DZP)[1]) 
BDPP_Fst_Pst_DZP$location <-  rep(c("England"), dim(BDPP_Fst_Pst_DZP)[1])
SAGR_Fst_Pst_EDZP$location <-  rep(c("Wales"), dim(SAGR_Fst_Pst_EDZP)[1]) 
BDPP_Fst_Pst_EDZP$location <-  rep(c("England"), dim(BDPP_Fst_Pst_EDZP)[1]) 
SAGR_Fst_Pst_ECC$location <-  rep(c("Wales"), dim(SAGR_Fst_Pst_ECC)[1] ) 
BDPP_Fst_Pst_ECC$location <-  rep(c("England"), dim(BDPP_Fst_Pst_ECC)[1]) 


DZP_Fst_Pst <- rbind(SAGR_Fst_Pst_DZP, BDPP_Fst_Pst_DZP)
EDZP_Fst_Pst <- rbind(SAGR_Fst_Pst_EDZP, BDPP_Fst_Pst_EDZP)
ECC_Fst_Pst <- rbind(SAGR_Fst_Pst_ECC, BDPP_Fst_Pst_ECC)


##count number of genes over the 5% threshhold for ECC and DZP ----

#DZP
dim(SAGR_DZP_Pst_out_2[which(SAGR_DZP_Pst_out_2$Value >= SA_GR_top_95[3]),])[1] #116
dim(SAGR_DZP_Pst_out_2[which(SAGR_DZP_Pst_out_2$Value < SA_GR_top_95[3]),])[1] #27 

dim(BDPP_DZP_Pst_out_2[which(BDPP_DZP_Pst_out_2$Value >= BD_PP_top_95[3]),])[1] #106
dim(BDPP_DZP_Pst_out_2[which(BDPP_DZP_Pst_out_2$Value < BD_PP_top_95[3]),])[1] #37

#ECC
dim(SAGR_ECC_Pst_out_2[which(SAGR_ECC_Pst_out_2$Value >= SA_GR_top_95[3]),])[1] #124
dim(SAGR_ECC_Pst_out_2[which(SAGR_ECC_Pst_out_2$Value < SA_GR_top_95[3]),])[1] #0

dim(BDPP_ECC_Pst_out_2[which(BDPP_ECC_Pst_out_2$Value >= BD_PP_top_95[3]),])[1] #124
dim(BDPP_ECC_Pst_out_2[which(BDPP_ECC_Pst_out_2$Value < BD_PP_top_95[3]),])[1] #0

#EDZP
dim(SAGR_EDZP_Pst_out_2[which(SAGR_EDZP_Pst_out_2$Value >= SA_GR_top_95[3]),])[1] #84
dim(SAGR_EDZP_Pst_out_2[which(SAGR_EDZP_Pst_out_2$Value < SA_GR_top_95[3]),])[1] #7

dim(BDPP_EDZP_Pst_out_2[which(BDPP_EDZP_Pst_out_2$Value >= BD_PP_top_95[3]),])[1] #80
dim(BDPP_EDZP_Pst_out_2[which(BDPP_EDZP_Pst_out_2$Value < BD_PP_top_95[3]),])[1] #11


###################################################-
##5.4 plot the results for the 3 gene sets ----
###################################################-


DZP_pstfst_plot <- ggplot(DZP_Fst_Pst, aes(x=Value, 
                                                     color=data_type, 
                                                     fill = data_type)) +
  geom_density(size = 0.6, alpha = 0.5)+
  scale_colour_manual(values = c(seaishgreen, "blue4"))+
  scale_fill_manual(values = c(seaishgreen, "blue"))+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 20),
        aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))+
  facet_wrap(~location, ncol=2)+
  geom_vline(data = location_cutoffs, aes(xintercept = Z))   
DZP_pstfst_plot 


EDZP_pstfst_plot <- ggplot(EDZP_Fst_Pst, aes(x=Value, 
                                           color=data_type, 
                                           fill = data_type)) +
  geom_density(size = 0.6, alpha = 0.5)+
  scale_colour_manual(values = c(seaishgreen, "blue4"))+
  scale_fill_manual(values = c(seaishgreen, "blue"))+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 20),
        aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))+
  facet_wrap(~location, ncol=2)+
  geom_vline(data = location_cutoffs, aes(xintercept = Z))   
EDZP_pstfst_plot 

ECC_pstfst_plot <- ggplot(ECC_Fst_Pst, aes(x=Value, 
                                           color=data_type, 
                                           fill = data_type)) +
  geom_density(size = 0.6, alpha = 0.5)+
  scale_colour_manual(values = c(seaishgreen, "darkred"))+
  scale_fill_manual(values = c(seaishgreen, "red"))+
  theme_classic()+
  theme(panel.background = element_blank(), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = 20),
        aspect.ratio = 1, axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))+
  facet_wrap(~location, ncol=2)+
  geom_vline(data = location_cutoffs, aes(xintercept = Z))   
ECC_pstfst_plot 

EDZP_DZP_ECC_combined <- ggarrange(ncol = 1, DZP_pstfst_plot, EDZP_pstfst_plot, ECC_pstfst_plot, align = "hv")

DZP_ECC_combined <- ggarrange(ncol = 1, DZP_pstfst_plot,  ECC_pstfst_plot, align = "hv")

pdf(width = 30, height = 20, paper = "a4", file = "DZP_ECC_combined_30_07_24.pdf")
DZP_ECC_combined
dev.off()

DZP_ECC_coopt_cuetransfer <- ggarrange(DZP_pstfst_plot, 
                                       ECC_pstfst_plot,
                                       cue_transfer_pstfst_plot, cooption_pstfst_plot, align = "hv")

# pdf(width = 40, height = 40, paper = "a4", file = "DZP_ECC_coopt_cuetransfer_30_07_24.pdf")
# DZP_ECC_coopt_cuetransfer
# dev.off()

setFplot_page(page = "a4", margins = "normal", units = "tw",pt = 20, w2h = 1.8, reset = FALSE)

pdf_fit(file = "DZP_ECC_coopt_cuetransfer_30_07_24_v2.pdf", pt = 16, width = 2.4, w2h = 1.6)
DZP_ECC_coopt_cuetransfer
dev.off()

setFplot_page(page = "a4", margins = "normal", units = "tw",pt = 12, w2h = 1.2, reset = FALSE)

pdf_fit(file = "DZP_ECC_combined_30_07_24_v2.pdf", pt = 12, width = 1.6, w2h = 1.2)
DZP_ECC_combined
dev.off()
