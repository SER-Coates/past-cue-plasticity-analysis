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
remove(list=ls())
#set working directory if needed
#setwd()

#0.1 Install and load required packages ####
#install.packages("Pstat")

library(Pstat)
library(tidyverse)
library(data.table)


##0.2. Load functions and useful code ----

###a) splitting strings:

spltVector <- function(vector, pattern, position) {
  spltEach <- function(vector, pattern, position){
    unlist( strsplit(vector, pattern) )[position]
  }
  return( as.vector( mapply(spltEach, vector, pattern, position) ) )
}

#############################################################-
#1.0 IMPORT RESULTS FROM DESEQ2 SCRIPT 1 with dds analysis  ----
###############################################################-

##1.1 Results from the combined data (DEalldds1) normcounts ----

##load counts file for combined experiment data (as out put from script step 1)
norm_counts <- setDT(read.csv("CC_norm_counts3_20_06_24.csv", sep = ","))
cooption_gene_set <- setDT(read.csv("cooption_gene_set_05_07_24.csv", sep = ","))
cue_transfer_table <- setDT(read.csv("cue_transfer_05_07_24.csv"))


dim(norm_counts)

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

Pstat::Pst(control_test_subset, csh = 1, Pw = c("SA", "GR"))

SAGR_control_Pst <-  Pstat::Pst(control_subset, csh = 1, Pw = c("SA", "GR"))

BDPP_control_Pst <- Pstat::Pst(control_subset, csh = 1, Pw = c("BD", "PP"))

str(BDPP_control_Pst)

SAGR_cooption_Pst  <- Pstat::Pst(cooption_Pst_data, csh = 1, boot = 5, Pw = c("SA", "GR"))
BDPP_cooption_Pst <- Pstat::Pst(cooption_Pst_data, csh = 1, boot = 5, Pw = c("BD", "PP"))

min(BDPP_cooption_Pst$Pst_Values) 

#####################################-
##2.3 visualise the results ----
#####################################-


# Create a dual-variable histogram with transparency
minx <- min(SAGR_control_Pst$Pst_Values, SAGR_cooption_Pst$Pst_Values)
maxx <- max(SAGR_control_Pst$Pst_Values, SAGR_cooption_Pst$Pst_Values)

hist(SAGR_control_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(SAGR_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))

hist(BDPP_control_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(BDPP_cooption_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topright", legend=c("control", "coopted_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))


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

CZ_delta_expr <- CZ_counts_wide[,c(1,2,5)]

CZ_Pst_data <- pivot_wider(CZ_delta_expr, names_from = c(gene_name), 
  values_from = zinc_expr_diff)

dim(CZ_Pst_data)

PopCZ <- c(paste(spltVector(CZ_Pst_data$individual, "", 1), spltVector(CZ_Pst_data$individual, "", 2), sep = ""))

CZ_Pst_data$Pop <- PopCZ

dim(CZ_Pst_data)

#put the Pop column in first position:
CZ_Pst_data <- CZ_Pst_data[,c(23095, 1:23094)]

CZ_Pst_data <- CZ_Pst_data[,-2]  

CZ_Pst_data

cue_transfer_Pst_data <- CZ_Pst_data[, c(1, which(colnames(CZ_Pst_data) %in% cue_transfer_genes))]

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

Z_Pst_data_subset <- Z_Pst_data[, c(1, subsetting_vec)] 

#################################-
##3.2 Conduct the analysis ----
#################################-

CZ_Pst_control_subset <- CZ_Pst_data[, c(1, subsetting_vec)] 

SAGR_CZ_Pst <-  Pstat::Pst(CZ_Pst_control_subset, csh = 1, Pw = c("SA", "GR"))

BDPP_CZ_Pst <- Pstat::Pst(CZ_Pst_control_subset, csh = 1, Pw = c("BD", "PP"))

SAGR_cue_transfer_Pst <-  Pstat::Pst(cue_transfer_Pst_data, csh = 1, Pw = c("SA", "GR"))

BDPP_cue_transfer_Pst <- Pstat::Pst(cue_transfer_Pst_data, csh = 1, Pw = c("BD", "PP"))

### 3.2.1 try ZZ comparison instead ----

SAGR_Z_Pst <-  Pstat::Pst(Z_Pst_data_subset, csh = 1, Pw = c("SA", "GR"))

BDPP_Z_Pst <- Pstat::Pst(Z_Pst_data_subset, csh = 1, Pw = c("BD", "PP"))

SAGR_cue_transfer_Z_Pst <-  Pstat::Pst(cue_transfer_Z_Pst_data, csh = 1, Pw = c("SA", "GR"))

BDPP_cue_transfer_Z_Pst <- Pstat::Pst(cue_transfer_Z_Pst_data, csh = 1, Pw = c("BD", "PP"))


#####################################-
##3.3 visualise the results ----
#####################################-

hist(SAGR_CZ_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(SAGR_cue_transfer_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topright", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))

hist(BDPP_CZ_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(BDPP_cue_transfer_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topright", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))

### 3.3.1 plot Z dataset for comparison ----

hist(SAGR_Z_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(SAGR_cue_transfer_Z_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topleft", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))

hist(BDPP_Z_Pst$Pst_Values, xlab="Value",
     ylab="Frequency", 
     col=rgb(0, 0, 1, alpha=0.5))
hist(BDPP_cue_transfer_Z_Pst$Pst_Values, col=rgb(1, 0, 0, alpha=0.5), add=T)
legend("topleft", legend=c("random gene subset", "cue_transfer_genes"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))


