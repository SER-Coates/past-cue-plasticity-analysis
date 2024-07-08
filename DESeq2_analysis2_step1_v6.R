########################################################################-
#############                  SARAH COATES                     #################-
#############                 DEseq2 analysis                   #################-
############# Salt and Zinc expression DE comparisons and PCAs  #################-
########################################################################-

################################################-
#0. Set up workspace and load requirements ----
################################################-

#clear anything previous:
remove(list=ls())
#set working directory if needed
#setwd()
#0.1 Install and load required packages ####

#package installation:
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("DESeq2")
# install.packages(tidyverse)
# install.packages("data.table")
# install.packages("fplot")
# install.packages("pdftools")

#load packages
library("DESeq2")
library(tidyverse)
library(pdftools)
library(fplot)

##########################-
##0.2. Load functions ----
###########################-

#a) make notin operator from in operator:
`%notin%` <- Negate(`%in%`)

#b)modified plotPCA function for PCA 3 and 4 only
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[3:4]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

#########################################################-
#0.3 set up my colour palette for later use in figures: ----
##########################################################-

#load colours for graphs
col_coast = "#87CEFA" #lightblue for coast
col_mine = "orange"  #orange for mine
col_PPBD = "#046C9A" #dark blue for england
col_GRSA = "#F21A00" #bright_red for wales

###############################################-
#1.set-up the gene count matrices ----
###############################################-

####################################################################-
##1.1 Load gene count matrices and phenotypes and edit names ----
####################################################################-

#read in the count matrices
scountdata <- as.matrix(read.csv("R_input_files/salt_gene_count_matrix.csv", sep = ",", row.names = "gene_id"))
zcountdata <- as.matrix(read.csv("R_input_files/zinc_gene_count_matrix.csv", sep = ",", row.names = "gene_id"))

#read in the phenotype matrices
sphenotypes <- read.csv("R_input_files/salt_phenotypes.csv", sep =",", row.names = "Sample")
zphenotypes <- read.csv("R_input_files/zinc_phenotypes.csv", sep =",", row.names = "Sample")
allphenotypes <- read.csv("R_input_files/all_phenotypes.csv", sep = ",", row.names = "Sample")

snewnames <- c("BD03_C_17","BD03_S_18","BD05_C_9","BD05_S_10","BD07_C_19","BD07_S_20",
              "GR-RNA-10_C_13", "GR-RNA-10_S_14", "GR-RNA-12_C_23","GR-RNA-12_S_24", "GR-RNA-6_C_5","GR-RNA-6_S_6",
              "PP-RNA-1_C_1","PP-RNA-1_S_2", "PP12_C_21", "PP12_S_22", "PP1_C_3", "PP1_S_4", 
              "SA02_C_7", "SA02_S_8", "SA04_C_15", "SA04_S_16", "SA07_C_11", "SA07_S_12")
colnames(scountdata) <- snewnames
colnames(scountdata)

znewnames <- c("BD11Z", "BD12C", "BD1Z", "GR10C", "GR10Z", "GR2C", "GR8C", "GR8Z", "PP1C", "BD1C", "PP1Z", "PP4C",
                   "PP4Z", "PP8C", "SA6Z", "SA7C", "SA7Z", "SA8C", "SA8Z", "BD11C", "BD12Z", "GR2Z", "PP8Z", "SA6C")
length(znewnames)
colnames(zcountdata) <- znewnames
colnames(zcountdata)

#convert characters to factors in the phenotypes table for each S, z and all:
#convert characters to factors in the sphenotypes table
sphenotypes$Treatment = as.factor(sphenotypes$Treatment)
sphenotypes$Ecotype = as.factor(sphenotypes$Ecotype)
sphenotypes$Population = as.factor(sphenotypes$Population)
sphenotypes$Location = as.factor(sphenotypes$Location)
sphenotypes$Pop_treat = as.factor(sphenotypes$Pop_treat) #this one is important!
sphenotypes$Individual = as.factor(sphenotypes$Individual)

#check data structure
str(sphenotypes)

#check names are the same
all(colnames(scountdata) == rownames(sphenotypes))

#convert characters to factors in the zphenotypes table
zphenotypes$Treatment = as.factor(zphenotypes$Treatment)
zphenotypes$Ecotype = as.factor(zphenotypes$Ecotype)
zphenotypes$Population = as.factor(zphenotypes$Population)
zphenotypes$Location = as.factor(zphenotypes$Location)
zphenotypes$Pop_treat = as.factor(zphenotypes$Pop_treat)
zphenotypes$Individual = as.factor(zphenotypes$Individual)

#check data structure
str(zphenotypes)
#check names are the same
all(colnames(zcountdata) == rownames(zphenotypes))

#restructure the combined phenotype file:
allphenotypes$Treatment = as.factor(allphenotypes$Treatment)
allphenotypes$Ecotype = as.factor(allphenotypes$Ecotype)
allphenotypes$Population = as.factor(allphenotypes$Population)
allphenotypes$Location = as.factor(allphenotypes$Location)
allphenotypes$Pop_treat = as.factor(allphenotypes$Pop_treat)
allphenotypes$Individual = as.factor(allphenotypes$Individual)

#check data structure
str(allphenotypes)
rownames(allphenotypes)

############################################-
##1.2 Count matrix filtering and combining both matrices ----
#############################################-

#check initial dimensions are the same.
dim(zcountdata)[1] #41603
dim(scountdata)[1] #41603

#explore raw count matrix structure
#salt data
summary(rowSums(scountdata))
hist(rowSums((scountdata)), breaks = 5000, xlim = c(0,20000))

#zinc data
hist(rowSums((zcountdata)), breaks = 5000, xlim = c(0,20000))
summary(rowSums(zcountdata))

###1.2.1 create filters and filter s and z data by them, then combine ----

#check for appropriate filter
sum(rowSums(zcountdata) >= 10) #31779
sum(rowSums(scountdata) >= 10) #30714

#create filters:
keeps <- rowSums(scountdata) >= 10
keepz <- rowSums(zcountdata) >= 10

#filter the data:
scountdata_filt <- scountdata[keeps,]
dim(scountdata_filt)[1]
zcountdata_filt <- zcountdata[keepz,]
dim(zcountdata_filt)[1]
head(scountdata_filt)

### 1.2.2 combine the newly filtered s and the z count matrices into 1 big one  ----

#create combined dataset 
intersect(rownames(scountdata_filt), rownames(zcountdata_filt))
#merge both sets into 1
allcountdata <- merge(scountdata_filt, zcountdata_filt, by = "row.names")
str(scountdata_filt)

#checks on the structure of the data:
dim(allcountdata)[1]
head(allcountdata)
colnames(allcountdata)
str(allcountdata)

#convert the new column that merge created to row.names again rather than a column:
rownames(allcountdata) <- allcountdata$Row.names #create new row names
allcountdata1 <- allcountdata[,-1]
as.matrix(allcountdata1)
str(allcountdata1)

#check data distribution for the matrix:
hist(rowSums((allcountdata1)), breaks = 5000, xlim = c(0,20000))
summary(rowSums(allcountdata1))
dim(allcountdata1)


###################################################################-
# 2. DEseq2 set-up of dds objects ----
####################################################################-

############################################################################-
##2.1 Set up DEseq matrix with the combined population-treatment factor ----
############################################################################-

#set up dds object for all combined experimental data
alldds1 <- DESeqDataSetFromMatrix(countData = allcountdata1,
                                colData = allphenotypes,
                                design = ~ Pop_treat)
alldds1

#check number of rows in dataframe
nrow(alldds1) #30178

#set up dds object for salt data
sdds1 <- DESeqDataSetFromMatrix(countData = scountdata_filt,
                               colData = sphenotypes,
                               design = ~ Pop_treat)
sdds1

#check number of rows in dataframe
nrow(sdds1) #30714

#set up dds object for zinc data
zdds1 <- DESeqDataSetFromMatrix(countData = zcountdata_filt,
                                colData = zphenotypes,
                                design = ~ Pop_treat)
zdds1

#check number of rows in dataframe
nrow(zdds1) #31779

################################################################ -
##2.2 generate dds2 objects for salt and zinc: ----
#################################################################-

#salt dds2 object:
sdds2 <- DESeqDataSetFromMatrix(countData = scountdata_filt,
                               colData = sphenotypes,
                               design = ~ Population + Population:Individual + Population:Treatment)
dim(sdds2)

#zinc dds2 object:
zdds2 <- DESeqDataSetFromMatrix(countData = zcountdata_filt,
                                colData = zphenotypes,
                                design = ~ Population + Population:Individual + Population:Treatment)
dim(zdds2) 

sum(rowSums(counts(sdds2)) >= 10) #30714
sum(rowSums(counts(zdds2)) >= 10) #31779

##################################################################-
#3. run DE seq analysis for all comparisons for salt, zinc and the combined set ----
#####################################################################-

##3.1 run dds 1 objects for all, s and z
DEalldds1 <- DESeq(alldds1)
DEsdds1 <- DESeq(sdds1)
DEzdds1 <- DESeq(zdds1)

#run for dds2 objects for salt and zinc
DEsdds2 <- DESeq(sdds2)
DEzdds2 <- DESeq(zdds2)

#check output
table(mcols(DEalldds1)$betaConv, useNA="always")
table(mcols(DEsdds1)$betaConv, useNA="always")
table(mcols(DEzdds1)$betaConv, useNA="always")
table(mcols(DEsdds2)$betaConv, useNA="always")
table(mcols(DEzdds2)$betaConv, useNA="always")

############################################################-
##3.2 check results for groups that haven't converged ----
############################################################-

#look at names of those that didn't converge to have a look at their expression
zinc_non_converged_names <- c(row.names(DEzdds2[which(mcols(DEzdds2)$betaConv == FALSE), ]))
zinc_non_converged_names
#use the names to pull out the relevant loci from the countdata table and have a look at them
zcountdata_non_converged <- zcountdata_filt[which(row.names(zcountdata_filt) %in% zinc_non_converged_names == TRUE), ]
hist(rowSums(zcountdata_non_converged))
summary(rowSums(zcountdata_non_converged))
summary(rowSums(zcountdata_filt))

#check dimensions of all dds2 results for zinc
dim(DEzdds2)[1] #31779 

##########################################################-
#4. Get results for the contrasts of interest  ----
############################################################-

#######################################################################-
##4.1 Results from contrasts using alldds1 ----
#######################################################################-

#control comparisons to find those genes that are consistent in expression between experimental controls:
res_SA_C1_2 <- results(DEalldds1, contrast=c("Pop_treat", "SA_C1", "SA_C2"), alpha = 0.05)
res_SA_C1_2$Names <- row.names(res_SA_C1_2)
res_BD_C1_2 <- results(DEalldds1, contrast=c("Pop_treat", "BD_C1", "BD_C2"), alpha = 0.05) 
res_BD_C1_2$Names <- row.names(res_BD_C1_2)
res_GR_C1_2 <- results(DEalldds1, contrast=c("Pop_treat", "GR_C1", "GR_C2"), alpha = 0.05) 
res_GR_C1_2$Names <- row.names(res_GR_C1_2)
res_PP_C1_2 <- results(DEalldds1, contrast=c("Pop_treat", "PP_C1", "PP_C2"), alpha = 0.05) 
res_PP_C1_2$Names <- row.names(res_PP_C1_2)

lfc_SA_C1_2 = lfcShrink(DEalldds1, contrast=c("Pop_treat", "SA_C1", "BD_C1"), type="ashr")
lfc_SA_C1_2$Names = rownames(lfc_SA_C1_2)
view(lfc_SA_C1_2)

#compare salt treatment in coasts to zinc treatment in mines
res_SZ_grsa <- results(DEalldds1, contrast=c("Pop_treat", "GR_Z", "SA_S"), alpha = 0.05)
res_SZ_grsa$Names <- row.names(res_SZ_grsa)
res_SZ_ppbd <- results(DEalldds1, contrast=c("Pop_treat", "PP_Z", "BD_S"), alpha = 0.05)
res_SZ_ppbd$Names <- row.names(res_SZ_ppbd)

#compare salt vs zinc in individual coasts and individual mines
res_gr_SZ <- results(DEalldds1, contrast=c("Pop_treat", "GR_Z", "GR_S"), alpha = 0.05)
res_gr_SZ$Names <- row.names(res_gr_SZ )
res_pp_SZ <- results(DEalldds1, contrast=c("Pop_treat", "PP_Z", "PP_S"), alpha = 0.05)
res_pp_SZ$Names <- row.names(res_pp_SZ)
res_sa_SZ <- results(DEalldds1, contrast=c("Pop_treat", "SA_Z", "SA_S"), alpha = 0.05)
res_sa_SZ$Names <- row.names(res_sa_SZ)
res_bd_SZ <- results(DEalldds1, contrast=c("Pop_treat", "BD_Z", "BD_S"), alpha = 0.05)
res_bd_SZ$Names <- row.names(res_bd_SZ )

#remove any NAs by taking only non NA data entries
res_SA_C1_2 <- res_SA_C1_2[is.na(res_SA_C1_2$padj) == FALSE,]
res_BD_C1_2 <- res_BD_C1_2[is.na(res_BD_C1_2$padj) == FALSE,]
res_GR_C1_2 <- res_GR_C1_2[is.na(res_GR_C1_2$padj) == FALSE,]
res_PP_C1_2 <- res_PP_C1_2[is.na(res_PP_C1_2$padj) == FALSE,]

#S vs zinc comparisons within and between populations
res_SZ_grsa <- res_SZ_grsa[is.na(res_SZ_grsa$padj) == FALSE,] 
res_SZ_ppbd <- res_SZ_ppbd[is.na(res_SZ_ppbd$padj) == FALSE,]
res_gr_SZ <- res_gr_SZ[is.na(res_gr_SZ$padj) == FALSE,]
res_pp_SZ <- res_pp_SZ[is.na(res_pp_SZ$padj) == FALSE,]
res_sa_SZ <- res_sa_SZ[is.na(res_sa_SZ$padj) == FALSE,]
res_bd_SZ <- res_bd_SZ[is.na(res_bd_SZ$padj) == FALSE,]

#get those across C1 and C2 not sig DE. (p > 0.05)
SA_C1_2_no_DE  <- res_SA_C1_2[res_SA_C1_2$padj > 0.05,]
BD_C1_2_no_DE  <- res_BD_C1_2[res_BD_C1_2$padj > 0.05,]
GR_C1_2_no_DE  <- res_GR_C1_2[res_GR_C1_2$padj > 0.05,]
PP_C1_2_no_DE  <- res_PP_C1_2[res_PP_C1_2$padj > 0.05,]

#filter SZ comparisons for those not sig DE between S and Z (p > 0.05)
res_SZ_grsa_no_DE <- res_SZ_grsa[res_SZ_grsa$padj > 0.05,] 
res_SZ_ppbd_no_DE <- res_SZ_ppbd[res_SZ_ppbd$padj  > 0.05,]
res_gr_SZ_no_DE <- res_gr_SZ[res_gr_SZ$padj  > 0.05,]
res_pp_SZ_no_DE <- res_pp_SZ[res_pp_SZ$padj  > 0.05,]
res_sa_SZ_no_DE <- res_sa_SZ[res_sa_SZ$padj  > 0.05,]
res_bd_SZ_no_DE <- res_bd_SZ[res_bd_SZ$padj  > 0.05,]

#filter SZ comparisons for those sig DE between S and Z (p < 0.05)
res_SZ_grsa_sig <- res_SZ_grsa[res_SZ_grsa$padj < 0.05,] 
res_SZ_ppbd_sig <- res_SZ_ppbd[res_SZ_ppbd$padj  < 0.05,]
res_gr_SZ_sig <- res_gr_SZ[res_gr_SZ$padj  < 0.05,]
res_pp_SZ_sig <- res_pp_SZ[res_pp_SZ$padj  < 0.05,]
res_sa_SZ_sig <- res_sa_SZ[res_sa_SZ$padj  < 0.05,]
res_bd_SZ_sig <- res_bd_SZ[res_bd_SZ$padj  < 0.05,]

#control comparisons for all populations
dim(SA_C1_2_no_DE)[1] 
dim(BD_C1_2_no_DE)[1]  
dim(GR_C1_2_no_DE)[1]  
dim(PP_C1_2_no_DE)[1] 
#not significant DE results:
dim(res_SZ_grsa_no_DE)[1]
dim(res_SZ_ppbd_no_DE)[1]
dim(res_gr_SZ_no_DE)[1]
dim(res_pp_SZ_no_DE)[1]
dim(res_sa_SZ_no_DE)[1]
dim(res_bd_SZ_no_DE)[1]
#sig. DE results:
dim(res_SZ_grsa_sig)[1]
dim(res_SZ_ppbd_sig)[1]
dim(res_gr_SZ_sig)[1]
dim(res_pp_SZ_sig)[1]
dim(res_sa_SZ_sig)[1]
dim(res_bd_SZ_sig)[1]

#write out the individual genes that are not significantly differentially expressed between the 2 controls for each pop:
write.csv(SA_C1_2_no_DE, "SA_C1_2_no_DE_18_05_23.csv")
write.csv(BD_C1_2_no_DE, "BD_C1_2_no_DE_18_05_23.csv")
write.csv(GR_C1_2_no_DE, "GR_C1_2_no_DE_18_05_23.csv")
write.csv(PP_C1_2_no_DE, "PP_C1_2_no_DE_18_05_23.csv")

#write out genes that are DE and not DE from sets analysed in may 23 (zinc/salt comparisons)
#no significant DE:
write.csv(res_SZ_grsa_no_DE, "res_SZ_grsa_no_DE_18_05_23.csv")
write.csv(res_SZ_ppbd_no_DE, "res_SZ_ppbd_no_DE_18_05_23.csv")
write.csv(res_gr_SZ_no_DE, "res_gr_SZ_no_DE_18_05_23.csv")
write.csv(res_pp_SZ_no_DE, "res_pp_SZ_no_DE_18_05_23.csv")
write.csv(res_sa_SZ_no_DE, "res_sa_SZ_no_DE_18_05_23.csv")
write.csv(res_bd_SZ_no_DE, "res_bd_SZ_no_DE_18_05_23.csv")
#significant DE:
write.csv(res_SZ_grsa_sig, "res_SZ_grsa_sig_18_05_23.csv")
write.csv(res_SZ_ppbd_sig, "res_SZ_ppbd_sig_18_05_23.csv")
write.csv(res_gr_SZ_sig, "res_gr_SZ_sig_18_05_23.csv")
write.csv(res_pp_SZ_sig, "res_pp_SZ_sig_18_05_23.csv")
write.csv(res_sa_SZ_sig, "res_sa_SZ_sig_18_05_23.csv")
write.csv(res_bd_SZ_sig, "res_bd_SZ_sig_18_05_23.csv")

################################################################################-
##4.2 Salt experiment DEsdds1 results analysis (EC in control, salt) ----
################################################################################-

#Use the coastal a denominator for geographical pair contrasts
#contrast = c('factorName','numeratorLevel','denominatorLevel'),

#checking the contents of the table of results to see what info is included.
#mcols(resGR_SA_C, use.names = TRUE)

#mines vs coasts in control and salt - EC in control and EC salt
resGR_SA_C <- results(DEsdds1, contrast=c("Pop_treat", "GR_C", "SA_C"), alpha = 0.05)
resGR_SA_C$Names <- row.names(resGR_SA_C)
summary(resGR_SA_C)
resPP_BD_C <- results(DEsdds1, contrast=c("Pop_treat", "PP_C", "BD_C"),  alpha = 0.05)
resPP_BD_C$Names <- row.names(resPP_BD_C)
summary(resPP_BD_C)
resGR_SA_S <- results(DEsdds1, contrast=c("Pop_treat", "GR_S", "SA_S"),  alpha = 0.05)
resGR_SA_S$Names <- row.names(resGR_SA_S)
summary(resGR_SA_S)
resPP_BD_S <- results(DEsdds1, contrast=c("Pop_treat", "PP_S", "BD_S"),  alpha = 0.05)
resPP_BD_S$Names <- row.names(resPP_BD_S)
summary(resPP_BD_S)

#Remove NAs: removes p-values for genes where counts are too low to get reliable estimate, as specified in user guide.
resGR_SA_C <- resGR_SA_C[is.na(resGR_SA_C$padj) == FALSE,]
resPP_BD_C <- resPP_BD_C[is.na(resPP_BD_C$padj) == FALSE,]
resGR_SA_S <- resGR_SA_S[is.na(resGR_SA_S$padj) == FALSE,]
resPP_BD_S <- resPP_BD_S[is.na(resPP_BD_S$padj) == FALSE,]

#filter results for those with significant DE (P < 0.05)
resGR_SA_C_sig <- resGR_SA_C[resGR_SA_C$padj < 0.05,]
resPP_BD_C_sig <- resPP_BD_C[resPP_BD_C$padj < 0.05,]
resGR_SA_S_sig <- resGR_SA_S[resGR_SA_S$padj < 0.05,]
resPP_BD_S_sig <- resPP_BD_S[resPP_BD_S$padj < 0.05,]

resGR_SA_C_noDE <- resGR_SA_C[resGR_SA_C$padj > 0.05,]
resPP_BD_C_noDE <- resPP_BD_C[resPP_BD_C$padj > 0.05,]
resGR_SA_S_noDE <- resGR_SA_S[resGR_SA_S$padj > 0.05,]
resPP_BD_S_noDE <- resPP_BD_S[resPP_BD_S$padj > 0.05,]

## write out results to files:
# mine comparison with coasts in control and salt:
write.csv(resGR_SA_C_sig, "resGR_SA_C_sig_18_05_23.csv")
write.csv(resPP_BD_C_sig, "resPP_BD_C_sig_18_05_23.csv")
write.csv(resGR_SA_S_sig, "resGR_SA_S_sig_18_05_23.csv")
write.csv(resPP_BD_S_sig, "resPP_BD_S_sig_18_05_23.csv")

write.csv(resGR_SA_C_noDE, "resGR_SA_C_noDE_18_05_23.csv")
write.csv(resPP_BD_C_noDE, "resPP_BD_C_noDE_18_05_23.csv")
write.csv(resGR_SA_S_noDE, "resGR_SA_S_noDE_18_05_23.csv")
write.csv(resPP_BD_S_noDE, "resPP_BD_S_noDE_18_05_23.csv")

################################################################################-
##4.3 Salt experiment sdds2 results analysis (salt plasticity) ----
################################################################################-

#results for differential expression for pairwise comparisons of population and treatment
#for each population test the difference in expression for salt vs control*population
resSA_c_s <- results(DEsdds2, contrast = list("PopulationSA.TreatmentS"), alpha = 0.05) #comparison i
resSA_c_s$Names <- row.names(resSA_c_s)
resBD_c_s <- results(DEsdds2, contrast = list("PopulationBD.TreatmentS"), alpha = 0.05) # comparison j
resBD_c_s$Names <- row.names(resBD_c_s)
resGR_c_s <- results(DEsdds2, contrast = list("PopulationGR.TreatmentS"), alpha = 0.05) #comparison k
resGR_c_s$Names <- row.names(resGR_c_s)
resPP_c_s <- results(DEsdds2, contrast = list("PopulationPP.TreatmentS"), alpha = 0.05) # comparison l
resPP_c_s$Names <- row.names(resPP_c_s)

#dimensions check
dim(resSA_c_s)[1]
dim(resBD_c_s)[1]
dim(resGR_c_s)[1]
dim(resPP_c_s)[1]

#summary of results to check:
summary(resSA_c_s)
summary(resBD_c_s)
summary(resGR_c_s)
summary(resPP_c_s)

#get significant and non-significant gene sets.

#remove any NAs by taking only non NA data entries - removing low counts!
resSA_c_s <- resSA_c_s[is.na(resSA_c_s$padj) == FALSE,]
resBD_c_s <- resBD_c_s[is.na(resBD_c_s$padj) == FALSE,]
resGR_c_s <- resGR_c_s[is.na(resGR_c_s$padj) == FALSE,]
resPP_c_s <- resPP_c_s[is.na(resPP_c_s$padj) == FALSE,]

#dimensions check
dim(resSA_c_s)[1]
dim(resBD_c_s)[1]
dim(resGR_c_s)[1]
dim(resPP_c_s)[1]

#Take those results that are less than p adjusted of < 0.05
resSA_c_s_sig <- resSA_c_s[resSA_c_s$padj < 0.05,]
resBD_c_s_sig <- resBD_c_s[resBD_c_s$padj < 0.05,]
resGR_c_s_sig <- resGR_c_s[resGR_c_s$padj < 0.05,]
resPP_c_s_sig <- resPP_c_s[resPP_c_s$padj < 0.05,]

#dimensions check
dim(resSA_c_s_sig)[1]
dim(resBD_c_s)[1]
dim(resGR_c_s)[1]
dim(resPP_c_s)[1]

#filter results for non sign. ones P > 0.05
resSA_c_s_noDE <- resSA_c_s[resSA_c_s$padj > 0.05,]
resBD_c_s_noDE <- resBD_c_s[resBD_c_s$padj > 0.05,]
resGR_c_s_noDE <- resGR_c_s[resGR_c_s$padj > 0.05,]
resPP_c_s_noDE <- resPP_c_s[resPP_c_s$padj > 0.05,]

#dimensions check
dim(resSA_c_s_sig)[1]
dim(resBD_c_s_sig)[1]
dim(resGR_c_s_sig)[1]
dim(resPP_c_s_sig)[1]
dim(resSA_c_s_noDE)[1]
dim(resBD_c_s_noDE)[1]
dim(resGR_c_s_noDE)[1]
dim(resPP_c_s_noDE)[1]

#wrote the DE genes to files:
#comment out when run:
write.csv(resSA_c_s_sig, "resSA_c_s_sig_18_05_23.csv")
write.csv(resBD_c_s_sig, "resBD_c_s_sig_18_05_23.csv")
write.csv(resGR_c_s_sig, "resGR_c_s_sig_18_05_23.csv")
write.csv(resPP_c_s_sig, "resPP_c_s_sig_18_05_23.csv")
write.csv(resSA_c_s_noDE, "resSA_c_s_noDE_18_05_23.csv")
write.csv(resBD_c_s_noDE, "resBD_c_s_noDE_18_05_23.csv")
write.csv(resGR_c_s_noDE, "resGR_c_s_noDE_18_05_23.csv")
write.csv(resPP_c_s_noDE, "resPP_c_s_noDE_18_05_23.csv")

################################################################################-
##4.4 Zinc experiment zdds1 results analysis (zinc population comparisons) ----
################################################################################-

#### C vs C and S vs S contrast denominator is the coast  - see how mine changes relative to this 
#1. GR/SA in control
resZ_GR_SA_C <- results(DEzdds1, contrast=c("Pop_treat", "GR_C", "SA_C"), alpha = 0.05)
resZ_GR_SA_C$Names <- row.names(resZ_GR_SA_C)
#2. PP vs BD in control
resZ_PP_BD_C <- results(DEzdds1, contrast=c("Pop_treat", "PP_C", "BD_C"),  alpha = 0.05)
resZ_PP_BD_C$Names <- row.names(resZ_PP_BD_C)
#3. GR vs SA in zinc
resZ_GR_SA_Z <- results(DEzdds1, contrast=c("Pop_treat", "GR_Z", "SA_Z"),  alpha = 0.05)
resZ_GR_SA_Z$Names <- row.names(resZ_GR_SA_Z)
#4. PP vs BD in zinc
resZ_PP_BD_Z <- results(DEzdds1, contrast=c("Pop_treat", "PP_Z", "BD_Z"),  alpha = 0.05)
resZ_PP_BD_Z$Names <- row.names(resZ_PP_BD_Z)

#check dimensions
dim(resZ_GR_SA_C)[1]
dim(resZ_PP_BD_C)[1]
dim(resZ_GR_SA_Z)[1]
dim(resZ_PP_BD_Z)[1]

#summaries of the results:
summary(resZ_GR_SA_C)
summary(resZ_PP_BD_C)
summary(resZ_GR_SA_Z)
summary(resZ_PP_BD_Z)

#Remove NAs: 
#removes p-values for genes where counts are too low to get reliable estimate, as specified in user guide.
resZ_GR_SA_C <- resZ_GR_SA_C[is.na(resZ_GR_SA_C$padj) == FALSE,]
resZ_PP_BD_C <- resZ_PP_BD_C[is.na(resZ_PP_BD_C$padj) == FALSE,]
resZ_GR_SA_Z <- resZ_GR_SA_Z[is.na(resZ_GR_SA_Z$padj) == FALSE,]
resZ_PP_BD_Z <- resZ_PP_BD_Z[is.na(resZ_PP_BD_Z$padj) == FALSE,]

#check dimensions
dim(resZ_GR_SA_C)[1]
dim(resZ_PP_BD_C)[1]
dim(resZ_GR_SA_Z)[1]
dim(resZ_PP_BD_Z)[1]

#Significantly different genes: Adjusted p-value < 0.05. 
resZ_GR_SA_C_sig <- resZ_GR_SA_C[resZ_GR_SA_C$padj < 0.05,]
resZ_PP_BD_C_sig <- resZ_PP_BD_C[resZ_PP_BD_C$padj < 0.05,]
resZ_GR_SA_Z_sig <- resZ_GR_SA_Z[resZ_GR_SA_Z$padj < 0.05,]
resZ_PP_BD_Z_sig <- resZ_PP_BD_Z[resZ_PP_BD_Z$padj < 0.05,]

#check dimensions
dim(resZ_GR_SA_C_sig)[1]
dim(resZ_PP_BD_C_sig)[1]
dim(resZ_GR_SA_Z_sig)[1]
dim(resZ_PP_BD_Z_sig)[1]

#not DE sig genes: Adjusted p-value > 0.05. 
resZ_GR_SA_C_noDE <- resZ_GR_SA_C[resZ_GR_SA_C$padj > 0.05,]
resZ_PP_BD_C_noDE <- resZ_PP_BD_C[resZ_PP_BD_C$padj > 0.05,]
resZ_GR_SA_Z_noDE <- resZ_GR_SA_Z[resZ_GR_SA_Z$padj > 0.05,]
resZ_PP_BD_Z_noDE <- resZ_PP_BD_Z[resZ_PP_BD_Z$padj > 0.05,]

#check dimensions
dim(resZ_GR_SA_C_noDE)[1]
dim(resZ_PP_BD_C_noDE)[1]
dim(resZ_GR_SA_Z_noDE)[1]
dim(resZ_PP_BD_Z_noDE)[1]

# #write to files:
write.csv(resZ_GR_SA_C_sig, "resZ_GR_SA_C_sig_18_05_23.csv")
write.csv(resZ_PP_BD_C_sig, "resZ_PP_BD_C_sig_18_05_23.csv")
write.csv(resZ_GR_SA_Z_sig, "resZ_GR_SA_Z_sig_18_05_23.csv")
write.csv(resZ_PP_BD_Z_sig, "resZ_PP_BD_Z_sig_18_05_23.csv")

write.csv(resZ_GR_SA_C_noDE, "resZ_GR_SA_C_noDE_18_05_23.csv")
write.csv(resZ_PP_BD_C_noDE, "resZ_PP_BD_C_noDE_18_05_23.csv")
write.csv(resZ_GR_SA_Z_noDE, "resZ_GR_SA_Z_noDE_18_05_23.csv")
write.csv(resZ_PP_BD_Z_noDE, "resZ_PP_BD_Z_noDE_18_05_23.csv")

################################################################################-
##4.5 Zinc experiment zdds2 results analysis (Zinc plasticity) ----
################################################################################-

resSA_c_z <- results(DEzdds2, contrast = list("PopulationSA.TreatmentZ"), alpha = 0.05)
resSA_c_z$Names <- rownames(resSA_c_z)
resBD_c_z <- results(DEzdds2, contrast = list("PopulationBD.TreatmentZ"), alpha = 0.05)
resBD_c_z$Names <- rownames(resBD_c_z)
resGR_c_z <- results(DEzdds2, contrast = list("PopulationGR.TreatmentZ"), alpha = 0.05)
resGR_c_z$Names <- rownames(resGR_c_z)
resPP_c_z <- results(DEzdds2, contrast = list("PopulationPP.TreatmentZ"), alpha = 0.05)
resPP_c_z$Names <- rownames(resPP_c_z)

#Summary of the results:
#could use these values for a results table of total numbers of loci up and down regulated in zinc vs control in each population.
summary(resSA_c_z)  #large stress response to zinc in coasts
summary(resBD_c_z) 
summary(resGR_c_z) 
summary(resPP_c_z) 

#check dimensions
dim(resSA_c_z)[1]
dim(resBD_c_z)[1]
dim(resGR_c_z)[1]
dim(resPP_c_z)[1]

#remove any NAs by taking only non NA data entries
resSA_c_z <- resSA_c_z[is.na(resSA_c_z$padj) == FALSE,]
resBD_c_z <- resBD_c_z[is.na(resBD_c_z$padj) == FALSE,]
resGR_c_z <- resGR_c_z[is.na(resGR_c_z$padj) == FALSE,]
resPP_c_z <- resPP_c_z[is.na(resPP_c_z$padj) == FALSE,]

#check dimensions
dim(resSA_c_z)[1]
dim(resBD_c_z)[1]
dim(resGR_c_z)[1]
dim(resPP_c_z)[1]

#filter p values for < 0.05 sig DE
resSA_c_z_sig <- resSA_c_z[resSA_c_z$padj < 0.05,]
resBD_c_z_sig <- resBD_c_z[resBD_c_z$padj < 0.05,]
resGR_c_z_sig <- resGR_c_z[resGR_c_z$padj < 0.05,]
resPP_c_z_sig <- resPP_c_z[resPP_c_z$padj < 0.05,]

#check dimensions
dim(resSA_c_z_sig)[1]
dim(resBD_c_z_sig)[1]
dim(resGR_c_z_sig)[1]
dim(resPP_c_z_sig)[1]

#filter p values for those > 0.05 no sig DE
resSA_c_z_noDE <- resSA_c_z[resSA_c_z$padj > 0.05,]
resBD_c_z_noDE <- resBD_c_z[resBD_c_z$padj > 0.05,]
resGR_c_z_noDE <- resGR_c_z[resGR_c_z$padj > 0.05,]
resPP_c_z_noDE <- resPP_c_z[resPP_c_z$padj > 0.05,]

#check dimensions
dim(resSA_c_z_noDE)[1]
dim(resBD_c_z_noDE)[1]
dim(resGR_c_z_noDE)[1]
dim(resPP_c_z_noDE)[1]

#write the results out
write.csv(resSA_c_z_sig, "resSA_c_z_sig_18_05_23.csv")
write.csv(resBD_c_z_sig, "resBD_c_z_sig_18_05_23.csv")
write.csv(resGR_c_z_sig, "resGR_c_z_sig_18_05_23.csv")
write.csv(resPP_c_z_sig, "resPP_c_z_sig_18_05_23.csv")

write.csv(resSA_c_z_noDE, "resSA_c_z_noDE_18_05_23.csv")
write.csv(resBD_c_z_noDE, "resBD_c_z_noDE_18_05_23.csv")
write.csv(resGR_c_z_noDE, "resGR_c_z_noDE_18_05_23.csv")
write.csv(resPP_c_z_noDE, "resPP_c_z_noDE_18_05_23.csv")


############################################################################-
#5. getting shrunken log fold change estimates for groups of interest ----
############################################################################-

#salt experiment coastal and mine population salt plasticity:

#estimate LFCs with LFCshrink:
# lfcSA_c_s <- lfcShrink(DEsdds2, contrast = list("PopulationSA.TreatmentS"), type = "ashr", alpha = 0.05)
# lfcSA_c_s$Names <- rownames(lfcSA_c_s)
# head(lfcSA_c_s)
# lfcBD_c_s <- lfcShrink(DEsdds2, contrast = list("PopulationBD.TreatmentS"), type = "ashr", alpha = 0.05)
# lfcBD_c_s$Names <- rownames(lfcBD_c_s)
# head(lfcBD_c_s)
# lfcGR_c_s <- lfcShrink(DEsdds2, contrast = list("PopulationGR.TreatmentS"), type = "ashr", alpha = 0.05)
# lfcGR_c_s$Names <- rownames(lfcGR_c_s)
# head(lfcGR_c_s)
# lfcPP_c_s <- lfcShrink(DEsdds2, contrast = list("PopulationPP.TreatmentS"), type = "ashr", alpha = 0.05)
# lfcPP_c_s$Names <- rownames(lfcPP_c_s)
# head(lfcPP_c_s)
# 
# #check certain gene membership:
# lfcSA_c_s[lfcSA_c_s$Names == "evm.TU.s_11645.12",]
# resSA_c_s[resSA_c_s$Names == "evm.TU.s_11645.12",]
# 
# #filter to remove FALSE results
# lfcSA_c_s <- lfcSA_c_s[is.na(lfcSA_c_s$padj) == FALSE,]
# lfcSA_c_s <- lfcSA_c_s[is.na(lfcSA_c_s$padj) == FALSE,]
# lfcSA_c_s <- lfcSA_c_s[is.na(lfcSA_c_s$padj) == FALSE,]
# lfcSA_c_s <- lfcSA_c_s[is.na(lfcSA_c_s$padj) == FALSE,]
# 
# #filter to remove those with p adj. below 0.05
# lfcSA_c_s_sig <- lfcSA_c_s[lfcSA_c_s$padj < 0.05,]
# lfcSA_c_s_sig <- lfcSA_c_s[lfcSA_c_s$padj < 0.05,]
# lfcSA_c_s_sig <- lfcSA_c_s[lfcSA_c_s$padj < 0.05,]
# lfcSA_c_s_sig <- lfcSA_c_s[lfcSA_c_s$padj < 0.05,]
#
#export results
# write.csv(lfcSA_c_s, "lfcSA_c_s_sig_20_06_23.csv")
# write.csv(lfcBD_c_s, "lfcBD_c_s_sig_20_06_23.csv")
# write.csv(lfcGR_c_s, "lfcGR_c_s_sig_20_06_23.csv")
# write.csv(lfcPP_c_s, "lfcPP_c_s_sig_20_06_23.csv")
#
#
#write.csv(lfcSA_c_s_sig, "lfcSA_c_s_sig_20_06_23.csv")


###########################################################################-
#6. PCAs for combined data and individual data ----
###########################################################################-

#########################################-
##6.1 transform data ready for PCAs ----
#########################################-

allvst1 <- vst(alldds1)
svst1 <- vst(sdds1)
#zvst1 <- vst(zdds1)

#subset vsts for certain combinations of treatments: 
allvst1_controls <- allvst1[,allvst1$Treatment %in% c("C1", "C2")]
allvst1_controls$Treatment
allvst1_treatments <- allvst1[,allvst1$Treatment %in% c("S", "Z")]
allvst1_treatments$Treatment

##############################################-
##6.2 generate pca data and pca axis labels ----
##############################################-

#overall dataset all groups PCA:
pca_data_allvst1 <- plotPCA(allvst1, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
pca_data_allvst1
pca_data_allvst1_controls <- plotPCA(allvst1_controls, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
pca_data_allvst1_controls
pca_data_allvst1_treatments <- plotPCA(allvst1_treatments, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
pca_data_allvst1_treatments

#zinc and salt PCA data for separate results:
pca_data_svst1 <- plotPCA(svst1, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
pca_data_svst1
#pca_data_zvst1 <- plotPCA(zvst1, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
#pca_data_zvst1

# #get percentage variance for each PCA (PCA1 and 2:
allvst1percentVar_all <- round(100 * attr(pca_data_allvst1, "percentVar"))
allvst1percentVar_C <- round(100 * attr(pca_data_allvst1_controls, "percentVar"))
allvst1percentVar_T <- round(100 * attr(pca_data_allvst1_treatments, "percentVar"))
svst1percentVar_all <- round(100 * attr(pca_data_svst1, "percentVar"))
#zvst1percentVar_all <- round(100 * attr(pca_data_zvst1, "percentVar"))

#check 3 and 4 for the controls:
pca_data_allvst1_controls_34 <- plotPCA.san(allvst1_controls, intgroup = c("Ecotype", "Location", "Treatment"), returnData=TRUE)
allvst1percentVar_C_34 <- round(100 * attr(pca_data_allvst1_controls_34, "percentVar"))



###################################-
##6.3 plot PCAs ----
###################################-

# #plot graph of PCA 1 and 2 for all data
# PCA_allvst1_all <- ggplot(pca_data_allvst1, aes(PC1, PC2, shape=Treatment, color=Location, fill=Ecotype)) +
#   geom_point(size = 2, stroke = 1) +
#   xlab(paste0("PC1:", allvst1percentVar_all[1], "% variance")) +
#   ylab(paste0("PC2:", allvst1percentVar_all[2], "% variance")) +
#   scale_shape_manual(values= c(21, 23, 22, 24)) +
#   scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
#   scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
#   #coord_fixed() +
#   theme(legend.position = "none")+
#   theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
#       axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 20), axis.text = element_text(size = 18),
#       legend.text = element_text(size = 18), legend.title = element_text(size = 20),
#       aspect.ratio = 1)
#PCA_allvst1_all

# #plot graph of PCA 1 and 2 for all controls
PCA_allvst1_C <- ggplot(pca_data_allvst1_controls, aes(PC1, PC2, shape=Treatment, 
                                                       color=Location, 
                                                       fill=Ecotype,
                                                        label=name)) +
  geom_point(size = 7, stroke = 1) +
  #geom_text(check_overlap = F, size = 2)+
  xlab(paste0("PC1:", allvst1percentVar_C[1], "% variance")) +
  ylab(paste0("PC2:", allvst1percentVar_C[2], "% variance")) +
  scale_shape_manual(values= c(21, 23)) +
  scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
  scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
  coord_fixed() +
  #theme(legend.position = "none")+
  theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
        axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.title = element_text(size = 22),
        aspect.ratio = 1)
PCA_allvst1_C


#axes 3 and 4. need to plot 3 against 1 ideally for extra info.
# PCA_allvst1_C_34 <- ggplot(pca_data_allvst1_controls_34, aes(PC3, PC4, shape=Treatment, color=Location, fill=Ecotype)) +
#   geom_point(size = 7, stroke = 1) +
#   xlab(paste0("PC1:", allvst1percentVar_C_34[1], "% variance")) +
#   ylab(paste0("PC2:", allvst1percentVar_C_34[2], "% variance")) +
#   scale_shape_manual(values= c(21, 23)) +
#   scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
#   scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
#   coord_fixed() +
#   #theme(legend.position = "none")+
#   theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
#         axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 18), axis.text = element_text(size = 18),
#         legend.text = element_text(size = 18), legend.title = element_text(size = 22),
#         aspect.ratio = 1)
# PCA_allvst1_C_34

PCA_allvst1_T <- ggplot(pca_data_allvst1_treatments, aes(PC1, PC2, shape=Treatment, color=Location, fill=Ecotype)) +
  geom_point(size = 7, stroke = 1) +
  xlab(paste0("PC1:", allvst1percentVar_T[1], "% variance")) +
  ylab(paste0("PC2:", allvst1percentVar_T[2], "% variance")) +
  scale_shape_manual(values= c(22, 24)) +
  scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
  scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
  coord_fixed() +
  #theme(legend.position = "none")+
  theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
        axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.title = element_text(size = 22),
        aspect.ratio = 1)
PCA_allvst1_T

PCA_svst1 <- ggplot(pca_data_svst1, aes(PC1, PC2, shape=Treatment, color=Location, fill=Ecotype)) +
  geom_point(size = 7, stroke = 1) +
  xlab(paste0("PC1:", svst1percentVar_all[1], "% variance")) +
  ylab(paste0("PC2:", svst1percentVar_all[2], "% variance")) +
  scale_shape_manual(values= c(21, 22)) +
  scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
  scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
  coord_fixed() +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(colour = "gray40", fill = "transparent"), panel.grid = element_blank(),
        axis.ticks.length = unit(.20, "cm"), axis.title = element_text(size = 22), axis.text = element_text(size = 18),
        legend.text = element_text(size = 18), legend.title = element_text(size = 22),
        aspect.ratio = 1)
PCA_svst1

# PCA_zvst1 <- ggplot(pca_data_zvst1, aes(PC1, PC2, shape=Treatment, color=Location, fill=Ecotype)) +
#   geom_point(size = 5, stroke = 1) +
#   xlab(paste0("PC1:", zvst1percentVar_all[1], "% variance")) +
#   ylab(paste0("PC2:", zvst1percentVar_all[2], "% variance")) +
#   scale_shape_manual(values= c(21, 24)) +
#   scale_fill_manual(values = c(col_coast, col_mine), guide = 'none') +
#   scale_color_manual(values = c(col_PPBD,  col_GRSA), guide = "none") +
#   #coord_fixed() +
#   theme(legend.position = "none")+
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), panel.background = element_blank(), aspect.ratio = 1, axis.title = element_text(size = 22), axis.text = element_text(size = 15))
# PCA_zvst1

####export the PCAs to PDF format:
# 
# #all PCA 
# pdf(file = "PCA_allvst1_all_18_05_23.pdf")
# PCA_allvst1_all
# dev.off()
#controls only
# pdf(file = "PCA_allvst1_C_18_05_23.pdf")
# PCA_allvst1_C
# dev.off()

# pdf(file = "PCA_allvst1_C_27_02_24.pdf")
# PCA_allvst1_C
# dev.off()
# # #treatments only
# pdf(file = "PCA_allvst1_T_27_02_24.pdf")
# PCA_allvst1_T
# dev.off()
# # #salt exp only
# pdf(file = "PCA_svst1_27_02_24.pdf")
# PCA_svst1
# dev.off()

setFplot_page(page = "a4", margins = "normal", units = "tw", pt = 20, reset = FALSE, w2h = 1)

# pdf_fit(file = "PCA_svst1_27_02_24_v3.pdf", pt =20, width = 1.15, w2h = 1)
# PCA_svst1
# fit.off()
# 
# pdf_fit(file = "PCA_svst1_27_02_24_v3.pdf", pt =20, width = 1.15, w2h = 1)
# PCA_svst1
# fit.off()

# #zinc exp only
# pdf(file = "PCA_zvst1_18_05_23.pdf")
# PCA_zvst1
# dev.off()


####################################################################################-
#7. Plot the normalised counts for all the genes in the dataset ----
###############################################################################-

#run DEseq estimate size factors
alldds1_SF <- estimateSizeFactors(alldds1)
alldds1_norm_counts <- counts(alldds1_SF, normalized=TRUE)

#check the counts to see what they look like :)
dim(alldds1_norm_counts)[1] #30178 dimension
str(alldds1_norm_counts)
hist(rowSums((alldds1_norm_counts)), breaks = 5000, xlim = c(0,20000))

#30178 genes are filtered in both and common across both comparisons
#write to a file:
write.csv(alldds1_norm_counts, "alldds1_norm_counts_18_05_23.csv")
# 
# #####
# #7.0checking some genes for counts ----
# 
# plotCounts(alldds1, gene = "evm.TU.s_36.759", intgroup = "Pop_treat")
# plotCounts(alldds1, gene = "evm.TU.s_49.16", intgroup = "Pop_treat")
# 
# #evm.TU.s_49.16



