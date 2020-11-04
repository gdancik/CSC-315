################################################################
# Lab 8: GEO Lab - Raw Data
# In this lab you will analyze two probes from a gene expression
# study of Alzheimer's Disease (AD). The dataset is
# available from: 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1297
################################################################

library(affy)
library(dplyr)
library(ggplot2)

#####################################################################
# 1) Download the raw data (GSE1297_RAW.tar file) and extract the 
# files. Read in the data using the ReadAffy() function after 
# setting the celfile.path argument appropriately below
################################################################

GSE1297 <- ReadAffy(celfile.path = "/Users/dancikg/Downloads/GSE1297_RAW")

###################################################################
# 2) Process the gene expression data using the Robust Multi-Array 
# Average (RMA) method and extract the expression data.
###################################################################

# Use the appropriate R function on the expression matrix to 

# a) find the number of probes 


# b) find the number of samples

# c) generate a boxplot of the expression values for each sample
#    to confirm that the data has been normalized
###################################################################

###################################################################
# We will see how to get the phenotype data from GEO in a later
# class; for this lab, the data has already been processed
# and can be read in using the statement below. The data includes
# the following information (though we will not use these):
# MMSE.Score = mini–mental state examination score for 
#   cognitive impairment (low scores indicate impairment)
# NFT.Score = protein markers for AD
################################################################

GSE1297.p <- read.delim("http://bioinformatics.easternct.edu/BCBET2/GSE1297.p.xlsx")

########################################################################
# The code below adds a column to the data frame that includes the
# Alzheimer's Disease status for each sample, which it gets by deleting
# the sample number (e.g., "Control 1008" becomes "Control". This is 
# done by using a regular expression, which removes everything in 
# the string beginning with the space). Note that your pheno data
# must be stored in GSE1297.p for the code to work.
########################################################################

GSE1297.p <- GSE1297.p %>%  
    mutate(AD.status = gsub(" .*", "", Sample)  )

########################################################################
# The code below generates a boxplot of MMSE scores across
# the different levels of AD severity
########################################################################

ggplot(GSE1297.p, aes(AD.status, MMSE.Score, fill = AD.status)) +
  geom_boxplot() +
  theme_classic() +
  guides(fill = "none") +
  ggtitle("Relationship between MMSE score and AD severity") 

#####################################################
# 3) A gene called APOE is associated with late onset
#    Alzheimer's disease. One of the probes for 
#    APOE is 203381_s_at, and one of the probe set 
#    sequences for this probe is
#   'AGGCCAAGGTGGAGCAAGCGGTGGA'. For this question,
#   we will focus on the first 5 nucleotides, 
#   'AGGCC'. During the gene expression profiling
#   process, if mRNA containing 'AGGCC' exists in
#   the sample, it will be converted to cDNA, which will
#   hybridize to the probe. What is the cDNA sequence
#   that will hybridize (bind) to 'AGGCC'?
#####################################################

############################################################
# 4) Construct side-by-side boxplots showing the expression
# of the probe 203381_s_at for CONTROL patients and 
# patients with SEVERE AD. (The boxplot must be constructed
# using ggplot, and should only include Control and Severe
# patients, since these are the groups we want to compare)
############################################################


###########################################################
# Perform a two sample t-test to evaluate whether or not
# expression is significantly different between CONTROL
# patients and patients with SEVERE AD. Report the 
# fold change and the p-value and state your conclusion.
###########################################################


#####################################################
## Repeat the boxplot and t.test for the gene PSEN1 
##  using the probe 207782_s_at
#####################################################

########################################################
## Based on the above analyses, what is your conclusion
## about the association between the genes APOE and
## PSEN1 and Alzheimer's Disease / cognitive 
## impairment?
###############################################s######

###########################################################
## If you are interested, more information about
## Alzheimer's Disease and these genes can be found
## at: http://ghr.nlm.nih.gov/condition/alzheimer-disease
##########################################################
