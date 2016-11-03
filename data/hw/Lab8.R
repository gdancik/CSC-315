################################################################
# Lab 8: GEO Lab - Raw Data
# In this lab you will analyze two probes from a gene expression
# study of Alzheimer's Disease (AD). The dataset is
# available from: 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1297
################################################################

library(affy)

#####################################################################
# Download the raw data (GSE1297_RAW.tar file) and extract the files  
# Read GSE1297 using ReadAffy() and the celfile.path argument
# which must be set appropriately below
################################################################

GSE1297 = ReadAffy(celfile.path = "/Users/dancikg/Downloads/GSE1297_RAW/")

################################################################
# Process the gene expression data using the Robust Multi-Array 
# Average (RMA) method and extract the expression data.
# How many probes and samples does this dataset contain?
# Generate a boxplot of the expression values of each sample
# to confirm that the data has been normalized
################################################################


################################################################
# Read the phenotype data:
# MMSE.Score = mini–mental state examination score for 
#   cognitive impairment (low scores indicate impairment)
# NFT.Score = protein markers for AD
################################################################

GSE1297.p = read.delim("http://www1.easternct.edu/dancikg/files/2015/10/GSE1297.p.xlsx")


##################################################################
# The code below gets the group names from the sample names of the 
# pheno table and constructs a scatterplot of MMSE and NFT
# scores with points color-coded by AD severity. You should 
# understand how to use the strsplit function and the how the 
# colors are produced since they will be needed later 
# (note: code assumes that the pheno table is stored in GSE1297.p)
###################################################################

##################################################
# get group names from sample names, which have
# format Group SampleNumber
##################################################
sample.names = as.character(GSE1297.p$Sample)
## split each row name based on a blank space ##
s = strsplit(sample.names, " ")
# function to return the 1st element of a list
get.first <- function(x) {return(x[[1]])}  
# apply this function to get the group name
groups = sapply(s, get.first)

# set colors corresponding to each level: 
# CTL, Incipient, Moderate, Severe
colors = c("darkblue", "orange", "purple", "darkred")
group.num = as.integer(factor(groups))

# col vector contains appropriate color for each sample
col = colors[group.num]

## plot MMSE score and NFT score, colored by severity ##
plot(GSE1297.p$MMSE.Score, GSE1297.p$NFT.Score, 
     main = "MMSE and NFT scores in AD patients",
     xlab = "MMSE Score", ylab= "NFT Score",
     pch = 19, col = col)
l = levels(factor(groups)) # possible values
legend("right", pch = 19, legend = l, col = colors)

########################################################
# Describe the relationship between MMSE and NFT score.
# Would you expect a person with a high MMSE score to
# have Alzheimer's Disease?
#####################################################

#####################################################
# A gene called APOE is associated with late onset
# Alzheimer's disease. One of the probes for 
# APOE is 203381_s_at
#####################################################

###########################################################
# Construct side-by-side boxplots showing the expression
# of the probe 203381_s_at across the 4 patient groups,
# using the color coding from above.
# Perform a two sample t-test to evaluate whether or not
# expression is significantly different between CONTROL
# patients and patients with SEVERE AD. Report the 
# fold change and the p-value and state your conclusion.
###########################################################

#####################################################
# Construct a scatterplot of gene expression of
# the probe 203381_s_at on the x-axis and MMSE score
# on the y-axis, using the color coding above. Give 
# the graph an appropriate title, and axis labels, and 
# legend, and also add the regression line. What is 
# the correlation between MMSE score and expression? 
#####################################################

#####################################################
# The cor.test function can be used to evaluate the
# null hypothesis that the correlation between two
# variables is equal to 0. The function is called
# using cor.test(x,y), where x and y are the vectors 
# of observations. Find the p-value, report the 
# correlation, and state whether or not the 
# correlation between expression and MMSE score is 
# statistically significantly 
#####################################################


#####################################################
## Repeat the boxplot, t.test, scatterplot, and 
## cor.test for the gene PSEN1 using the probe 
## 207782_s_at
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
