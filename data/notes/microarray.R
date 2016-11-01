#############################################
# Microarray Basics
#############################################

library(affy)

############################################################
# Look at the Dilution dataset: this contains 4 liver
# tissue hybridized in concentrationsof 10 and 20 micrograms
# to scanner 1(A) and scanner 2(B). Note that Dilution is an
# object of type AffyBatch (or an extension of an 
# ExpressionSet (eSet) object)
######################################################
library(affydata)
data(Dilution)
Dilution

####################################################
# Methods (functions) for AffyBatch objects
####################################################
sampleNames(Dilution)     # the names of the samples
experimentData(Dilution)  # experiment information
annotation(Dilution)      # annotation (the microarray used)
phenoData(Dilution)       # get phenotypic (clinical) data
pData(Dilution) # phenotypic data in table form
varMetadata(Dilution) # description of phenotypic data

####################################################
# Let's look at each microarray
####################################################
image(Dilution, col = heat.colors(500))

####################################################
# Produces a boxplot of log base 2 intensities # 
# Does it seem fair to compare gene expression 
# across concentrations or scanners?
####################################################
boxplot(Dilution, col = 1:4, ylab = "log2 expression")

####################################################
# Normalization of microarray data involves
# 1. Background correction to remove noise
# 2. Normalization of samples
# 3. Estimation of "average" probe intensity values
####################################################

################################################################
# We will use Robust Multi-array average (RMA) which involves 
# 1. Background correction of probe level intensity values
# 2. Quantile Normalization
# 3. Estimation of average probe set values on log2 scale
################################################################

Dilution.rma = rma(Dilution)        # perform RMA
Dilution.expr = exprs(Dilution.rma) # extract the expression values

################################################################
# confirm that samples are quantile normalized
################################################################
boxplot(Dilution.expr, col = 2:5, ylab = "log2 expression",
        main = "Small part of dilution study",
        xlab = "Sample")


################################################################
# Let's look at a leukemiadataset, and compare 
# acute lymphoblastic leukemia (ALL) samples with healthy, 
# non-leukemia (noL) bone marrow samples
################################################################
library(leukemiasEset)
data(leukemiasEset)

################################################################
# Extract phenotype data. How many samples of each leukemia 
# type is there?
################################################################
leukemia.p = pData(leukemiasEset)


################################################################
# This data is already processed using RMA, so it is sufficient
# to extract the expression data
################################################################
leukemia.expr = exprs(leukemiasEset) 

# confirm that data has been processed #
boxplot(leukemia.expr, main = "Leukemia samples", ylab = "log2 expression")


######################################################################
# Let's compare expression between ALL and healthy bone marrow 
# samples for for the following probes:

# ENSG00000171960 - corresponds to gene PPIH 
#       (https://www.ncbi.nlm.nih.gov/gene/10465)
# ENSG00000135679 - corresponds to gene MDM2
#       (https://www.ncbi.nlm.nih.gov/gene/4193)

# We want to calculuate the fold change (see below) and the p-value
# evaluating whether or not any difference in expression is statistically 
# significant
######################################################################

# find expression for desired probe
m = match("ENSG00000171960", rownames(leukemia.expr))
s.all = split(leukemia.expr[m,], leukemia.p$LeukemiaType)
boxplot(s.all)

# let's only look at ALL and NoL
s = list(ALL = s.all$ALL, Normal = s.all$NoL)
boxplot(s, ylab = "log2 expression",
        main = "PPIH (ENSG00000171960) expression",
        col = c("darkred", "darkblue"))


########################################################################
# Calculate fold change (FC) which is average expression in the 
# first group divided by the average expression in the second group,
# Since data is on the log2 scale, we must convert back to normal scale
########################################################################

l = lapply(s, mean)
fc = 2**abs(l[[1]]-l[[2]])  ## data is on log2 scale
main = paste("PPIH (ENSG00000171960) expression,\nFC = ", round(fc,2))
boxplot(s, col = c("darkred", "darkblue"), main = main)

########################################################################
# Is the difference in fold change statistically significant??
# H0: mu_ALL - mu_normal = 0
# HA: mu_ALL - mu_normal != 0
# where mu_ALL is the mean expression of ALL samples for the 
# probe of interest and mu_normal is the mean expression of
# normal samples
# How do we carry out a hypothesis test comparing two population means?
########################################################################

########################################################################
# Repeat analysis for probe ENSG00000135679 (MDM2)
########################################################################
m = match("ENSG00000135679", rownames(leukemia.expr))
s = split(leukemia.expr[m,], leukemia.p$LeukemiaType)

