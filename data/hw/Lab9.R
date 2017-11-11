##############################################################
# Name:
# CSC-315
# Lab #9: Limma, heatmaps, and analyzing processed GEO data
#############################################################

##########################################################################
# Add R code to the script below and create a Notebook which to complete
# the steps and explicitly answer the following questions
##########################################################################

library(limma)
library(ggplot2)
library(dplyr)

##################################################################################
# 1.The code below reads in the survey data and performs a 
#   2-sample t-test to evaluate whether there is a statistically
#   significant difference in alcohol consumption between 'Cat' vs. 'Dog'
#   people. Based on the code below, (a) find the p-value and state 
#   your conclusion regarding the null hypothesis of H0: mu_cat - mu_dog = 0;
#   and (b) calculate the difference in mean Alcohol consumption between
#   groups, using the formula: 
#     mean consumption of dog people - mean consumption of cat people
##################################################################################
survey <- read.csv("https://gdancik.github.io/CSC-315/data/datasets/survey_fall-2017.csv")
s <- split(survey$Alcohol, survey$CatOrDogPerson)
res <- t.test(s$Cat, s$Dog, var.equal = TRUE)


##################################################################################
# 2.Fit a linear model that predicts Alchohol consumption based on 
#   whether an individual is a cat or a dog person. You should use
#   the treatment contrast where 'cat person' is the reference (x = 0) and 
#   'dog person' is the treatment (x = +1)
#    
# (a) Find and interpret the y-intercept of the regression line in the
#      context of this problem.
# (b) Find and interpret the slope in the context of this problem
# (c) What is the p-value for the hypothesis test that there is a
#     significant difference in Alcohol consumption between the two groups?
#     (show this result in R, based on the linear model)
##################################################################################


###############################################################
# 3. Get the processed data for GSE19143 and pull out the 
#    expression data and phenotype data. Note that this
#    dataset contains gene expression samples from children
#    with Acute Lymphoblastic Leukemia (ALL), a cancer of
#    the bone marrow. Tumor samples were treated with
#    the anti-inflammatory drug prednisolone, and determined 
#    to be either sensitive (responsive) or resistant 
#    (non-responsive) to this drug. 
###############################################################

# (a) How many samples had their gene expression values profiled?

# (b) How many probes are on the array?



#####################################################################
# 4.How many individuals are resistant to prednisolone and
# how many are sensitive? 
#####################################################################

#####################################################################
# 5. Find the top differentially expressed probes, with a FDR of 10%,
# between individuals that are resistant vs. sensitive to prednisolone.
# Note: there should be 16 probes total. How many of these probes 
# are up-regulated (i.e., have higher expression) in resistant 
# individuals and how many are down-regulated (i.e., have lower 
# expression) in resistant individuals. 
#####################################################################

########################################################################
# 6. Construct a heatmap of these top 16 probes, with individuals 
# color-coded by response to prednisolone (with green=sensitive and 
# red = resistant). (Note: if you are unable to complete question 5), 
# you may do this with the first 16 probes in the expression matrix).
########################################################################

########################################################################
# 7. If you answered question 5 correctly, the SECOND hit 
# should be for the probe 209374_s_at. Show that this probe
# corresponds to the gene IGHM, by first downloading the 
# correct platform data from GEO, and then finding the gene
# associated with this probe. 
#######################################################################



#####################################################################
# 8. How many probes are there for the gene IGHM on the platform
# in this study? Note: you must search for this gene using the
# regular expressions covered in the GEO-and-limma.R script. Your 
# code must also output the number of probes. 
####################################################################



########################################################################
# Final Notes: the heatmap in question 6 provides a candidate list
# of probes associated with prednisolone response in children with 
# leukemia. Although much additional work and testing needs to be done, 
# this kind of gene signature could ultimately be used to determine whether a 
# child with leukemia would benefit from prednisolone treatment, or 
# whether an alternative treatment might be more effective.

# The IGHM finding is also interesting. IGHM is a gene that codes
# for an antibody protein involved in the immune reponse; the 
# fact that this gene is differentially expressed beween responders and 
# non-respnoders suggests that a patient's immune system may play a
# role in how they respond to prednisolone)
########################################################################