##############################################################
# Name:
# CSC-315
# Lab #9: Limma, heatmaps, and analyzing processed GEO data
#############################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps below and to explicitly answer the following questions
##########################################################################

library(limma)
library(GEOquery)
library(ggplot2)
library(dplyr)

# 1) The code below reads in our class survey data and performs a 
#   2-sample t-test to evaluate whether there is a statistically
#   significant difference in the hours of sleep between 
#   'Cat' and 'Dog' people. Based on the code below, 

survey <- read.csv("https://gdancik.github.io/CSC-315/data/datasets/CSC315_survey_Fall_2020.csv")
s <- split(survey$Sleep, survey$CatOrDogPerson)
res <- t.test(s$Cat, s$Dog, var.equal = TRUE)

#   (a) find the p-value and state your conclusion regarding the null 
#       hypothesis of H0: mu_cat - mu_dog = 0


#   (b) calculate the difference in mean hours of sleep between
#   groups, using the formula: 
#     mean hours of sleep for dog people - mean hours of sleep for cat people


# 2) Fit a linear model that predicts Hours of Sleep based on 
#   whether an individual is a cat or a dog person. You should use
#   the treatment contrast where 'cat person' is the reference (x = 0) and 
#   'dog person' is the treatment (x = +1)

  
# (a) Find and interpret the y-intercept of the regression line in the
#      context of this problem.


# (b) Find and interpret the slope of the regression line in the context of 
#     this problem


# (c) What is the p-value for the hypothesis test that there is a
#     significant difference in Hours of Sleep between the two groups?
#     (show this result in R, based on the linear model). Note: the 
#     p-value from the linear model should match the p-value from the
#     two-sample t-test from problem 1(a) above.

# 3) We can also fit a linear model using the 'sum' contrasts, which 
#    is done below. This model has the form:
#
#    y = a + bx, where x = 1 for a Cat person and -1 for a Dog person

fit <- lm(Sleep ~ CatOrDogPerson, data = survey, 
          contrasts = list(CatOrDogPerson = 'contr.sum'))

#  a) Show that the y-intercept is the average (mean) of the 
#     the group means, i.e.,the y-intercept is equal to 
#     [mean(sleep_cat_person) + mean(sleep_dog_person)] / 2
#     by calculating and displaying this value, and also 
#     displaying the y-intercept.

#  b) Find and interpret the slope of the regression line in the 
#     context of this problem


# 4) Get the processed data for GSE19143 and extract the 
#    expression data and phenotype data. Note that this
#    dataset contains gene expression samples from children
#    with Acute Lymphoblastic Leukemia (ALL), a cancer of
#    the bone marrow. Tumor samples were treated with
#    the anti-inflammatory drug prednisolone, and were
#    determined to be either sensitive (responsive) or 
#    resistant (non-responsive) to this drug. 

# (a) Use the expression matrix to find the number of samples that 
#     had their gene expression values profiled.


# (b) Use the expression matrix to find the number of probes on the array.


# (c) Take the log2 of the expression data, and generate a boxplot
#     to show that the samples are properly processed and normalized.
#     The analysis beginning with question 6 must use the log2 data; 
#     otherwise the results will not be correct.


# 5) How many individuals are resistant to prednisolone and
# how many are sensitive? (Hint: you will need to determine
# which column contains this information; it is one of the last
# columns)


# 6) Find the most differentially expressed probes, using an FDR of 
# 10%, for probes that are differentially expressed between 
# individuals who are resistant and individuals who are sensitive 
# to prednisolone. Note: there should be 16 probes total. 

# (a) How many of these probes are up-regulated (i.e., have higher
#     expression) in resistant individuals 

# (b) How many are down-regulated (i.e., have lower expression) 
#     in resistant individuals. 

# (c) How many of these probes do you expect to be false positives?


# 7) Construct a heatmap of these 16 probes, with individuals 
# color-coded by response to prednisolone (with green=sensitive and 
# red = resistant). For visualizing gene expression values, it is 
# common to use blue for high expression and yellow for low expression,
# but you may use any color combination that you wish. Note: if you 
# are unable to complete question 5), you may do this with the first 
# 16 probes in the expression matrix.

# 8) If you answered question 6 correctly, the SECOND hit 
# should be for the probe 209374_s_at. Show that this probe
# corresponds to the gene IGHM, by first downloading the 
# correct platform data from GEO, and then finding the gene
# associated with this probe. 


# 9) How many probes are there for the gene IGHM on the platform
# in this study? Note: you must search for this gene using the
# regular expressions covered in the GEO-and-limma.R script. Your 
# code must also output the number of probes. 



########################################################################
# Final Notes: the heatmap in question 6 provides a candidate list
# of probes associated with prednisolone response in children with 
# leukemia. Although much additional work and testing would need to be 
# done to validate these findings, this kind of gene signature could 
# ultimately be used to determine whether a child with leukemia would 
# benefit from prednisolone treatment, or whether an alternative 
# treatment might be more effective.

# The IGHM finding is also interesting. IGHM is a gene that codes
# for an antibody protein that is involved in the immune reponse; the 
# fact that this gene is differentially expressed between responders 
# and non-responders suggests that a patient's immune response may 
# play a role in how well they respond to prednisolone.
########################################################################
