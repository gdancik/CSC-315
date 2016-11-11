##############################################################
# Name:
# CSC-315
# Lab #9: Limma, heatmaps, and analyzing processed GEO data
#############################################################

##########################################################################
# Add R code to the script below and create a Notebook which to complete
# the steps and explicitly answer the following questions
##########################################################################


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
survey = read.delim("http://pastebin.com/raw/QDSga7qF")
s = split(survey$Alcohol, survey$CatOrDog)
cat = s[[1]] # alcohol consumption of cat people
dog = s[[2]] ## alcohol consumption of dog people
res = t.test(cat, dog, var.equal = TRUE)


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

#    How many samples had their gene expression values profiled?

#    How many probes are on the array?


##############################################################
# 4. From the clinical data,  how many males are there and 
#    how many females are there? 
##############################################################


#####################################################################
# 5. We have seen previously that the gene XIST is associated 
# with sex. How many probes for XIST are on the platform used
# in this study. This can be found usng the hgu133a.db library.
# Construct a boxplot and find the p-value corresponding to 
# whether the first probe for this gene is differentially expressed 
# across males and females in this study 
####################################################################


#####################################################################
# 6.How many individuals are resistant to prednisolone and
# how many are sensitive? (Note: by coincidence, the numbers
# are the same as the numbers in (4)).
#####################################################################

#####################################################################
# 7. Find the top 10 differentially expressed probes between
# individuals that are resistant vs. sensitive to prednisolone.
# How many of these probes are up-regulated (i.e., have higher
# expression) in resistant individuals and how many are down-
# regulated (i.e., have lower expression). 
# What percentage of these probes do you expect to be false 
# positives? (Note: this is the adjusted p value for the last probe
# in the list)
#####################################################################

########################################################################
#8. Construct a heatmap of these top 10 probes, with individuals 
# color-coded by response to prednisolone (with green=sensitive and 
# red = resistant). 
########################################################################

########################################################################
# 9. If you answered question 7 correctly, the SECOND hit 
# should be for the probe 209374_s_at. What gene does this probe 
# correspond to? According to genecards, what is the "Molecular Function"
# of this gene?  
#######################################################################


########################################################################
# Final note: the heatmap in question 8 provides a candidate list
# of genes for predicting prednisolone response in ALL 
# patients. Although much additional work and testing needs to be
# done, a similar gene signature could ultimately be used to determine
# whether a child with leukemia would benefit from prednisolone
# treatment, or whether an alternative treatment might be more 
# effective.
########################################################################