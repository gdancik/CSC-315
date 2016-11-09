################################################################
# Relationship between the t-test and linear regression
################################################################

## We will compare the GPAs of drinkers to non-drinkers in our class ##
survey = read.delim("http://pastebin.com/raw/QDSga7qF")
s = split(survey$College.GPA, survey$Alcohol > 0)
nondrinkers = s[[1]] ## 1st group corresponds to survey$Alcohol > 0 being FALSE
drinkers = s[[2]] ## 2nd group corresponds to survey$Alcohol > 0 being TRUE

## test for difference in means, using two-sample t-test 
## by default, this does not assume equal variances between groups
t.test(nondrinkers, drinkers) 

## we can assume equal variances
result = t.test(nondrinkers, drinkers, var.equal = TRUE) 

####################################################
# let's formulate this in terms of a linear model 
# We want to predict GPA based on drinking status
####################################################
drinker = as.integer(survey$Alcohol > 0)
plot(drinker, survey$College.GPA, pch = 19, 
     xlab = "Drinker (0 = No, 1 = Yes)",
     ylab = "College GPA")

fit1 = lm(survey$College.GPA ~ drinker)
abline(fit1, col = "red")

#################################################
# Now compare p-values of t-test to the p-values
# from the linear model which tests against
# H0: b = 0 in y = a + bx
# HA: b != 0
#################################################
result$p.value
summary(fit1)

#################################################
## These analyses are the same! P-values are the
## same and the slope of the linear model is
## equal to the difference between means
#################################################

# compare difference in means
diff(result$estimate)
fit1$coefficients

########################################################
## There are several ways of formulating a linear model
## for this problem, in terms of the explanatory
## variable (we need one value for drinker and one
## value for non-drinker)
#######################################################

drinking.status = factor(survey$Alcohol > 0)
levels(drinking.status) = c("NonDrinker", "Drinker")
groups = levels(drinking.status)

##########################################################
# Treatment contrasts: explanatory variable is 0 for the
# reference group and 1 for the second group
##########################################################

# NonDrinker (reference) = 0, Drinker = 1
contr.treatment(groups)

# note that the explanatory variable (drinking.status) is a 
# factor, and we specify how to define the contrasts
fit.treatment = lm(survey$College.GPA ~ drinking.status,
                    contrasts = list(drinking.status = "contr.treatment"))

summary(fit.treatment)

## the slope is equal to the difference in means
coefficients(fit.treatment)[2]

######################################################
# Sum contrasts: explanatory variable is -1 for the 
# first group and +1 for the second group (these 
# values sum to 0)
######################################################

# NonDrinker (1st level) = 1, Drinker (2nd level) = -1
contr.sum(groups)
fit.sum = lm(survey$College.GPA ~ drinking.status,
                   contrasts = list(drinking.status = "contr.sum"))

# What does the slope represent? And why? #
coefficients(fit.sum)[2]

###################################################
# Use indicator (1 or 0) for group1 AND indicator 
# (1 or 0) for group2. Note that this model has 2 
# 'slope' variables and no 'intercept', i.e.,
# y = b1x1 + b2x2, where b1 is the mean for group 1,
# and b2 is the mean for group 2
###################################################
fit.indicator = lm(survey$College.GPA ~ -1 + drinking.status)
coefficients(fit.indicator)

#####################################################
# The p-values for each coefficient (b) test against
# H0: b = 0, but this is not of interest to us.
# The hypothesis test of interest is whether 
# b1-b2 = 0, which can be calculated manually or by 
# using additional packages (we will not worry about
# this for now). This linear model formulation
# is what we will use identify differentially expressed
# probes
#####################################################

