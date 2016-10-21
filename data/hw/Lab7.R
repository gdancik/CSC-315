####################################################################
# Name:
# Lab 7: Hypothesis testing for population means.
# For these questions, we will assume that the Central 
# Limit Theorem applies (that the populations are normally 
# distributed or n is sufficiently large) and that the samples 
# are representative of the population of interest. 
####################################################################

# Read in our survey data
survey = read.delim("http://pastebin.com/raw/QDSga7qF")

##########################################################################
# 1. Is there evidence that a student's college GPA differs from his/her
# high school GPA? We can test this by evaluating whether or
# not the mean difference between GPAs differs from 0. 
diff = survey$College.GPA-survey$HS.GPA

# (a) State the null and alternative hypotheses
# (b) Calculate / find the test statistic (and specify the degrees of freedom)
# (c) Find the p-value using t.test
# (d) Find the p-value 'manually' based on the test statistic and
#     appropriate degrees of freedom (if calclating the df manually,
#     be careful for missing values)
# (e) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.
##########################################################################

##########################################################################
# 2. Is there evidence that the college GPA of a 'cat' person 
# differs from that of a 'dog' person?
# (a) State the null and alternative hypotheses
# (b) Create side-by-side boxplots showing College GPA for 'cat' and
#     'dog' people. Make sure to label the y-axis and give the 
#     chart a title.
# (c) We will now formally test the hypotheses that mean college GPA
#     is different between 'cat' and 'dog' people. The command 
#     t.test(x,y) will perform a two-sample t-test for the null 
#     hypothesis that the 'x' and 'y' populations have the same mean.
#     Use the t.test function to find the test statistic and the 
#     degrees of freedom
# (d) Find the p-value using t.test
# (e) Find the p-value 'manually' based on the test statistic and
#     appropriate degrees of freedom
# (f) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.
# (g) What would it mean in the context of this problem if a Type I 
#     error occurred?
##########################################################################

##########################################################################
# 3:  Find the p-values associated with the following t test 
# statistics, for a one-sample t-test, and state whether you would reject 
# or fail to reject the null hypothesis at α=0.05:
# (a) t = 2.78, n = 45
# (b) t = -3.3, n = 51
# (c) t = 1.11, n = 100
##########################################################################


# use the cereal data to complete the last question
cereal = read.delim("http://pastebin.com/raw/0G6DrHyC")

#################################################################################
# 4:  The 'sugars' column contains the sugar content (in grams), while
#     the 'shelf' column contains the shelf in which the cereal is 
#     shelved on, with 1 = lower shelf, 2 = middle shelf (which is at
#     eye level for children), and 3 = top shelf.
#     You will use this data to test whether mean sugar content differs
#     between the lower shelf and the middle shelf.
# (a) State the null and alternative hypotheses
# (b) Calculate / find the test statistic (and specify the degrees of freedom)
# (c) Find the p-value
# (d) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.
#################################################################################

