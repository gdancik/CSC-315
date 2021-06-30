########################################################
# Name:
# CSC-315
# Lab #6: Hypothesis Tests -- Proportions 
########################################################
  
##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions
##########################################################################

# 1) A CBS news poll that surveyed 2,226 registered voters towards the 
#    end of August 2020 found that 1,158 (52%) answered 'Yes' to the question 
#    "Are you better off today than you were 4 years ago". Let p.hat be the
#    proportion of all registered voters that would answer 'Yes' to that question.
#    Link: # https://www.cbsnews.com/news/republicans-economy-coronavirus-opinion-poll-cbs-news-battleground-tracker/

# a) What is the mean and standard deviation of the distribution of p.hat
#    under the null hypothesis that p = 0.50. 

# b) Graph the distribution of p.hat and draw a vertical line at 
#    p.hat = 1158/2226, the proportion answering 'Yes' to the question.

# c) Calculate the z test statistic and graph its distribution under the 
#    null hypothesis as was done in class, drawing a vertical line at the 
#    z statistic. Find the p-value based on this test statistic. The 
#    z statistic should be around 1.907568 and the p-value should be about 
#    0.056557.

# d) Use the prop.test function to conduct the hypothesis test without
#    the continuity correction. Calculate the z test statistic from the 
#    prop.test object and extract the p-value (Note: these should match the 
#    test statistic and p-value from parts (b) and (c).
                                          
# e) Your p-value should be about 0.0564. Based on this p-value, state the 
#    conclusion regarding whether registered voters believe that
#    they are better off today than they were 4 years ago, using a level
#    of significance value of 0.05. Is there evidence that a majority agree 
#    or disagree with the statement?


# 2) A person who claims to be a psychic says that she can correctly 
#    predict the outcome of a roll of a standard die (with numbers 1-6) 
#    more times than what would be expected by chance. When you roll a 
#    die 50 times, she correctly predicts the outcome 12 times. 
#    Use the prop.test function to complete (b) and (c). 

# a) State the null and alternative hypothesis corresponding to this claim.

#       H0:
#       HA:

# b) Find the z test statistic 


# c) Find the p-value


# d) The p-value should be around 0.2295. Using this p-value, state the conclusion
#    regarding whether or not the person's ability to predict the outcome of the 
#    die is different than what would be expected by chance.


# 4) Find the p-values associated with the following z test statistics, and 
#    state whether you would reject or fail to reject the null hypothesis 
#    at alpha = 0.05.

# a)   z = -1.15	

# b)   z = 2.94

# c)   z = 1.05


# 5) For question (2), the null hypothesis is that the "psychic" does not do
#    better or worse than random chance; the alternative hypothesis is that
#    the psychic's predictive abilities are different than what is expected
#    from random chance. Complete the following to specify what would it mean
#    if a Type I or Type II error occured.

#    A Type I error means that we conclude that ___________________________
#    but in reality, _______________________________.

#    A Type II error means that we conclude that ___________________________
#    but in reality, _______________________________.

