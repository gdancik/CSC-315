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

# b) Graph the distribution ofp.hat and draw a vertical line at 
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


# 5) In May of 2020, a survey of 10,957 American adults found that 72% of them
#    (7,889) would "definitely" or "probably" get a coronavirus vaccine if one 
#    were available. In September, a similar survey of 10,093 American
#    adults found that only 51% (5,147) would get a vaccine.

#    Ref: https://www.pewresearch.org/science/2020/09/17/u-s-public-now-divided-over-whether-to-get-covid-19-vaccine/

#    For this question, you will carry out a hypothesis test to evaluate whether
#    the proportion of American adults willing to get a vaccine has changed.
#    Note that the null and alternative hypotheses are as follows:

#    H0: pMay - pSept = 0
#    H1: pMay - pSept != 0

#   where pMay and pSept are the proportion of American adults, in May and September,
#   who would "definitely" or "probably" get a coronavirus vaccine.

# a) Show that the difference between population proportions is about 0.210039

# b) Show that the estimate of the common population proportion 
#     (i.e., the value of p on slide 17) is about 0.6193.

# c) Show that the standard deviation of the difference between the two 
#    proportions (i.e., the value of the denominator of Z on slide 17)
#    is approximately 0.006699.
  
# d) Calculate the z test statistic by dividing your answer to (b) by your answer
#    to (c)

# e) Find the z test statistic by carrying out the hypothesis test WITHOUT the
#    continuity correction, and calculate the z test statistic from the prop.test
#    object (this should match your answer to (d)

# f) To be a little more accurate, use prop.test to find the p-value WITH the
#    continuity correction. What is the p-value?
  
# g) The p-value should be very close to 0 (well below 0.05). Assuming that this 
#    is the case, state the conclusion regarding whether the proportion of 
#    Americans willing to get a coronavirus vaccine has changed between May and 
#    September

# 6) For question (2), the null hypothesis is that the "psychic" does not do
#    better or worse than random chance; the alternative hypothesis is that
#    the psychic's predictive abilities are different than what is expected
#    from random chance. Complete the following to specify what would it mean
#    if a Type I or Type II error occured.

#    A Type I error means that we conclude that ___________________________
#    but in reality, _______________________________.

#    A Type II error means that we conclude that ___________________________
#    but in reality, _______________________________.

# 7) For question (5), the null hypothesis is that there is no change in 
#    the proportion of American adults willing to get a coronavirus vaccine 
#    between May and September; the alternative hypothesis is that the 
#    proportion has changed. Complete the following to specify what would it mean
#    if a Type I or Type II error occured.

#    A Type I error means that we conclude that ___________________________
#    but in reality, _______________________________.

#    A Type II error means that we conclude that ___________________________
#    but in reality, _______________________________.


