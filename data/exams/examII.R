# CSC 315, Exam II

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. 

# When you are finished, use Knit to create an HTML file (don't forget to specify 
# the library path using the .libPaths() function if necessary. Once you produce
# the HTML file, it should be submitted through Blackboard at the provided link:
# https://ct-ecsu.blackboard.com/webapps/login/


# 1. (Basic probability) Suppose a basket contains 2 green balls, 1 red
#    ball, and 15 black balls. If a ball is randomly selected, what is the 
#    probability that it is green?


# 2. The function below returns TRUE if the vector 'x' has 3 (and only 3)
#    identical values, and FALSE otherwise

three.equal <-function(x) {
  t = table(x)                   
  if (length(t) == 1 & t[1] == 3) { 
    return (TRUE)
  }
  return(FALSE)
}

#    Using the above function, find the classical probability that if a person
#    rolls three dice, then all values will be the same. The code below generates
#    the sample space for 3 die rolls.

library(gtools)
rolls = permutations(6, 3, 1:6, repeats.allowed = TRUE)


# 3. Suppose that X ~ N(73,2). 
#   (a) Find P(X < 74)
#   (b) Find P(X > 70)
#   (c) Find the probability that the sample mean > 74, in a sample of size 20
   

# 4. In a population that is normally distributed, calculate the probability 
#    that a randomly selected observation is more than 2.8 standard deviations 
#    from the mean.

# 5. Consider a null hypothesis about a population proportion, and the 
#    following Z test statistics. Find the p-value, and
#    state whether you would REJECT the null hypothesis, or FAIL TO REJECT
#    the null hypothesis.

# (a) Z = -1.76
# (b) Z = 2.9

# 6. Consider a null hypothesis about a population mean and the following
#    t test statistics and sample sizes. Find the p-value, and
#    state whether you would REJECT the null hypothesis, or FAIL TO REJECT
#    the null hypothesis.

# (a) t = -3.11, n = 23
# (b) t = .29, n = 97

# 7. A study is conducted to determine whether or not the proportion of females in
#   the United States is 50%. The proportion of females in the sample is 0.54. 
#   A hypothesis test is carried out, for the following null and alternative 
#   hypotheses: 

#   H0: p = 0.50, where p = the proportion of adult females in the U.S. 
#   HA: p !- 0.50

# (a) If the p-value for the hypothesis test is 0.61, what is the conclusion
#     regarding the null and alternative hypotheses in the context of this problem?

# (b) If the p-value for the hypothesis test is 0.023, what is the conclusion
#     regarding the null and alternative hypotheses in the context of this problem

# (c) What would it mean in the context of this problem if a Type I error occured?

# (d) What would it mean in the context of this problem if a Type II error occured?


#8. A study is conducted to determine whether a new drug reduces flu symptoms in
#   a group of patients. The new drug is compared to an old one. State the
#   null and alternative hypotheses for this problem.


#9. In 2013, the proportion of adults who smoke in the U.S. was 0.18. 
#    A 2015 study involving 1000 adults found that 163 of them smoked.
#    Is there evidence that the smoking rate has changed?

#    The null and alternative hypotheses for this problem are:
#    H0: p = 0.18, where p = the proportion of adults who smoke
#    HA: p != 0.18

# (a) Calculate / find the test statistic 

# (b) Find the p-value 

# (c) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.

# 10. The code below reads in height data (in inches) for males (coded 0) 
#     and females (coded 1). Complete the steps below to test whether 
#     the heights of males and females differe, on average.

heights = read.csv("http://pastebin.com/raw/g7UdTFKG")

#    The null and alternative hypotheses for this problem are:
#    H0: mu_males - mu_females = 0
#    HA: mu_males - mu_females != 0

# (a) Define mu_males and mu_females as used in the hypotheses above

# (b) Calculate / find the test statistic 

# (c) Find the p-value 

# (d) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.




# Extra Credit. Sample size calculation.

# This problem looks at determining the appropriate number of coin tosses needed
# to determine whether or not a biased coin is biased (when the true probability
# of heads is 0.60)

# (a) Complete the function below which simulates flipping a biased coin 
# 'n' times. The function should then carry out a hypothesis test against 
# H0: p = 0.50, and return TRUE if the NULL hypothesis would be rejected (i.e.,
# if there is enough evidence to conclude that the coin is biased)

reject.it <-function(n) {
  
  # flip coin n times
  flips = sample(c("H","T"), n, replace = TRUE, prob = c(0.60, 0.40))
  
  # find the number of heads
  x = sum(flips=="H")
  
  # TO DO: carry out hypothesis test and return TRUE if H0 is rejected; 
  # otherwise return FALSE
  
}

# (b) Using the function above, as well as the replicate function, find 
# the empirical probability that you would reject the null hypothesis if 
# you flipped the biased coin 10 times, 100 times, and 1000 times. 
# Based on your empirical probabilities, if you were carrying out this study
# would it be better to flip the coin 100 times or 1000 times? Why?




