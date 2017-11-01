library(gtools)
library(reshape2)
library(ggplot2)

# CSC 315, Exam II 

# Name: _____________________

# Directions: Answer the questions below.

# 1. Suppose that the majors of 20 students in a particular class is given by
#    the following frequency table:

class <- data.frame(Major = c("Computer Science", "Mathematics", "Biology"),
                    Frequency = c(13, 5, 2))
class

# (a) If a student is selected at random, what is the probability that
#     the student is a computer science major? 



# (b) If a student is selected at random, what is the probability that the 
#     student is NOT a Mathematics major?



# 2. (a) Find the empirical probability that a person rolls two dice with 
#    a sum of 11, by considering 5000 simulations where the two dice are
#    rolled. Answer this question by using the function below, which
#    simulates two die rolls, and returns TRUE if the sum of the two 
#    rolls is equal to 11; and returns FALSE otherwise 

roll.11 <- function() {
  rolls <- sample(1:6, 2, replace = TRUE)
  return (sum(rolls) == 11)
}

# (b). You will now use classical probability to find the same probability, 
#      by using the sample space and the sum of each outcome below.

sample.space <- permutations(6, 2, 1:6, repeats.allowed = TRUE)

# show the first 8 possible outcomes
sample.space[1:8,]

# find the sum of each possible outcome
r <- rowSums(sample.space)

# display the first 8 sums
r[1:8]








# 3. Normal probability questions 
#    (Write the R code that will output the correct probability)

# (a) Suppose that X ~ N(21,2). Find P(X > 25). 



# (b) In a population that is normally distributed, calculate the probability 
#     that a randomly selected observation is more than 1.4 standard 
#     deviations from the mean (in either direction)



# (c) Suppose that X ~ N(86,2), and answer the following questions.
#    
#      (i) In a sample of size 20, what is the expected value of the
#         sample mean?



#      (ii) In a sample of size 20, what is the standard deviation of the 
#           sample mean?



#      (iii) Find the probability that the sample mean is at least 87.




# 4. In 2013, the proportion of adults who smoke in the U.S. was 0.18. 
#   A 2015 study involving 1000 adults found that 163 of them smoked.
#   Is there evidence that the smoking rate has changed? The hypothesis 
#   test is carried out using the code below.

results <- prop.test(163, 1000, p = 0.18)
results

# (a) State the null and alternative hypotheses, making sure to define the
#     relevant parameters




# (b) Extract and calculate the Z test statistic from the results object, and 
#     assign this to the variable Z.


# (c) Extract the p-value from the results object



# (d) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.






# 5. Consider a null hypothesis about a population proportion or comparing two 
#   proportion, and the following Z test statistics. Write the R code that would
#   calculate the p-value for the hypothesis test

# (a) Z = -1.76


# (b) Z = 2.9


# 6. Consider a null hypothesis about a population mean 
#   and the following t test statistics and sample sizes. Write the R code that would
#   calculate the p-value for the hypothesis test.

# (a) t = -3.11, n = 23


# (b) t = .29, n = 97


# For questions (7) - (8), state the null and alternative hypotheses 
# (these must include parameters such as p and mu, and their descriptions)


# 7. A study is conducted to determine whether or not the proportion of females in
#   the United States is 50%.






# 8. A study is conducted to determine whether or not the average age of an 
#    adult male in the U.S. is different than the average age of an adult female.






# 9. Consider the hypotheses given by the following, to determine if a student
#    has a preference for Coca-Cola or Pepsi
        
#       H0: p = 0.50
#       HA: p != 0.50, where p is the probability that a person prefers Coca-Cola 
#           over Pepsi
#                    
#   (i) What would it mean in the context of this problem if a Type I error 
#       occured?





#   (ii) What would it mean in the context of this problem if a Type II error 
#        occured?
    




# 10. In Lab 6, we saw that there was sufficient evidence that Americans thought 
#     that  “strengthening law and order” should be a bigger priority than 
#     “reducing bias against minorities” for the U.S. criminal justice system, 
#     based on an October 2016 survey. In this question, we will look at whether 
#     or not this preference was associated with political affiliation 
#     (Democrats vs. Republicans). The data and conditional probabilities are
#     provided below. Complete the steps below to determine whether is an 
#     association between preference and political affiliation.

counts <- cbind(LawAndOrder = c(300, 197), ReducingBias = c(66, 370))
rownames(counts) <- c("Republicans", "Democrats")
counts

# total number of Republicans and Democrats
rowSums(counts)

conditional <- prop.table(counts, 1)
conditional

# Complete the steps below regarding the hypothesis test given below:

# H0: p_democrat - p_republican = 0
# HA: p_democrat - p_republican != 0, where p_democrat is the proportion of 
#     democrats choosing "Law and Order" as a priority, and p_repbulican 
#     the proportion of republicans

# a) Write the appropriate R code that carries out the hypothesis test 
#    above, and assigns the result to the object 'res'



# b) The p-value for this hypothesis test is < 2.2e-16, which is < 0.05. Based 
#    on the p-value, what is your conclusion regarding the null and alternative 
#    hypothesis, in the context of this problem?






# Extra Credit

# We have discussed that for a two-sample t-test, the degrees of freedom is 
# calculated using the Satterthwaite formula, but can be obtained by extracting
# the relevant object after using the t.test function is called. Write a function 
# called 'df' that returns the degrees of freedom from a two-sample t-test when 
# the vectors of observations are stored in 'x' and 'y'.

