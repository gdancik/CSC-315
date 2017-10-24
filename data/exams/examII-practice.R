# CSC 315, Exam II Practice Problems

# Note: This is not a comprehensive review, but contains exercises covering some
# of concepts that may appear on Exam II. In addition to these exercises,
# make sure you understand concepts covered in lecture and on the previous labs.

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. 


# 1. Blackjack time! A player has a blackjack if dealt two cards,
#    and one is an ace (a 1) and the other is a 10-13 (10,J,Q,K)
#    What is the probability that a player is dealt a blackjack?

# (a). Answer this question by first writing a function that 
#    randomly draws a 2-card hand from a deck, and returns TRUE
#    if the hand is a blackjack. Then simulate 5000 blackjack games
#    and find the empirical probability that the player is 
#    dealt a blackjack. Note: a full deck is given by the following
#    code:

deck <- rep(1:13,4)

# (b). You will now use classical probability to answer this queestion. The code
#     below generates the sample space for all 2-card hands
  
library(gtools)
hands <- combinations(52, 2, deck, repeats.allowed = FALSE, set = FALSE)

#   (i) How many possible 2-card hands are there?
#   (ii) Find the probability that a player is dealt a blackjack
#   (ii) Find the probabiliy that a person has a hand totalling
#        16 or fewer points (in which case the dealer has to hit). 
#       For this problem, we'll assume that an Ace is worth 1 point.
#       10,J,Q, and K (10-13) are worth 10 points    


#2. In 2013, the proportion of adults who smoke in the U.S. was 0.18. 
#   A 2015 study involving 1000 adults found that 163 of them smoked.
#   Is there evidence that the smoking rate has changed?

# (a) State the null and alternative hypotheses
# (b) Calculate / find the test statistic (and specify the degrees of freedom)
# (c) Find the p-value 
# (d) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.
# (e) What would it mean in the context of this problem if a Type I error occured?
# (f) What would it mean in the context of this problem if a Type II error occured?


#3. Consider a null hypothesis about a population proportion or comparing two 
#   proportion, and the following Z test statistics. Find the p-value, and
#   state whether you would REJECT the null hypothesis, or FAIL TO REJECT
#   the null hypothesis.

# (a) Z = -1.76
# (b) Z = 2.9
# (c) Z = 0 (additional question: why would Z be equal to 0)? 
# (d) Z = 2.2
# (e) Z = -0.19


#4. Consider a null hypothesis about a population mean 
#   and the following t test statistics and sample sizes. Find the p-value, and
#   state whether you would REJECT the null hypothesis, or FAIL TO REJECT
#   the null hypothesis.

# (a) t = -3.11, n = 23
# (b) t = .29, n = 97
# (c) t = -1.3, n = 348
# (d) t = 5.3, n = 45
# (e) t = 0.23, n = 62


# For questions (5) - (7), state the following:
# (a) The null and alternative hypotheses (these must include parameters 
#     such as p and mu)
# (b) What does it mean in the context of this problem if a Type I error occured?
# (c) What does it mean in the context of this problem if a Type II error occured?

#5. A study is conducted to determine whether or not the proportion of females in
#   the United States is 50%.

#6. A study is conducted to determine whether or not the average age of an adult
#   male in the U.S. is different than the average age of an adult female.

#7. A study is conducted to determine whether a new drug reduces flu symptoms in
#   a group of patients. The new drug is compared to an old one.

#8. Suppose that for the new drug, 43 / 100 individuals have reduced flu symptoms, and for those
#   on the old drug, 37 / 111 individuals have reduced flu symptoms. Find the test statistic
#   and p-value for whether or not one drug is more effective than the other, and state your
#   conclusion based on the p-value.
    
