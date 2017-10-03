########################################################
# Name:
# CSC-315
# Lab #4: Probability 
########################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions
##########################################################################

# Basic Probability Questions -- Use R as a calculator and specify the
# answers below

# A standard deck of cards contains 52 cards, with 13 cards from each 
#   suit (hearts, clubs, diamonds, and spades).
 
#1. If one card is selected at random, what is the probability that it
#   is the ace of spades?


#2. What is the probability that it is NOT the ace of spades?


#3. What is the probability that it is an ace?


#4. What is the probability that it is an ace or a 4?


# Use R to answer the remaining questions. You MUST use R to 
# enumerate and analyze the sample space or to carry out 
# probability experiments (simulations) in order to find the 
# requested probability.


#5. If you roll two dice (each with the values 1-6), what is the 
#   probability that the sum of the numbers is 7? Answer this 
#   question by first using the 'permutations' function from 'gtools'
#   to enumerate the sample space obtained by rolling two dice.
#   Note: the correct sample space has 36 outcomes.


#6. Calculate the answer to (5) by finding the empirical probability.
#   In order to do this, first write a function that randomly rolls 
#   two dice and then calculates the sum. Then use the 'replicate' 
#   function to find the empirical probability based on 5000 experiments.



# A probability distribution of a discrete random variable gives
# the probability for each value of that variable. 


#7. Use the 'permutations' function to enumerate the sample space 
#   obtained from flipping a coin 3 times. Using this sample
#   space, find the total number of heads (we will call this X) in each 
#   set of 3 coin flips. Then create a relative frequency table (i.e.,
#   the probability distribution table for X = the number of heads in 3
#   coin tosses). Also create a bar graph of the relative frequencies,
#   making sure to appropriately label the x-axis, y-axis, and title 
#   of the graph. Note: your results should match the theoretical 
#   probabilities which are P(X=0) = 0.125, P(X=1) = 0.375, P(X=2) = 0.375,
#   and P(X = 3) = 0.125



#8. Find the empirical distribution of X = the number of heads in 3
#   coin tosses by first using the replicate function and the 'count.heads'
#   function below to create a vector of X values in 1000 experiments. 
#   Then create a relative frequency table (i.e., the empirical 
#   probability distribution table). You do not need to create a bar 
#   graph for this problem.


################################################################
# this function counts the number of heads in 'n' coin tosses
################################################################
count.heads <-function(n) {
  s = sample(c("H", "T"), n, replace = TRUE)
  num = sum(s=="H")
  return(num)
}



#   Poker Time! The commands below generate all possible poker hands.
#   Here we ignore the suit since this is not needed for the questions 
#   below.  We also use combinations instead of permutations to create
#   the sample space. With combinations the order (e.g., of cards in 
#   a hand) does not matter. The cards are also sampled without 
#   replacement (repeats.allowed = FALSE by default).
#   We also specify 'set = FALSE' to allow for duplicate values in the deck 
#   vector. Finally, each combination is equally likely, so classical 
#   probability can be used.

deck = rep(1:13,4)
hands = combinations(52, 5, deck, set = FALSE)

#9. How many possible poker hands are there?



#10. Show that the probability of being dealt a four-of-a-kind is
#   approximately 0.00024 (or 1/4165). Note: You MUST use the
#   apply function and the four.of.a.kind function below to find this. 
#   Because the hands matrix contains more than 2.5 million rows, 
#   this calculation will take time (~5 minutes). 
#   You should therefore test your code on a subset of the 
#   hands matrix first. There are two 4-of-a-kinds
#   if you look at the first 20,000 rows.

##############################################################
# this function returns true if a hand contains a 4-of-a-kind
##############################################################
four.of.a.kind <- function(x) {
  t = table(x)  # frequency table for cards in the hand
  m = max(t)    # how frequent is the most common card?
  if (m == 4) return (TRUE)
  return (FALSE)
}


#11. Create a function that determines whether a vector 'x' contains 
#   a full house (i.e., 3 of a kind and 1 pair). You can assume
#   that 'x' includes exactly 5 cards.


#12. Show that the probability of being dealt a full house is 
#    approximately 0.00144 (roughly 1/694). Note for testing purposes,
#    that there are 18 full houses in the first 20,000 rows of the
#    hands matrix



