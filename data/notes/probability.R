######################################
## Module 4: Probability
######################################

# Note: the 'gtools' library is needed to enumerate
# permutations and combinations
# make sure to change the .libPaths if necessary
# you can install 'gtools' by typing install.packages("gtools")

library(gtools)

######################################
## simulate rolling a die 10000 times
######################################
roll = sample(1:6, 10000, replace = TRUE)
barplot(table(roll), 
        main = "Outcome from rolling die 10000 times",
        ylab = "Frequency")

## proportion of 6s in first 'i' elements of 'x' ##
prob.6 <- function(i,x) {
  sum(x[1:i]==6) / i
}

## apply this function to all integers from 1- 10000
probs = sapply(1:10000, prob.6, x=roll)

## plot cumulative proportion of 6s rolled ##
plot(probs, xlab = "# of rolls", type = "l", ylim = c(0,1),
     ylab = "cumulative proportion", 
     main = "Empirical Probability of Rolling a six")
abline(h = 1/6, col = "red")


## simulate flipping a fair coin 1000 times
coins = sample(c("H", "T"), 1000, replace=TRUE)
t = table(coins)
barplot(t/sum(t), main = "Outcome of 1000 coin flips",
        ylab = "Proportion")
abline(h = 0.50, col = "red")

## simulate flipping a biased coin 1000 times
coins = sample(c("H", "T"), 1000, prob = c(.9,.1), replace=TRUE)
t = table(coins)
barplot(t/sum(t), main = "Outcome of 1000 coin flips",
        ylab = "Proportion", ylim = c(0,1))
abline(h = 0.90, col = "red")


#########################################################
## general method for repeating a probability experiment
## ex: flip a coin 2 times, count heads
#########################################################

flip.two.heads <- function() {
  f = sample(c("H", "T"), 2, replace = TRUE)
  count = sum(f=="H")
  return (count == 2)
}
two.heads = replicate(1000, flip.two.heads())
prop.heads = sum(two.heads) / length(two.heads)

#########################################################
## Question: flip a coin 3 times and determine the 
## probability that the coin lands on  heads at least
## 2 times
## Copy and modify the above code to find the empirical
## probability using the 'replicate' function
#########################################################



#########################################################
## classical probability definitions - when all outcomes
## are equally likely
#########################################################

##########################################################
## useful logical function: 
#   all(x) returns true if all elements in x are TRUE
#   any(x) returns true if any element in x is TRUE
##########################################################

x = rep(1,3)
any(x == 1)
all(x == 1)

x = 1:3
any(x == 1)
all(x == 1)



##########################################################
# We can use the 'permutation' function from the gtools
# library to do classical probability calculations
# the arguments for the permuation function include 
# (1) number of outcomes / trial, 
# (2) number of trials,
# (3) sample space for each trial
# (4) repeats = TRUE if outcomes can repeat across trials
# Permutations gives all possible arrangements where
#    order matters
##########################################################

S = permutations(2,3, c("C", "I"), repeats = TRUE)

# A = answer all questions correctly
correct = S == "C"
S[ apply(correct, 1, all),   ]

# A = answer at least 2 questions correctly
index =   rowSums(  S == "C" ) >=2
S[index,]


##########################################################
# Sample space S for flipping a coin 3 times  #
##########################################################

S = permutations(2,3, c("H", "T"), repeats = TRUE)
# Probability of getting all Heads, P(H = 3)
all.heads = apply(S== "H", 1, all)
sum(all.heads) / length(all.heads)

# Probability of getting at least 1 Head, P(H >= 1)
any.heads = apply(S== "H", 1, any)
sum(any.heads) / length(any.heads)


# Note that the complement of at least one head (H >=1)
# is no heads (H=0) (or all tails)

# by the rule of complements, 
# P(at least one head) = 1 - P(no heads)
no.heads = !apply(S== "H", 1, any)
p.no.heads = sum(no.heads) / length(no.heads)
1 - p.no.heads


