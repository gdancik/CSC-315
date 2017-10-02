######################################
## Module 4: Probability
######################################

# Note: the 'gtools' library is needed to enumerate
# permutations and combinations

# make sure to change the .libPaths if necessary
# you can install 'gtools' by typing install.packages("gtools")
library(gtools)
library(ggplot2)

######################################
## simulate rolling a die 10000 times
######################################
roll.num <- sample(1:6, 10000, replace = TRUE)

## treat the results as a factor (qualitative variable) and construct a barplot
roll <- factor(roll.num)
ggplot() + geom_bar(aes(roll, fill = roll)) + 
  ggtitle("Outcome from rolling die 10000 times") +
  theme(legend.position = "none")


## function to find the proportion of sixes occuring in first 'i' elements of 'x' 
proportion.sixes <- function(i,x) {
  count <- sum(x[1:i]==6)
  return (count/i)
}

## apply this function to all integers from 1- 10000
props <- sapply(1:length(roll.num), proportion.sixes, x=roll.num)


## plot empirical and theoretical probabilities
six.plot <- ggplot() + geom_point(aes(x = 1:length(props), y=props), color = "blue") +
  ggtitle("Empirical probability of rolling a 6 with a fair die") +
  labs(x = "# rolls", y = "proportion of sixes rolled") +
  geom_hline(aes(yintercept=1/6, linetype = "theoretical probability"),color = "red") +
  theme_linedraw() +
  scale_linetype_manual(name = "", values = 2) # this adds the legend (values = 2 for dotted line)


# plot only first 100 rolls
six.plot + xlim(0,100)

# look at complete plot
six.plot                



## simulate flipping a fair coin 1000 times
coins = sample(c("H", "T"), 1000, replace=TRUE)

# generate bar graph of frequencies; the guides() function is
# used to suppress the legend for the 'fill' elements (the coins)
ggplot() + geom_bar(aes(coins, fill = coins)) + 
  ggtitle("Outcome of flipping a fair coin 1000 times") +
  labs(x = "outcome", y = "Frequency") +
  guides(fill = FALSE)
  

# calculate relative frequency table, then 
# generate bar graph of frequencies; the guides() function is
# used to suppress the legend for the 'fill' elements (the coins)
t <- table(coins)
p <- prop.table(t, margin = NULL) # with no margin, divide by total
df <- data.frame(p)

ggplot(df) + geom_bar(aes(x = coins, y = Freq, fill = coins), stat = "identity") + 
  ggtitle("Outcome of flipping a fair coin 1000 times") +
  labs(x = "outcome", y = "Frequency") +
  geom_hline(aes(yintercept=1/2, linetype = "theoretical probability"),color = "black") +
  scale_linetype_manual(name = "", values = 2) +
  guides(fill = FALSE) + ylim(0,1)


## simulate flipping a biased coin 1000 times 
##    (90% probability of Heads, 10% probability of tails)
coins <- sample(c("H", "T"), 1000, prob = c(.9,.1), replace=TRUE)

# calculate relative frequency table, then 
# generate bar graph of frequencies; the guides() function is
# used to suppress the legend for the 'fill' elements (the coins)

t <- table(coins)
p <- prop.table(t, margin = NULL) # with no margin, divide by total
df <- data.frame(p)

ggplot(df) + geom_bar(aes(x = coins, y = Freq, fill = coins), stat = "identity") + 
  ggtitle("Outcome of flipping a fair coin 1000 times") +
  labs(x = "outcome", y = "Frequency") +
  geom_hline(aes(yintercept=9/10, linetype = "theoretical probability"),color = "black") +
  scale_linetype_manual(name = "", values = 2) +
  guides(fill = FALSE) + ylim(0,1)



############################################################
## general method for repeating a probability experiment --
## write a function to do the experiment once, then
## use the replicate function to repeat the experiment
############################################################

## example function: flip a fair coin 2 times, 
#     determine if we get 2 heads
# returns TRUE if we get 2 heads, FALSE otherwise
flip.two.heads <- function() {
  f = sample(c("H", "T"), 2, replace = TRUE)
  count = sum(f=="H")
  return (count == 2)
}

# flip a coin 2 times (or flip 2 coins), repeat 1000 times
two.heads <- replicate(1000, flip.two.heads())

# find the empirical probability of getting 2 heads = 
# number of times we get 2 heads / number of experiments
prop.heads <- sum(two.heads) / length(two.heads)
prop.heads


#########################################################
## Question: flip a coin 3 times and determine the 
## probability that the coin lands on heads at least
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

x <- 1:3
any(x == 4)
any(x == 1)
all(x == 1)

x <- c(1,1,1)
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

S <- permutations(2,3, c("C", "I"), repeats = TRUE)

# A = answer all questions correctly
correct <- S == "C"
index <- apply(correct, 1, all)
S[index, ]

# A = answer at least 2 questions correctly
index <-  rowSums(  S == "C" ) >=2
S[index,]


##########################################################
# Sample space S for flipping a coin 3 times  #
##########################################################

S = permutations(2,3, c("H", "T"), repeats = TRUE)
# Probability of getting all Heads, P(H = 3)
all.heads <- apply(S== "H", 1, all)
sum(all.heads) / length(all.heads)

# Probability of getting at least 1 Head, P(H >= 1)
any.heads <- apply(S== "H", 1, any)
sum(any.heads) / length(any.heads)


# Note that the complement of at least one head (H >=1)
# is no heads (H=0) (or all tails)

# by the rule of complements, 
# P(at least one head) = 1 - P(no heads)
no.heads <- !apply(S== "H", 1, any)
p.no.heads <- sum(no.heads) / length(no.heads)
1 - p.no.heads


