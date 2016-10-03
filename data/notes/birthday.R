############################################################
## the birthday problem
############################################################

## anyDuplicated function 
x = 1:10
anyDuplicated(x) # returns 0 (FALSE) if no elements are duplicated
x[1:2] = 1
# returns index of first duplicated element (TRUE) if any are duplicated
anyDuplicated(x) 

############################################################
## randomly select 2 birthdays, assuming 365 days / year
############################################################
s = sample(1:365, 2, replace = TRUE)

############################################################
## function to randomly select 'n' birthdays, 
## assuming 365 days / year
############################################################
birthdays <- function(n) {
  s = sample(1:365, n, replace = TRUE)  
  return(s)
}


#######################################################
## randomly generates 'n' birthdays
## returns TRUE if any are the same; FALSE otherwise
#######################################################
duplicate.birthdays <- function(n) {
  b = birthdays(n)
  dups = anyDuplicated(b)
  dups = any(dups>0)
  return(dups)
}

##################################################################
## returns the empirical probability (based on 5000 trials)
## that >1 person has the same birthday in a room with 'n' people
##################################################################
prob.duplicate <-function(n) {
    r = replicate(5000, duplicate.birthdays(n))
    prob = sum(r) / length(r)  
    return(prob)
}
  

##################################################################
## plot empirical probabilities for n = 1,2, ... 100 people
##################################################################
n = 1:100
probs = sapply(n, prob.duplicate)
plot(n, probs, type = "l", main = "The Birthday Problem", xlab = "# of people in room",
     ylab = "Probability at least 2 have the same birthday")
abline(h = 0.5, col = "red")
abline(v = 23, col = "red")
