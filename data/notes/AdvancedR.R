###################################
## writing functions in R
###################################

## example with two arguments ##
divide <-function(x,y) {
    return(x/y)
}
divide(3,4)

## example with a default argument ##
divide <-function(x,y = 1) {
  return(x/y)
}

### arguments are matched in order unless they are named
divide(3)         # 3 / 1 (since 1 is default y-value)
divide(3,2)       # 3/2 = 1.5
divide(2,3)       # 2/3 = .666
divide(y = 2, 3)  # 3/2 = 1.5

###################################################
# return is implicit if the function ends with
# an expression (but NOT an assignment)
##################################################

# same divide function as above, with implicit return
divide <-function(x,y = 1) {
  x/y
}

divide(3,2) ## 3/2 = 1.5

#################################################
## for loops - loop for each element in a 
##        vector or list
#################################################

# basic format: loop for each element in a vector (or list) 
xvalues <- c(1,3,2,9)
# 'traditional approach', loop over index
for (i in 1:length(xvalues)) {
    cat (xvalues[i], " ")
}

# loop over each value of x
for (x in xvalues) {
    cat (x, " ")
}

######################################################################
## Note: In most cases, functional approaches (apply, lapply, sapply)
##    are preferred over using 'for loops' in R since functional
##    approaches are simpler and more readable (once understood).
##    In R, a 'functional' is a function that takes another function
##    as an argument
#####################################################################

######################################################################
## apply - applies a function to a row (MARGIN = 1) 
## or column (MARGIN = 2) of a matrix. 
######################################################################

###########################################
# example: find max of each row of 'm'
###########################################
m <- matrix(1:20, ncol=5, byrow=TRUE)
rowMaxes <- apply(m, 1, max)
rowMaxes

###############################################################
## add the 'n' smallest numbers in a vector
## Note: could also return sum(sort(x)[1:n])
###############################################################
add.smallest <- function(x, n=2) {
  x.sorted <- sort(x)
  sum(x.sorted[1:n])
}
ans.row <- apply(m, 1, add.smallest) ## for each row
ans.row

ans.col <- apply(m, 2, add.smallest) ## for each col
ans.col

ans.row.4 <- apply(m, 1, add.smallest, n = 4)
ans.row.4

## alternative format using an inline (anonymous) function
ans.row <- apply(m, 1, function(x) {
    sum(sort(x)[1:2])
})

## another alternative format using an inline (anonymous) function
ans.row <- apply(m, 1, function(x) sum(sort(x)[1:2]))

## lapply applies a function to each object in a list, and returns a list
person <- list(name = "Bob", sibling.ages = c(43,21), pet.ages = c(8,3))
lapply(person, length)

## sapply does the same but returns a vector
sapply(person, length)

########################################
# if statements - follows C/C++/Java format
########################################
compare <-function(x, ref = 0) {
  if (x < ref) {
    cat("the number", x, "is less than", ref, "\n")
  } else if (x > ref) {
    cat("the number", x, "is greater than", ref, "\n")
  } else {
    cat("the number", x, "is equal to", ref, "\n")
  }
}

compare(5, 3)

#######################################################################
# Question set (use the grades matrix below):
#######################################################################

grades <- matrix(c(71,86,82,93,87,92,85,85,98,99,100,92),ncol=3,byrow=T)
rownames(grades) <- c("Steve", "Joe", "Jane", "Andrea")

# 1. Using the apply function, find the mean grade for each student, and
#    the mean grade for each assignment


# 2. The following function takes a vector of grades and returns the
#    corresponding letter grade, based on the average (mean). Use 
#    this function and apply to find the letter grade for each 
#    student
letterGrade  <- function(x) {
  m <- mean(x, na.rm = TRUE)
  if (m >=  90) {
    return('A')
  } else if (m >= 80) {
    return('B')
  } else if (m >= 70) {
    return('C')
  } else if (m >= 60) {
    return('D')
  }
  return('F')
}

# 3. Write a function called 'is.A' that takes a vector, and returns TRUE 
#    if the mean value of the vector is >= 90 (in the A range). Then use
#    this function and 'apply' to identify the names of
#    the students with an A average. Can you write code that outputs 
#    only the names of the students with As?

