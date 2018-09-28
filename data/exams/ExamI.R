# CSC 315, Exam I 
# Name: 

# Directions: Modify this script by adding R code to answer the questions and/or 
# complete the steps below. When you are finished, create a Notebook and e-mail
# your Notebook to dancikg@easternct.edu with the subject "CSC-315, Exam I".

# Note: all graphs must have an appropriate title and axis labels, and must
#       be generated using ggplot()


# Set your library path if necessary
library(ggplot2)
library(dplyr)

# 1. Create a vector that contains the following numbers: 
#     1-100, and the number 200

# 2. In the vector you created in (1), change the first three
#    values to 0

# The next set of question are based on data from
# the Apple's App Store (collected July 2017).
# Questions will focus on the following columns:
#     price - the price of the app
#     user_rating - the average user rating (0-5 scale)
#     prime_genre - the category of the app
#     cont_rating - the content rating, which is the 
#               recommended age level 
# The dataset is read into R using the statement below
AppleStore <- read.csv("http://bioinformatics.easternct.edu/data/AppleStore.csv")

# 3. How many apps have been rated (e.g., how many rows are 
#       in the table)?

# 4. Construct a relative frequency table for the category of the 
#    app (prime_genre). Which app category has the most apps?

# 5. Find the mean user_rating for "Games" apps (where prime_genre is 
#     equal to "Games")

# 6. Run the code below to generate a histogram for the app price.
#    Describe the shape of the histogram. 

hist(AppleStore$price, breaks = 30)

# 7. Find the mean and median app price. For this data, which
#    is a better measure of average, and why?

# 8. Is there an association between app category (prime_genre) and 
#    recommended age (cont_rating)? Answer this question by 
#    creating a stacked bar graph that shows the conditional 
#    proportions for content rating for each prime_genre. Then
#    describe what the graph shows you regarding whether or not
#    there is an association.


# The next set of questions will use a dataset containing nutritional
# information for menu items from McDonald's. We will focus on the following
# columns:
#     Calories - the total number of calories
#     Fat - the total amount of fat, in grams

mcdonalds <- read.csv("http://bioinformatics.easternct.edu/data/mcdonalds.csv")

# 9. Construct a scatterplot that predicts calories from fat content
#    and add the corresponding regression line.


# 10. Find the linear regression line that predicts Calories from Fat. 
#     Find and interpret the y- intercept. Find and interpret the slope.

    

# Extra Credit: 

# 1. Complete the function below so that it returns the names of the 
#    maximum x-values. For example, if x = c(1,3,1) 
#    and x.names = c("A","B","C"), then max.item(x,x.names) would 
#    return "B" (since "B" has the largest x value of 3)

max.item <- function(x, x.names) {
  
}

# test that function call returns "B"
x <- c(1,3,1)
x.names <- c("A", "B", "C")
max.item(x, x.names)

# 2. In a single statement, find the McDonald's items (Item column) that 
#    have the most Calories, Cholesterol, and Sugars (columns 4,6, and 10).
#    These column indices are stored in the variable below.
index <- c(4,6,10)

