# CSC 315, Exam I Practice Problems

# set your library path if necessary
library(dplyr)
library(ggplot2)
library(reshape2)
library(readr)

# Note: This is not a comprehensive review, but contains exercises
# covering some of the concepts that will be on the first exam. 
# In addition to these exercises, make sure you understand concepts 
# covered in lecture and on the previous labs.

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. 

# When you are finished, create a Notebook and submit this assignment
# through Blackboard (NOTE: you will not submit these practice problems, 
# but will for the actual Exam)


# 1. Create a vector called ages that contains the numbers 21, 24, 32, and 19

# 2. Create a vector called evens that contains all even numbers 
#     between 1 and 100, and the number 200.

# 3. Write a function in R called min.max which takes a vector x as an argument 
#    and returns a list containing two named elements, the minimum of x and the
#    maximum of x. Use this function to find the minimum and maximum of the 
#    vector ages from problem (1). 

# 4. Include the following code in your script to create a matrix filled with 
#    5 columns and 20 rows, that is filled with random numbers between 0 and 1.

      m = matrix(runif(100), ncol = 5, nrow = 20)

#   Find the median of each row and the median of each column.

# 5. The questions below are based on a dataset that is available here:
#    http://pastebin.com/raw/1csmBawE
#    Read this data into R and include the statement or statements for doing 
#    so below.      
      
      
# 6. Find the correlation between Alcohol Consumption and College GPA and 
#     describe the trend in the data 
      

# 7. Find the mean and standard deviation of College GPA for those who agree
#    with same sex marriage legalization, and those that disagree.
      
# 8. Construct side-by-side boxplots for FB usage based on whether or not
#    a person agrees or disagrees with same sex marriage 
#    Make sure to label the y-axis and give the chart a meaningful title. 
#    Is there an association between FB usage and views on same sex marriage?

# 9. Construct a histogram of alcohol consumption, and describe the shape of the histogram. 
#   Is it unimodal, bimodal, or flat. Is it skewed right, skewed left, or symmetric?

# 10.  Construct a bar graph showing the proportion of females who agree
#      with marijuana legalization and the proportion that disagree.
#      Give the bar graph a meaningful title and y-axis label. 
