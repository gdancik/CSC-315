########################################################
# Name:
# CSC-315
# Lab #3: Associations
########################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions.
# Your Notebook should include output showing the requested tables and
# graphs, and answers to questions should be provided in comments.
# Also, all graphs must be given an appropriate title, x-axis label, and
# y-axis label. The ggplot2 library must be used to generate all
# graphs unless stated otherwise. When you are finished, submit an
# HTML Notebook through Blackboard as described previously.
#
# Note: DO NOT delete / modify any of the questions / comments below!
##########################################################################

# 1) The code below creates a contingency table for Happiness as it
#    relates to income level. Create the appropriate conditional 
#    proportion table, where Income is the explanatory variable. Does
#    there appear to be a relationship between Income and happiness? 
#    Why or why not? 

happiness <- cbind(`Not Too Happy` = c(26, 117,172),
                   `Pretty Happy` = c(233, 473, 383),
                   `Very Happy` = c(164, 293, 132))
rownames(happiness) <- c('Above average', 'Average', 'Below average')

# 2) Import our class survey data, which is available at:
#    https://gdancik.github.io/CSC-315/data/datasets/CSC315_survey_Fall_2020.csv

# 3) Construct a stacked bar graph that shows whether there is an 
#    association between whether someone is a Cat or Dog person, and  
#    their preferred Modality for this class. Does there appear to
#    be an association between these two variables in our class?
#    Why or why not?

# 4) Construct a stacked bar graph that shows whether there is an
#    association between whether someone is a Cat or Dog person, and
#    their preferred Video Conferencing platform. Does there appear to
#    be an association between these two variables in our class?
#    Why or why not?

# 5) Generate a contingency table for the variables CatOrDogPerson and
#    preferred Video Conferencing platform, and calculate the 
#    corresponding conditional proportions, conditional on whether
#    someone is a Cat or Dog Person. In a comment, specify the
#    proportion of Cat people that prefer Zoom, and the proportion
#    of Dog people that prefer Zoom. What is wrong with cat people?
#    (just kidding, don't answer that)

# 6) Construct a scatterplot of HS GPA vs. College GPA, so that 
#    College GPA would be predicted from HS GPA, and add the 
#    regression line from the corresponding linear model. 

# 7) Calculate the correlation and describe the association between 
#    HS and College GPA, based on your answers to (6) and (7).


# For the remaining questions, you will use the 'mtcars' dataset,
# which contains data on 32 cars extracted from the 1974 Motor 
# Trend US magazine. This dataset is available in R in the 
# data.frame 'mtcars', and can be viewed using the code below:
  
# Note: this code should be commented out before compiling the
# Notebook, as it can cause errors on some systems

  View(mtcars)

# The two variables we will examine are wt, the weight of the car in 
# thousands of pounds, and mpg, the gas mileage in miles per gallon
# from road tests. Additional information about the dataset can be 
# found by typing ?mtcars in the R console. 
  
# 8) Construct a scatterplot that predicts gas mileage from the 
# vehicle’s weight, and add the corresponding regression line. 
# Describe the relationship between weight and miles per gallon 
# based on this graph.

# 9) Find the linear regression line that predicts miles per gallon 
#    from weight. Find and interpret the y-intercept. Find and 
#    interpret the slope.

# 10) Based on this set of cars (in 1974), what would you predict 
#     the miles per gallon to be for a car that weighed 3000 pounds? 
#     What would you predict the miles per gallon to be for a car 
#     that weighed 7000 pounds? (Remember that if the prediction 
#     would be an extrapolation, you should say so and not make 
#     this prediction). 

