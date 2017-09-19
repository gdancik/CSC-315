########################################################
# Name: In-class Review (partial)
# CSC-315
# Lab #2: Graphical and Numerical Summaries of Data
########################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions

# Note: all graphs must be given an appropriate title, x-axis label, and
#   y-axis label
##########################################################################

# 1.load our classes survey data 
#   (available at https://gdancik.github.io/CSC-315/data/datasets/survey_fall-2017.csv)
#   and add the code for this to the script. 

survey <- read_csv("https://gdancik.github.io/CSC-315/data/datasets/survey_fall-2017.csv")


# 4. Construct a frequency bar graph for the response to "Are you a cat or a dog person?"
#    Remove the legend by adding the following component to the end of your 
#    ggplot() code: theme(legend.position = "none")


# 8. Out of the "Cat" people in this class, what proportion list "Python" as their 
#    favorite programming language? Answer this question by first creating a new
#    data.frame for "Cat" people only. Then generate a relative frequency table
#    for favorite programming language. Then answer the same question for "Dog" 
#    people. What do you conclude about programming language preference (for 
#    students in this class) based on this data?


# 13. Create side-by-side boxplots showing the College GPA based on
#     a person's favorite programming language.  Does there appear 
#     to be a significant difference in the GPAs between these 
#     groups? Are there any outliers? If so, how many?


