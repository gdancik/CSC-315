########################################################
# Name:
# CSC-315
# Lab #2: Graphical and Numerical Summaries of Data
########################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions.
# Your Notebook must include output showing the requested tables and
# graphs, and answers to questions should be provided in comments.
# Also, all graphs must be given an appropriate title, x-axis label, and
# y-axis label. The ggplot2 library must be used to generate all
# graphs unless stated otherwise. When you are finished, submit an
# HTML Notebook through Blackboard as described previously.
#
# Note: DO NOT delete / modify any of the questions / comments below!
##########################################################################

# 1) Load our classes survey data (available at:
#  https://gdancik.github.io/CSC-315/data/datasets/csc315_survey_fall_2021.csv)
#   and add the code for this to the script. Note: I suggest to 
# change the name of the survey data frame (such as to 'survey') so
# that it is easier to type!


# 2) How many students completed the survey (i.e., how many rows are there)?

# 3) How many questions were asked (i.e., how many columns are there)?

# 4) Display the 'CatOrDogPerson' value for the 4th person (and only this person)

# 5) (a) Construct a frequency table for a person's favorite meal.

#    (b) Construct a relative frequency table for a person's favorite meal. 

#    (c) Based on this data, what meal does this class like the most?

# 6) Construct a frequency bar graph for the "Mobile" data, which contains the 
#    response to the question "What type of mobile phone do you prefer?".
#    The bars should be colored according to mobile phone preference, using the default 
#    colors, and thegraph should not have a legend. What do you conclude
#    about mobile phone preference in this class?


# 7) Construct a Pareto Chart of a person's favorite season
#    (you may display either the frequency or relative frequency). What
#    do you conclude about the choice of favorite season in this class?


# The code below generates a relative frequency table for whether
# a student consumes alcohol (i.e., alcohol > 0).

# Note: this assumes the survey data frame is stored in 'survey';
# you must edit the next line if this is not the case
consumes_alcohol <- table(survey$Alcohol > 0)

# Now FALSE means that Alcohol is not > 0; TRUE means Alcohol > 0
names(consumes_alcohol)  

# re-name values
names(consumes_alcohol) <- c('No', 'Yes')  

# convert to relative frequencies and a data frame
consumes_alcohol <- prop.table(consumes_alcohol)
df_alcohol <- data.frame(consumes_alcohol)
df_alcohol

# 8) Using the df_alcohol data frame, generate a relative frequency bar graph 
#    describing alcohol consumption in this class. If a student is selected at random
#    are you more or less likely to select a student who drinks?

# 9) Do Android users prefer fruits or vegetables? Answer this
#    question by using dplyr's 'filter' function to create a new data frame for
#    "Android" people, then generate a relative frequency table for the 
#    "FruitsOrVeggies" column.


# 10) Construct a histogram for alcohol consumption using the 'hist' 
#     function with default parameters. Note: you should be looking
#     at the original Alcohol data, not the Yes/No values from 
#     question (8). Describe the shape of its distribution. 
#     In particular, describe whether it is unimodal, bimodal, or flat, and 
#     whether it is skewed right, skewed left, or symmetric?


# 11) Calculate the mean and median for Alcohol consumption. 
#     Which is a better measure of averages (based on the shape of the 
#     distribution)? 

# 12) Create side-by-side boxplots of College GPA
#     based on whether a person is a "cat" or "dog" person, 
#     and answer the questions below:
#     (a) Does there appear to be a different in College GPA
#         between "cat" and "dog" people in this class?
#     (b) Are there any outliers? If so, how many, and for which group?



# 13) (a) Find the variance and standard deviation for College GPA using
#         the R functions 'var' and 'sd'. 
#
#     (b) Use the 'sqrt' function to confirm that the standard deviation
#         is the square root of the variance, by writing a logical statement
#         that evaluates to TRUE.
