########################################################
# Name:
# CSC-315
# Lab #2: Graphical and Numerical Summaries of Data
########################################################

##########################################################################
# Add R code to the script below and create a Notebook to complete
# the steps and explicitly answer the following questions.
# Your Notebook include output showing the requested tables and
# graphs, and answers to questions should be provided in comments.
# Also, all graphs must be given an appropriate title, x-axis label, and
# y-axis label. The ggplot2 library must be used to generate all
# graphs unless stated otherwise. When you are finished, submit an
# HTML Notebook through Blackboard as described previously.
#
# Note: DO NOT delete / modify any of the questions / comments below!
##########################################################################

# 1) Load our classes survey data (available at:
#   https://gdancik.github.io/CSC-315/data/datasets/CSC315_survey_Fall_2020.csv)
#   and add the code for this to the script. 


# 2) How many students completed the survey?

# 3) How many questions were asked (i.e., how many columns are there)?

# 4) Construct a frequency table for the response to whether someone is a
#    'Cat' or 'Dog' person.

# 5) Construct a frequency bar graph for the response to 
#    "What is your preferred modality for this class?",
#    where the bars are colored according to modality using the default 
#    colors. Remove the legend by adding the following component to the end of your 
#    ggplot() code: theme(legend.position = "none")


# 6) Construct a relative frequency table for preferred videoconferencing software. 
#     What proportion of students said that they preferred Microsoft Teams?


# 7) Construct a Pareto Chart using the frequencies for preferred videoconferencing
#    software (you may display either frequency or relative frequency). 


# 8) Construct a relative frequency table for whether or not a student consumes alcohol
#   (i.e., consumes alcohol > 0 days per week).
#    Do this by first creating a logical vector where TRUE corresponds to consuming
#    alcohol and FALSE corresponds to does not consume alcohol. Then create a relative
#    frequency of these TRUE and FALSE values. Tables and relative frequency tables
#    are stored as named vectors (e.g., x <- c(item1 = 1, item2 = 2)). Use the 'names' 
#    function to change the names from FALSE and TRUE to "Consumes alcohol" and
#    "Does not consume alcohol"


# 9) Out of the "Cat" people in this class, generate a relative frequency
#    table for their preferred videoconferencing platform. Answer this
#    question by using dplyr's 'filter' function to create a new data frame for
#    "Cat" people, then generate a relative frequency table for the "VideoConferencing"
#    column. Repeat the analysis to answer the same question for "Dog" people. 
#    What do you conclude about a person's choice of video conferencing platform?


# 10) Construct a histogram for Alcohol consumption, by using the hist() function with 
#     the argument breaks = 14 to set the number of groupings. Describe the shape of its 
#     distribution. In particular, describe whether it is unimodal, bimodal, or flat, and 
#     whether it is skewed right, skewed left, or symmetric?


# 11) Calculate the mean and median for Alcohol consumption. 
#     Which is a better measure of averages (based on the shape of the distribution)? 


# 12) What is the 75th percentile for the average number of hours of sleep a student gets??

# 13) Create side-by-side boxplots showing the average hours of sleep 
#     based on whether a person is a "cat" or "dog" person, 
#     and answer the questions below:
#     (a) Does there appear to be a different in the amount of sleep
#         between "cat" and "dog" people in this class?
#     (b) Are there any outliers? If so, how many, and for which group?



# 14) For college GPA, what is the variance and standard deviation?


# 15) Create a vector with 20 values that has a standard deviation of 0.

