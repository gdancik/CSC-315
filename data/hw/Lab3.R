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

# 1) The code below creates a contingency table of the results from 
#    the clinical trial evaluating the effectiveness of the
#    Pfizer-BioNTech vaccine (actually named BNT162b2) for preventing
#    COVID-19 in individuals aged 16 years and older. Note: a 'placebo'
#    is a 'fake' version of the treatment, such as a salt solution or sugar
#    pill. The placebo is necessary because individuals may respond
#    positively if they believe they are receiving the appropriate 
#    treatment, even if they are not.
#
#    Data source: https://www.nejm.org/doi/full/10.1056/nejmoa2034577

vaccine_trial <- cbind(`No Infection` = c(21566, 21712),
                   `Covid-19 Infection` = c(162, 8))
rownames(vaccine_trial) <- c('Placebo', 'Vaccine')

#  (a) Create the appropriate conditional proportion table, 
#    where Vaccination Status is the explanatory variable. 
# 
#  (b) The conclusion of the clinical trial was that the vaccine was 95%
#      effective. Let's verify this. First, calculate the relative risk,
#      which has the formula 

#      RR = proportion of vaccinated individuals who are infected / 
#           proportion of unvaccinated individuals who are infected            
#       
#     Note: To carry out this calculation, you should directly access 
#     elements of your conditional proportion table (so that your 
#     calculation will be correct even if the data changes)
#
#     Then calculate vaccine efficiency, which is given by:

#     Vaccine efficiency = (1-RR) * 100%
#
#     Vaccine efficiency measures the reduction of cases expected from the 
#     vaccine (i.e., an efficiency of 95% means vaccination reduces the
#     number of cases by 95%) 

#    (c) Based on the above analysis, does there appear to be an association
#        between this vaccine and COVID-19 infection? Why or why not?


# For the next several questions, we will analyze our class survey data, 
# which is read in using the following code: 

library(readr)
survey <- read_csv('https://gdancik.github.io/CSC-315/data/datasets/csc315_survey_fall2022.csv')

# 2) Construct a table of conditional proportions for whether someone is a 
#    Cat or Dog person, based on their Favorite Season. Does there appear
#    to be a relationship between these variables? Why or why not?

# 3) Construct a stacked bar graph that shows whether there is an 
#    association between someone's preferred mode (light or dark),
#    and whether someone is a drinker. Does there appear to
#    be an association between these two variables in our class?
#    Why or why not? Note: the line below creates a factor variable for
#    whether each person is a drinker or not.

drinker <- factor(survey$Alcohol > 0, labels = c('non-drinker', 'drinker'))


# 4) Construct a scatterplot of HS GPA vs. College GPA, so that 
#    College GPA would be predicted from HS GPA, and add the 
#    regression line from the corresponding linear model. 

# 5) Calculate the correlation and describe the association between 
#    HS and College GPA, using the language given by the following site:
#    https://www.dummies.com/education/math/statistics/how-to-interpret-a-correlation-coefficient-r/
#    Based on these results, is HS GPA a good predictor of College GPA
#    in this class?

# For the remaining questions, you will use the 'mtcars' dataset,
# which contains data on 32 cars extracted from the 1974 Motor 
# Trend US magazine. This dataset is available in R in the 
# data.frame 'mtcars', and can be viewed using View(mtcars)
  
# The two variables we will examine are wt, the weight of the car in 
# thousands of pounds, and mpg, the gas mileage in miles per gallon
# from road tests. Additional information about the dataset can be 
# found by typing ?mtcars in the R console. 
  
# 6) Construct a scatterplot that predicts gas mileage from the 
# vehicle’s weight, and add the corresponding regression line. 
# Describe the relationship between weight and miles per gallon 
# based on this graph.

# 7) Find the linear regression line that predicts miles per gallon 
#    from weight. Find and interpret the y-intercept. Find and 
#    interpret the slope.

# 8) Based on this set of cars (in 1974), what would you predict 
#     the miles per gallon to be for a car that weighed 3000 pounds? 
#     What would you predict the miles per gallon to be for a car 
#     that weighed 7000 pounds? (Remember that if the prediction 
#     would be an extrapolation, you should say so and not make 
#     this prediction). 

# For the remaining questions, you will construct graphs
# of COVID-19 infections using data available from the
# state of Connecticut (https://data.ct.gov/Health-and-Human-Services/COVID-19-County-Level-Data/ujiq-dy22)
# This data is up-to-date, but only contains records going back to June 2022.

library(dplyr)
covid <- read_csv('https://data.ct.gov/resource/ujiq-dy22.csv')

# 9) Create a scatterplot of the total number of COVID cases
#    over time for Windham county only. You will plot 'report_date'
#    on the x-axis and 'cases_7days' on the y-axis (which is the
#    number of cases reported over the past 7 days). Add a geom_smooth
#    layer with default parameters (you should get a curve, not a line).

# 10) Create a scatterplot of the total number of COVID cases
#    over time for all counties, after filtering out results where
#    the county is "not_available". Your code will be similar to (9) but
#    you should color code the points by county by setting the 'color' 
#    aesthetic to the county column. Also add a 'geom_smooth' layer 
#    using the default method. You should also do the following:
#     - Use 'case_rate_weekly' on the y-axis, which gives you the 
#       number of cases per 100,000 individuals, so that we can
#       compare case numbers across counties. 
#     - in the geom_smooth layer, set se = FALSE, which removes the
#       gray bands which reflect variability in the curve
#     - in the geom_point layer, set size to 0.5 to decrease the size
#       of the points

#    What two counties have been hit hardest by COVID, 
#    in terms of number of cases, based on this data?