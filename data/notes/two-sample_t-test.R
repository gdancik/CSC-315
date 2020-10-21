##########################################################################
# Example of two-sample t-test
##########################################################################

library(dplyr)

# The ChickWeight dataset includes data from an experiment that looked 
# at the effect of 4 different protein diets on chicken weight, where
# chickens were weighed every other day for 20 days, and at day 21.

# Here we will compare diet 1 and diet 2, and weights at day 21, which
# was the end of the experiment

df <- ChickWeight %>% filter(Time == 21, Diet %in% 1:2) %>%
  mutate(Diet = paste0('Diet_', Diet))

# Question: Is there evidence that diet 1 and diet 2 result in different
# average (mean) weights at day 21?

# (a) State the null and alternative hypotheses (done for you):


# (b) Create side-by-side boxplots (using ggplot) showing the 
#     distribution of weight for each diet. Make sure to label the 
#     y-axis and give the chart a title.

# (c) Carry out the two-sample t-test and report the p-value that
#     tests against H0.

# (e) Find the p-value 'manually' based on the test statistic and
#     appropriate degrees of freedom from the t.test result (which 
#     is stored in the $parameter object)

# (f) State the conclusion regarding the null and alternative hypotheses in 
#     the context of this problem.

