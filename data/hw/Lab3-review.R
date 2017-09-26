# TO DO: set library path if necessary
library(ggplot2)


# 1) Construct the contingency table in R, along with the table of 
#    conditional proportions. Does there appear to be a relationship 
#    between Income and happiness? Why or why not?

happiness <- cbind(Not.Too.Happy = c(26,117,172), Pretty.Happy = c(233, 473,383), Very.Happy = c(164,293,132))
rownames(happiness) = c("Above Average", "Average", "Below Average")

# TO DO: need to construct table of conditional proportions and answer the 
# above question


# 8) Construct a scatterplot that predicts gas mileage from the vehicle’s 
#    weight, and add the corresponding regression line. Describe the 
#    relationship between weight and miles per gallon based on these results.

ggplot(mtcars, aes(wt, mpg)) + geom_point() + 
         geom_smooth(method = "lm", color = "darkred") +
         ggtitle("Relationship between car weight and gas mileage in 1974") +
         labs(x = "Weight in 1000 pounds", y = "miles per gallon (mpg)")
                            


# 9) Find the linear regression line that predicts miles per gallon from weight. 
#    Find and interpret the y-intercept. Find and interpret the slope.



