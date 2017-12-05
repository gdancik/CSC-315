# CSC 315, Exam III

library(ggplot2)
library(dplyr)
library(class)

# Name:

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. 

# When you are finished, use Knit to create an HTML file (don't forget to specify 
# the library path using the .libPaths() function if necessary. Once you produce
# the HTML file, it should be submitted through Blackboard using the link provided
# (NOTE: you will not submit these practice problems, but will for the actual Exam)


# 1.Fit a linear model that predicts Petal Width based on the type of flower from
# the iris dataset (see below). You should code your explanatory variable
# using x = 0 if the sample is of type "setosa" or x = +1 if the sample is of type
# "versicolor" (type "virginica" will be ignored)
#
# (a) Find and interpret the y-intercept of the regression line in the
#      context of this problem.
# (b) Find and interpret the slope in the context of this problem
# (c) What is the p-value for the hypothesis test that there is a
#     significant difference in Petal Width between the two flower types?

#  The code below constructs a scatterplot of this data
ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) + geom_point() +
  theme_linedraw() + ggtitle("Petal Length and Petal Width of Iris data set") +
  theme(legend.background =  element_rect(colour = "black"))


## this loads the expression (GSE29561.expr) and clinical (GSE29561.p) data

load(url("https://gdancik.github.io/CSC-315/data/exams/GSE29561.RData"))

# 2. (a) How many samples are there? 
#    (b) How many probes are profiled?


# 3. In the phenotype data table, characteristics_ch1.1 indicates whether 
#    the sample responded to the therapy (treatment response = 1) or 
#    whether the sample did not respond (treatment response = 0).  
#    What proportion of samples responded to the therapy and what
#    proportion did not?


# 4.  Is the probe 202271_at differentially expressed between samples 
#     that respond to therapy and samples that do not? 
#     Construct a boxplot of the expression of this probe across the two 
#     groups, and report the fold-change and p-value. 

# 5. The statement below loads the platform data for this data set.
#    a) Show that the gene associated with this probe is FBXO28? (Note:
#    the third character is a capital letter 'O' and not a zero)

load(url("https://gdancik.github.io/CSC-315/data/notes/GPL96.RData"))

#    b) find how many probes in total there are for this gene, using
#    the appropriate regular expressions as discussed previously
#    What are the probes for FBXO28?



# 6.	Using the limma package, find the top 5 differentially expressed probes
#     (sorted by p-value).  What is the overall FDR for these 5 probes? Is 
#     there evidence that these probes are differentially expressed? Why or 
#     why not? (Note: limma requires that the column names of your design 
#     matrix NOT be numbers)




# 7. Based on the probes below, generate a heatmap and color the columns with "darkgreen" for
#    "treatment resonse in vivo: 1" and "darkred" for "treatment response in vivo: 0". In 
#    other words, resistant individuals are darkred and sensitive individuals are darkgreen

probes = c("201690_s_at","212611_at","218721_s_at",
           "219463_at","205014_at","216092_s_at",
           "207375_s_at","202271_at","221551_x_at","219053_s_at")

#8. Using a k-nearest neighbor (knn) classifier with k = 3, classify samples with
#   "treatment response in vivo: 1" vs. "treatment response in vivo: 0", by finding the optimal
#   number of probes to use based on the average sensitivity in a leave-one-out cross-
#   validation. Answer this question by finding the average sensitivities for the first 2, 4
#   6, 8, and all 10 probes above. You can calculate the average sensitivity by using
#   the function below. This is followed by a statement which scales the expression data.

## returns average sensitivity (this is simplified from our original function)
avg.sensitivity <-function(predicted, true) {
  predicted = factor(predicted)
  true = factor(true)
  levels(predicted) = levels(true)
  t = table(predicted = predicted, true = true)
  if (nrow(t) != 2) {
    stop("only 1 row in accuracy table")
  }
  acc = diag(t) / colSums(t)
  mean(acc)
}

# scale the expression data by probe (each probe will have mean 0 and sd of 1)
X <- t(scale(t(GSE29561.expr)))

