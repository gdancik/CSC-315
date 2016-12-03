# CSC 315, Exam III

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
# "versicolor" (type "virginica will be ignored)
#
# (a) Find and interpret the y-intercept of the regression line in the
#      context of this problem.
# (b) Find and interpret the slope in the context of this problem
# (c) What is the p-value for the hypothesis test that there is a
#     significant difference in Petal Width between the two flower types?

#  The code below constructs a scatterplot of this data (the iris data is available
#  by default)

col = as.integer(iris$Species)
plot(iris$Petal.Length, iris$Petal.Width, pch = 19, col = col)
legend("topleft", legend = levels(iris$Species), col = 1:3, pch = 19)


# 2.The dataset GSE29561 on the Gene Expression Omnibus (GEO) contains 
#   gene expression profiles of breast cancer tumors treated by chemotherapy
#   (epirubicin and docetaxel). Read in the processed data for GSE29561, 
#   and extract the expression data and phenotype data. Note: this data
#   is already on the log2 scale. (a) How many samples are there? 
#   (b) How many probes are profiled?


# 3. In the phenotype data table, characteristics_ch1.1 indicates whether 
#    the sample responded to the therapy (treatment response = 1) or 
#    not (treatment response = 0).  What proportion of samples responded 
#    to the therapy?


# 4.  Is the probe 202271_at differentially expressed between samples 
#     that respond to therapy and samples that do not? 
#     Construct a boxplot of the expression of this probe across the two 
#     groups, and report the fold-change and p-value. 

# 5. Download the appropriate platform data and determine the gene that is 
#    associated with this probe (Note: the correct gene is FBXO28)? Using this
#    platform data, find how many probes in total there are for this gene.
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
#   6, 8, and all 10 probes above. A more consise version of the avg. sensitivity function
#   is provided below

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

