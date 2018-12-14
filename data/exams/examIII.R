# CSC 315, Exam III

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. When you are finished, create an R
# Notebook and e-mail your Notebook to dancikg@easternct.edu with the
# subject CSC-315, Exam III

library(ggplot2)
library(limma)
library(class)

# 1. The statement below reads in the survey data for our class. Fit a 
# linear model that predicts Hours of Sleep based on whether a
# person is a "Cat" or "Dog" person. Code your explanatory variable 
# using x = 0 for a "Cat" person and x = +1 for a "Dog" person

survey <- read.csv("https://gdancik.github.io/CSC-315/data/datasets/csc-315_survey_cleaned.csv")

# (a) Find and interpret the y-intercept of the regression line in the
#      context of this problem.

# (b) Find and interpret the slope in the context of this problem

# (c) Based on your linear model, what is the p-value for the hypothesis 
#     test that there is a difference in the average amount of sleep between
#     "Cat" and "Dog" people? 

# (d) Verify the p-value from (c) is identical to the p-value obtained from
#     the two-sample t-test with equal variance


# 2.The code below loads in a data set that contains gene expression of human cells
#   that were exposed to anthrax, which is a lethal toxin. The expression data is
#   stored in 'GSE34407.expr' while the phenotype data is stored in 'GSE34407.p'
#   Note: this data is already on the log2 scale.

load(url("http://bioinformatics.easternct.edu/RData/GSE34407.RData"))

#   Answer the following questions based on the expression data:

#    (a) How many samples are there? 

#    (b) How many probes are profiled?

#    (c) What is the mean expression value of the 3rd probe?

#    (d) Find the probe that has the highest mean expression. Recall that 
#        'which.max' can be used to find the index of the maximum value of a 
#        vector



# 3.  The probe 206991_s_at corresponds to the gene CCR5, a gene that
#     is expressed on the surface of macrophages (an immune system cell).
#     Determine if this probe is differentially expressed between 
#     treated and control samples, by constructing a boxplot (using ggplot)
#     showing the expression of this probe across the two groups.
#     The title of the boxplot should contain the probe name, the 
#     fold-change, and the p-value from the two-sample t-test. 





# 4.  Using the limma package, find all probes differentially expressed
#     between treatment and controls, using a false discovery rate (FDR)
#     of 1%.

#     (a) Show that there are 269 probes differentially expressed with an FDR of 1%.

#     (b) How many of these probes are expected to be false positives?



# 5. The top 20 probes are given below. Based on these probes, generate a heatmap and 
#    color the columns with "darkgreen" for control cells and darkred for treated cells.

probes <-c("232181_at", "238520_at", "212686_at", "223741_s_at", "226820_at", 
          "225298_at", "205698_s_at", "212757_s_at", "235052_at", "209933_s_at", 
          "228230_at", "202341_s_at", "227450_at", "49306_at", "226272_at", 
          "221539_at", "219714_s_at", "202342_s_at", "221666_s_at", "223095_at")

# 6. The top 20 differentially expressed genes are given below. Using DAVID 
#    (https://david.ncifcrf.gov/home.jsp), identify the KEGG pathways associated
#    with these 20 genes.

genes <- c("PPARGC1B", "TRERF1", "PPM1H", "TTYH2", "ZNF362", "PNKD", "MAP2K6", 
           "CAMK2G", "ZNF792", "CD300A", "HELZ2", "TRIM2", "ERP27", "RASSF4", 
           "RCAN3", "EIF4EBP1", "CACNA2D3", "TRIM2", "PYCARD", "MARVELD1")


#7. The iris dataset contains sepal and petal measurements for various 
#   species of iris flowers. In the code below, X is a matrix that includes each 
#   measurement, with each row a different scaled measurement, and each 
#   column a different flower. Answer parts (a) - (d) below.

# extract and scale the feature data 
X <- t(scale(iris[,1:4]))

# extract and scale the class labels
Y <- iris[,5]


# (a) Use the avg.sensitivity function below to find the average sensitivity for predicting species, 
#     using a leave-one-out cross-validation, and using the k-nearest
#     neighbor classification method with k = 3. 


# returns the average sensitivity for a set of predictions and corresponding
# true values
avg.sensitivity <-function(predicted, true) {
  true <- factor(true)
  predicted <- factor(predicted, levels = levels(true))
  t <- table(predicted = predicted, true = true)
  acc <- diag(t) / colSums(t)
  mean(acc)
}

# (b) Find the optimal value of 'k', considering 'k' values of 1,3,5,7,9,11

# (d) The 'newX' data frame contains scaled sepal and petal measurements for an
#     iris flower where the species is not known. Note that this data frame 
#     contains measurements in each column, instead of in each row. Use a knn 
#     classifier with k = 3 to predict the species for this flower.

newX <- data.frame(Sepal.Length = 1.2, Sepal.Width = 0.3, Petal.Length=1.5, Petal.Width=1.1)

# (d) When talking about classifiation, we used the average sensitivity rather than the
#     overall accuracy. Why is average sensitivity a better measure than accuracy? 
#     Give an example where accuracy would be a misleading measure of performance.




# Extra Credit

# The next statement loads in a made-up platform dataset with gene names and
# their descriptions. You are interested in the (made up) gene with the name 'BIOINF'. 
# Note: You are only interested in genes named 'BIOINF', and not interested
# in other genes (such as 'BIOINFORMATICS' or 'BIOINF-1') that contain 'BIOINF'. 
# In addition, some probes correspond to multiple genes. For example,
# "CS /// BIOINF", corresponds to the gene 'BIOINF'.

# load in the platform data, which is stored in 'pl'
load(url("http://bioinformatics.easternct.edu/RData/platform.RData"))


# (a). How many "BIOINF" genes are on this platform?

# (b). Create a vector containing all numbers between 1 and 1301 that are 
#      divisible by either 3 or 5. 

# (c). Use the 22nd, 35th, 398th, 501st, and 600th elements of your vector 
#      from (b) as an index for the "BIOINF" genes you identified in (a), and output
#      the description of these 5 genes

