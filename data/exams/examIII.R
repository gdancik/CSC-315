# CSC 315, Exam III

# Name:

# Directions: Modify this script to add R code in order to answer the questions 
# and/or complete the steps below. 

# When you are finished, use Knit to create an HTML file (don't forget to specify 
# the library path using the .libPaths() function if necessary. Once you produce
# the HTML file, it should be submitted through Blackboard using the link provided


# The dataset GSE40791 on the Gene Expression Omnibus (GEO) contains gene expression
# profiles for lung tumor and normal samples. As was discussed on Piazza, the accompanying 
# article (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4612811/) is a great example
# of a real analysis that aligns with the work done in this course. For questions 1 - 13,
# you will use this dataset to reproduce some of the published findings. Note that for 
# various reasons not discussed, your results will not be identical to the results in 
# the article.

# load required packages; don't forget to set your library path if needed 
library(GEOquery)
library(limma)
library(class)

# the code below reads in the processed data and extracts the expression
# and phenotype data. Note that the data is alread on the log2 scale
GSE40791 = getGEO("GSE40791")
GSE40791.expr = exprs(GSE40791[[1]])
GSE40791.p = pData(GSE40791[[1]])

#1. How many samples were profiled?

#2. How many probes were profiled?

# 3. Generate a boxpot of expression values. In order to do this, color the 
#    boxes red, set the outline argument to FALSE (so the outliers are not 
#    plotted). In addition, give the plot a main title and an appropriate label
#    for the vertical axis. This boxplot will be similarto the boxplot on 
#    the bottom of Figure 1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4612811/figure/fig01/



# 4. In the phenotype data table, the column 'source_name_ch1' indicates whether 
#    the sample is from a normal lung or a tumor (lung adenocarcinoma) sample. 
#    Construct and output a frequency table for the number of normal and tumor
#    samples. Also construct and output a relative frequency table of the
#    corresponding proportions.


# 5.  Is the probe 205203_at differentially expressed between normal and 
#     tumor samples? Construct a boxplot of the expression of this probe 
#     across the two groups, and report the fold-change and p-value. Is
#     the probe differentially expressed? Why or why not?


# The code below uses the limma package to fit a linear model to 
# contrast expression between normal samples and tumor samples. 

# set up design matrix for tumor vs. nontumor
tumor = as.character(GSE40791.p$source_name_ch1)
design = model.matrix(~0+tumor)
colnames(design) = c("normal", "tumor")

# fit linear model, calculate contrasts, and apply eBayes step
fit = lmFit(GSE40791.expr, design)
contrast.matrix <- makeContrasts(tumor - normal,levels=design)
fit = contrasts.fit(fit, contrast.matrix)
fit = eBayes(fit)

# 6A. After running the code above, call the topTable function and then
# use the 'nrow' function to determine how many probes are differentially
# expressed with an FDR of 1%. 



# 6B. Output the top 3 probes, the logFC, and adjusted p-values only. 



# 7. Output the platform information, and download the appropriate platform 
#    data (which is GPL570). 


# 8. The top 3 differentially expressed probes are as follows:

probes = c("217046_s_at", "207547_s_at", "206209_s_at")

#    Using R and the platform data you just downloaded, find the gene 
#    symbols for these 3 probes. Note: these match the 1st, 5th, and 
#    2nd genes in Table 1 of the article -- 
#    (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4612811/table/tbl01/).


# 9. According to Table 1 , the top gene identified by the authors is AGER. 
#    Using the platform data you downloaded above, how many probes are there 
#    for the gene AGER (Note: The gene symbol information is in the 'Gene Symbol' 
#    column of the platform data? 


# 10. Based on the probes below, generate a heatmap and color the columns with
#     "darkblue" for "normal lung" and "darkred" for "lung adenocarcinoma".

probes = c("217046_s_at", "207547_s_at", "206208_at", "209469_at", 
            "209470_s_at", "215918_s_at", 
            "202200_s_at", "222416_at")


#11. Using a k-nearest neighbor (knn) classifier with k = 3, and the probes above,
#   find the average sensitivity based on a leave-one-out cross-validation
#   A simplified version of the avg. sensitivity function
#   is provided below:

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



#12. Optimize the value of k, by finding the average sensitivity for values
#   of k in 1,3,5, 7, ...., 19. Construct a scatterplot with values
#   of k on the x-axis and the average sensitivity on the y-axis, and
#   state which value of 'k' gives you the highest average sensitivity


#13. Let's now look at whether these same probes can be used to distinguish 
#    males from females (the gender information is in column 'characteristics_ch1.1' 
#    of the phenotype data). Select ONE of the options below:
#       (1) find the average sensitivity of a gender knn classifier (with k = 3),
#         based on a leave-one-out cross-validation; OR 
#       (2) generate a heatmap that is color-coded with pink for females and 
#           blue for males. 
#    Based on your results, what do you conclude about the ability of these
#    probes to distinguish between tumor vs. normal samples, and between males and 
#    females?


#14. A linear model that predicts height (in inches) by gender is fit, with
#  x = 1 corresponding to males and x = 0 corresponding to females. The slope
#  of this model is 4.32 and the y-intercept is 65.13. 
# (a) Interpret the y-intercept of the regression line in the context
#     of this problem
# (b) Interpret the slope of the regression line in the context of
#     this problem
# Note: There is no R code to show for this problem


# Extra Credit
# When analyzing gene expression data, there are several methods for working with
# genes that have more than one probe. In the lung cancer article, the 
# authors state that "[for] Probes matching more than one gene ... [the] 
# average [mean] value was used." 

# Write a function that takes a gene symbol name, the expression matrix,
# and the platform data, and returns the mean expression value of the gene 
# symbol for each sample. Use this function to find the average expression
# of all probes for the gene AGER


