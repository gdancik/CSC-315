###############################################################
# k-nearest neighbor (KNN) example
###############################################################

###############################################################
# Let's get the GSE1297 data again (modified from limma.R)
###############################################################
library(GEOquery)
GSE1297 = getGEO("GSE1297")

###########################################################
# Pull out gene expression data and pheno type data
###########################################################
# For this dataset, expression data must be logged -- see limma.R
GSE1297.expr = exprs(GSE1297[[1]])
GSE1297.expr = log2(GSE1297.expr)
GSE1297.p = pData(GSE1297[[1]])

################################################################
# Finding differentially expressed genes between males and
# females
################################################################

library(limma)
gender = as.character(GSE1297.p$characteristics_ch1.1)
levels(gender) = c("Female", "Male")
design = model.matrix(~0+gender)
colnames(design) = c("Female", "Male")

## limma package fits a linear model to each row of the expression matrix ##
fit = lmFit(GSE1297.expr, design)

## Contrasts need to match column names of design matrix ##
contrast.matrix <- makeContrasts(Female - Male,levels=design)

## fit model based on contrasts (e.g., Female - Male)
fit = contrasts.fit(fit, contrast.matrix)

# calculate moderated t-statistics by moderating standard errors
# toward a common value, which makes answers more robust
fit = eBayes(fit)

## get top probes, sorted by p-value, from all probes (note: if number is not set, then
## only 10 probes at most are returned)
tt.05 = topTable(fit,sort.by = "p", p.value = 0.05, number = nrow(GSE1297.expr))


##############################################################################
# In classification problems, it is often desirable to scale each probe,
# so that no single probe has a dominating effect. Below is an example
##############################################################################

#####################################################
# prints Euclidian distance matrix and plots cluster
#####################################################
plot.clust <-function(x) {
  d = dist(t(x))
  print(d)
  h = hclust(d)
  plot(h)
}

## generate data
r = c(1.1, 2, 3) #each row
M = rbind(p1=r,p2=r,p3=r,p4=r) # 4 rows
colnames(M) = c("A","B","C")

M

# A and B are the closest
plot.clust(M)


# if a probe with large expression is added, A and C are the 'closest', but this
# is completely determined by probe 5 (p5)
M = rbind(M, p5=c(10,15,10))
M
plot.clust(M)

# scale each row so that each probe has the same weight
row.scale <-function(x) {
  x.scale = t(scale(t(x)))
  return(x.scale)
}

# with row scaling, A and B are the most similar
# with microarray data, it is usually desirable to scale each row (probe). If not,
# then probes with higher expression will carry more weigh the classification
M.scale = row.scale(M) # scale each row
plot.clust(M.scale)


###############################################################
## back to microarray data 
###############################################################

## get probes that are differentially expressed
m = match(rownames(tt.05), rownames(GSE1297.expr))
X = GSE1297.expr[m,]

# visualization of probes - with knn, an unknown observation (?)
# is classified based on it's k-nearest neighbors. Note: only
# 2 dimensions are shown, though 12 probes are used for 
# classification
col.gender = as.integer(as.factor(gender))
col.gender = c("pink", "blue")[col.gender]

X.scale = row.scale(X)
plot(X.scale[1,], X.scale[2,], col = col.gender, pch = 19, main = "top 2 probes, scaled")
points(-.5, -.5, pch ="?") # add unknown probe to be classified




###################################################################
# Let's predict the gender of each individual, using the 12 probes
# with FDR < 5%
# Note that knn requires samples in rows and features in columns
###################################################################

library(class) # required for knn

################################################
# Make a prediction for a set of observations
###############################################
X.test <- matrix(c(-4, -1, 0, -1, 0, 0, -1, -1, -1, 5, 0, 1, 
              3, -3, 0, 0, 0, -2, 2, 1, -2, -2, -2, 1,
              -2, 0, 1, 0, -2, -1, -2, -2, -3, 3, -1, 3), 
              nrow = 12)

# testing data should have same scaling as training data
X.test = row.scale(X.test)

# knn(train, test, classes, k), where 
#   train - training data (probes in columns and samples in rows)
#   test - testing data (probes in columns and samples in rows)
#   classes - the known classes corresponding to the training data
#   k - value of k for the 'k' nearest neighbors
preds = knn(t(X.scale), t(X.test), gender, k = 3)
preds

###################################################################
# Leave one out cross-validation (loocv):
# for each sample, predict its class after removing it from the 
# training set. Use knn.cv function
###################################################################
preds = knn.cv(t(X.scale), gender, k = 3) 
table(predicted = preds, true = gender)

# overall accuracy (% correct) #
sum(preds == gender) / length(gender)

#######################################################################
# The overall accuracy is not a good measure of performance because it
# depends on how balanced the data is.
# Sensitivity(A) is the probability of correctly classifying samples
# from group A. We will use average sensitivity as a measure of
# accuracy. Note that this function assumes that there are two
# groups
#######################################################################

avg.sensitivity <-function(predicted, true) {
   t = table(predicted = predicted, true = true)
   if (nrow(t) != 2) {
     stop("only 1 row in accuracy table")
   }
   acc = diag(t) / colSums(t)
   sens1 = acc[1]
   sens2 = acc[2]
   list(table = t, sensitivity1 = sens1, sensitivity2 = sens2, avg.sensitivity = mean(acc))
}

avg.sens = avg.sensitivity(preds, gender)
avg.sens

############################################################
# Manually make leave one out prediction for the 1st sample
###########################################################

# remove test sample from training data 
i.test = 1
X.train = t(X.scale[,-i.test])   
Y.train = gender[-i.test] 

# isolate test sample so we can classify it
X.test = t(X.scale[,i.test, drop = FALSE])  # drop maintains the dimensions of the array
Y.test = gender[i.test]

## Use the 'knn' function to make predictions (Here we use k = 3)
p.test = knn(X.train, X.test, Y.train, k = 3)
p.test
Y.test


########################################################################
# General classification procedure:
#   Choose a classifier, and use leave-one-out cross-validation to 
#   estimate the number of probes/genes to use (based on number, or FDR), 
#   and classification parameters such as 'k' in knn. Then evaluate
#   the classfier on an independent dataset
########################################################################