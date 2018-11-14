library(limma)
library(class)

rm(list = ls()) ## this will delete all objects in the workspace

###################################################################
# Useful functions
###################################################################

## scale rows of given matrix to have mean 0 and sd of 1
row.scale <-function(x) {
  x.scale = t(scale(t(x)))
  return(x.scale)
}

# returns performance information, incuding average sensitivity
avg.sensitivity <-function(predicted, true) {
  predicted = factor(predicted)
  true = factor(true)
  levels(predicted) = levels(true)
  t = table(predicted = predicted, true = true)
  if (nrow(t) != 2) {
    stop("only 1 row in accuracy table")
  }
  acc = diag(t) / colSums(t)
  sens1 = acc[1]
  sens2 = acc[2]
  list(table = t, sensitivity1 = sens1, sensitivity2 = sens2, avg.sensitivity = mean(acc))
}


###################################################################
# Load the Challenge Data
###################################################################

# X.train - the log2 gene expression data for the training samples
# Y.train - the class labels (NMI or MI) for the training samples
# X.test - the log2 gene expression data for the test samples

load(url("https://gdancik.github.io/CSC-315/data/hw/Challenge.RData"))


# Note: to get your predictions in text format, with one prediction per line, 
# and with no other text, use the write.table function with row.names = FALSE


