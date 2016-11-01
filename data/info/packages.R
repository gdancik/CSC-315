############################################################
## CSC 315: Required packages 
## Note: make sure that your .libPaths() is set correctly
##    if working from a school computer
############################################################

## for knitting HTML documents ##
install.packages("knitr")
install.packages("yaml")
install.packages("htmltools")
install.packages("rmarkdown")

## for permutation function ##
install.packages("gtools")

## for colorRampPalette function used to color heatmaps ##
install.packages("RColorBrewer")

#################################################
# Bioconductor packages for Microarray Analysis
#################################################

## to install bioconductor packages, source the following once ##
source("http://bioconductor.org/biocLite.R")

biocLite("affy")
biocLite("affydata")
biocLite("leukemiasEset")
biocLite("GEOquery")
biocLite("limma")
biocLite("hgu133a.db")

