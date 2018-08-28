#####################################################################################
## CSC 315: Required packages
## Note #1: make sure that your .libPaths() is set correctly
##    if working from a school computer, e.g., using
## .libPaths("/Users/dancikg/OneDrive - Eastern Connecticut State University/Rlib")
## Note #2: it is recommended that these be installed from the R console
##    rather than from RStudio
#####################################################################################

## package for manipulating data frames ##
install.packages("dplyr")

## packages for knitting HTML documents ##
install.packages("knitr")
install.packages("yaml")
install.packages("htmltools")
install.packages("rmarkdown")

## package for plotting
install.packages("ggplot2")

## package for permutation function ##
install.packages("gtools")

## package for colorRampPalette function used to color heatmaps ##
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

