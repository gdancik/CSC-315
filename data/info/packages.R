n#####################################################################################
# CSC 315: Required packages
# Run this script to install the packages we will be using this semester
# These packages are already installed on our classroom computers
#####################################################################################

# package for data manipulation and visualization
install.packages('tidyverse')

# package to combine multiple ggplot2 plots into 1
install.packages('cowplot')

# package for knitting HTML documents
install.packages('rmarkdown')

#package for permutations
install.packages('gtools')

## package for colorRampPalette function used to color heatmaps ##
install.packages('RColorBrewer')

#################################################
# Bioconductor packages for Microarray Analysis
#################################################

# Bioconductor packages are now installed through BiocManager
install.packages('BiocManager')

# load BiocManager and install the packages
library(BiocManager)
BiocManager::install('affy')
BiocManager::install('affydata')
BiocManager::install('leukemiasEset')
BiocManager::install('GEOquery')
BiocManager::install('limma')
BiocManager::install('hgu133a.db')

