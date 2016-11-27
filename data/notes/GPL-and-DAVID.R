###############################################################
# GPL and DAVID
###############################################################

###############################################################
# Let's get the GSE1297 data again (modified from limma.R)
###############################################################
library(GEOquery)
GSE1297 = getGEO("GSE1297")

###########################################################
# Pull out gene expression data and pheno type data
###########################################################

# expression data from this dataset must be logged
GSE1297.expr = exprs(GSE1297[[1]])
GSE1297.expr = log2(GSE1297.expr)
GSE1297.p = pData(GSE1297[[1]])


################################################################
# Find differentially expressed probes between males and
# females
################################################################

library(limma)
gender = as.character(GSE1297.p$characteristics_ch1.1)
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

## get top 50 probes (for demonstration only), sorted by p-value
tt = topTable(fit,sort.by = "p", number = 50)

## get platform for this dataset (the 1st dataset in the list)
platform = annotation(GSE1297[[1]])

## get platform data ##
GPL = getGEO(platform)

#####################################################################
## you will need to view the platform data to determine which column
## contains the probe IDs (usually the ID column) and which column
## includes the gene symbol column (in this case, 'Gene Symbol')
#####################################################################
genes.db = Table(GPL)

## find the genes corresponding to the DE probes ##
m = match(rownames(tt), genes.db$ID)
genes = as.character(genes.db$'Gene Symbol'[m])

## only keep probes that correspond to genes ##
keep = genes!=""
genes = genes[keep]

## if more than one name is given, keep only the first ##
genes = strsplit(genes, " /// ")
get.first <-function(x) {
  x[1]
}
genes = sapply(genes, get.first)

# get a unique set of genes 
genes = unique(genes) 


## view genes for input into DAVID (we will use version 6.7), 
## select Functional Annotation -> Upload gene list, 
## set identifier to OFFICIAL_GENE_SYMBOL, select Gene List.
## genes may be saved to a file by setting the file argument (which 
## will save the file in your current working directory - see getwd() - 
## by default)
write.table(genes, row.names = FALSE, quote = FALSE)

#####################################################################
## we will look at GOTERM_BP_FAT under Gene_Ontology and KEGG_PATHWAYS
## under Pathways
#####################################################################