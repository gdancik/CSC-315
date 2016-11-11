################################################################
# limma.R. This script downloads processed data from GEO
# (using the GEOquery from Bioconductor), and uses the limma
# package to identify differentially expressed probes
################################################################

library(GEOquery)
###########################################################
# Get the processed data for GSE1297. The object
# returned by getGEO is a LIST of AffyBatch objects
# getGEO returns a list because each GEO Series
# may contain multiple platforms
###########################################################
GSE1297 = getGEO("GSE1297")
GSE1297
length(GSE1297)

###########################################################
# Pull out gene expression data and pheno type data
###########################################################
GSE1297.expr = exprs(GSE1297[[1]])
GSE1297.p = pData(GSE1297[[1]])

## is data normalized? If not, then take log2
boxplot(GSE1297.expr, main = "processed data")

GSE1297.expr = log2(GSE1297.expr)
boxplot(GSE1297.expr, main = "log2 processed data")



#################################################################
# How many males and how many females are there?
#################################################################
gender = as.character(GSE1297.p$characteristics_ch1.1)

################################################################
# Finding differentially expressed probes
################################################################

###############################################################
## DE between males and females: We will use limma, which
## requires us to design a model.matrix using indicator
## variables and to specify the contrasts we are interested in
## (e.g., Females - Males)
###############################################################

library(limma)



design = model.matrix(~0+gender)
head(design) # note that indicator variables are used

# let's change the column names
colnames(design) = c("Female", "Male")

## limma package fits a linear model to each row of the expression matrix ##
fit = lmFit(GSE1297.expr, design)

## for each probe, we now have the mean for each group as well as the 
## standard deviation
head(fit$coefficients)
head(fit$sigma)

## Specify the contrasts, which must match column names of design matrix ##
contrast.matrix <- makeContrasts(Female - Male,levels=design)

## fit model based on contrasts (e.g., Female - Male)
fit = contrasts.fit(fit, contrast.matrix)
head(fit$coefficients)
head(fit$sigma)

# calculate moderate t-statistics by moderating standard errors
# toward a common value, which makes answers more robust
fit = eBayes(fit)

## get top probes, sorted by p-value (gives top 10 genes by default)
tt = topTable(fit,sort.by = "p")
tt

###############################################################
# let's confirm the top probe
###############################################################
probe = rownames(tt)[1]
m = match(probe, rownames(GSE1297.expr))
s = split(GSE1297.expr[m,], gender)
boxplot(s, col = c("pink", "lightblue"), ylab = "log2 expression",
        main = probe)

## logFC should match limma table ##
means = sapply(s, mean)

means[1] - means[2] # 4.675664
tt[1,]

## convert to FC ##
logFC = tt[1,]$logFC
2**logFC


###############################################################
# How many genes have FDR < 5%?
###############################################################

## p.value cutoff of 0.05, based on all probes ##
## Note: the p.value argument corresponds to the adjusted p-value (FDR), and
##    NOT the actual p-value
tt.05 = topTable(fit,sort.by = "p", p.value = 0.05, number = nrow(GSE1297.expr))
nrow(tt.05)

###############################################################
# Create a heatmap using the top probes (FDR < 0.05)
###############################################################
m = match(rownames(tt.05), rownames(GSE1297.expr))
X = GSE1297.expr[m,]
col.heat = colorRampPalette(c("yellow", "blue"))(200)

# set colors for gender #
col.gender = as.integer(as.factor(gender))
col.gender = c("pink", "blue")[col.gender]

# clustering is done on original data, but rows are 
# scaled (converted to z-scores) by default for visualization
heatmap(X, ColSideColors = col.gender, col = col.heat)


d = dist(t(X)) # calculate distances between samples
h = hclust(d)  # cluster the samples

#plot the clusters
label = gsub("gender: ", "", gender)
plot(h, label = label)


###############################################################
# Let's find the gene associated with the top probe
# If you look at the raw data annotation this will be hgu133a
# for GEO processed data, this is GPL96 (same as hgu133a)
###############################################################
annotation(GSE1297[[1]])   

## load database of Affymetrix HGU133A annotation
library(hgu133a.db)

# genes are stored in a list, with element names corresponding to probes
genes.db = as.list(hgu133aSYMBOL)
head(genes.db)

# find the gene for the desired probe 
probe = rownames(tt.05)[1]
m = match(probe, names(genes.db))
genes.db[m]

###############################################################
# Sometimes we start with the gene and need to identify the
# corresponding probes. We will do this for XIST. The 'grep'
# function should be used with a '^' at the beginning and a '$' 
# at the end of the gene name to require exact match
# grep("XIST") would match any gene of the form "*XIST*"
# grep("^XIST$") finds exact matches of the gene
###############################################################
g = grep("^XIST$", genes.db)   
genes.db[g]
