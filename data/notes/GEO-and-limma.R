################################################################
# GEO-and-limma.R. This script demonstrates how to download 
# data from GEO and how to use the limma package to 
# identify differentially expressed probes
################################################################

library(GEOquery)

###########################################################
# Get the processed data for GSE1297. The object
# returned by getGEO is a LIST of AffyBatch objects
# getGEO returns a list because each GEO Series
# may contain multiple platforms
###########################################################
GSE1297 <- getGEO("GSE1297")

# GSE1297 is a list of length 1
typeof(GSE1297)
length(GSE1297)

# summary of GSE1297
GSE1297


###########################################################
# Pull out gene expression data and pheno type data
# (same functions as before, but GSE1297 is a list)
###########################################################

GSE1297.expr <- exprs(GSE1297[[1]])
GSE1297.p <- pData(GSE1297[[1]])

## is data normalized (is distribution similar across samples)? 

boxplot(GSE1297.expr, main = "processed data")

# If not, then take log2 (very important -- if data is not normalized,
# then downstream analysis will be WRONG. Processed data from GEO may
# or may not be normalized, so make sure to check)
GSE1297.expr <- log2(GSE1297.expr)
boxplot(GSE1297.expr, main = "log2 processed data")


#################################################################
# How many males and how many females are there?
#################################################################
gender <- as.character(GSE1297.p$characteristics_ch1.6)

################################################################
# Find differentially expressed (DE) probes between males
# and females. We will use limma, which requires us to design 
# a model.matrix using indicator
# variables and to specify the contrasts we are interested in
# (e.g., Females - Males)
###############################################################

library(limma)

# construct design matrix
design <- model.matrix(~0+gender)
head(design) # note that indicator variables are used

# let's change the column names
colnames(design) <- c("Female", "Male")

## limma package fits a linear model to each row of the expression matrix ##
fit <- lmFit(GSE1297.expr, design)

## for each probe, we now have the mean for each group as well as the 
## standard deviation
head(fit$coefficients)
head(fit$sigma)

## Specify the contrasts, which must match column names of design matrix ##
contrast.matrix <- makeContrasts(Male - Female,levels=design)

## fit model based on contrasts (e.g., Female - Male)
fit2 <- contrasts.fit(fit, contrast.matrix)
head(fit2$coefficients)
head(fit2$sigma)

# calculate moderate t-statistics by moderating standard errors
# toward a common value, which makes answers more robust
fit2 <- eBayes(fit2)

## get top probes, sorted by p-value (gives top 10 genes by default)
tt <- topTable(fit2,sort.by = "p")
tt

###############################################################
# let's confirm the top probe
###############################################################

# create data frame with expression values and gender
probe <- rownames(tt)[1]
m <- match(probe, rownames(GSE1297.expr))
df <- data.frame(expr = GSE1297.expr[m,], gender = gender)

# group data frame by gender, then summarize by expression value (new approach)
means <- df %>% group_by(gender) %>% summarize(mean = mean(expr))
means

diff(means$mean) # -4.675664
tt[1,]

## convert to FC ##
logFC <- tt[1,]$logFC
2**logFC

## visualize ##
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", probe, ", ", FC)

ggplot(df, aes(x = gender, y = expr, fill = gender)) + geom_boxplot() +
  ylab("log2 expression") + ggtitle(main) +
  scale_fill_manual(values = c("pink", "lightblue")) +
  theme_classic() + theme(legend.position = "none")


###############################################################
# How many genes have FDR < 5%?
###############################################################

## we need to set the following arguments:
## p.value - the adjusted p-value (FDR) cutoff (not the p-value)
## number - the maximum number of probes to return (should be
##          total number of probes in the dataset)
tt.05 <- topTable(fit2,sort.by = "p", p.value = 0.05, number = nrow(GSE1297.expr))
nrow(tt.05)

#######################################################################
# Aside: A closer look at the eBayes step above; eBayes 
# (empirical Bayes) pulls the estimated variance towards a common value
# (based on assumption that most genes are not differentially 
# expressed). This is visualized below.
#######################################################################

df <- data.frame(sigma.original = fit2$sigma, sigma.adjusted = sqrt(fit2$s2.post))
df <- arrange(df, sigma.original)
df <- mutate(df, x = 1:nrow(df))

ggplot(df) + geom_line(aes(x, sigma.original, color = "original")) +
             geom_line(aes(x, sigma.adjusted, color = "adjusted")) +
             theme_classic() +
             ggtitle("eBayes \"shrinks\" sigma towards a common value") +
             scale_color_manual(values = c("red", "black")) +
             geom_hline(yintercept = sqrt(fit2$s2.prior[1]), linetype = 2) +
             theme(legend.title = element_blank()) +
             xlab("index (smallest to highest)") + ylab("variance") + 
             annotate("text", x = 19000, y=.5, label = "common variance")



###############################################################
# Create a heatmap using the top probes (FDR < 0.05)
###############################################################
m <- match(rownames(tt.05), rownames(GSE1297.expr))
X <- GSE1297.expr[m,]
col.heat <- colorRampPalette(c("yellow", "blue"))(200)

# set colors for gender #
col.gender <- as.integer(as.factor(gender))
col.gender <- c("pink", "blue")[col.gender]

# clustering is done on original data, but rows are 
# scaled (converted to z-scores) by default for visualization
heatmap(X, ColSideColors = col.gender, col = col.heat)


# for clustering only, use 'dist' and 'hclust'
d <- dist(t(X)) # calculate distances between samples
h <- hclust(d)  # cluster the samples

#plot the clusters
label <- gsub("Sex: ", "", gender)
plot(h, label = label)


###############################################################
# Let's find the gene associated with the top probe
# This requires using microarray annotations available from
# GEO (see below) or bioconductor (not discussed) 
###############################################################
platform <- annotation(GSE1297[[1]])   

pl <- getGEO(platform)
pl <- Table(pl)


##########################################################
# Next, find the gene for the desired probe 
# The column names may vary by platform, but for GPL96,
# the probes are in the ID column and the genes are in 
# the "Gene Symbol" column
##########################################################
probe <- rownames(tt.05)[1]
m <- match(probe, pl$ID)
pl$`Gene Symbol`[m]

# get genes for top probes 
m <- match(rownames(tt.05), pl$ID)
pl$`Gene Symbol`[m]


###############################################################
# Sometimes we are interested in a gene and need to identify 
# corresponding probes. We will do this for RPS11. We use 
# 'grep' instead of 'match' because match only returns the index
# of the first match; grep returns the indices of all matches, 
# and also allows the use of regular expressions. The exact 
# nature of this search depends on the format of the Gene. 
###############################################################

# some basic regular expression syntax:
# "^a$" - an exact match of 'a' starting from the beginning ("^") and ending
#         at the end ("$") of a string
# "a|b" - matches either 'a' or 'b' anywhere in the string

x <- c("a", "a dog barked", "the cat meowed")
grep("a", x) # matches any occurence of 'a'
grep("^a$", x) # matches the single character 'a'
grep("dog|cat", x) # matches either 'dog' or 'cat'


# To find probes for RPS11, the code below does not work because it will
# also find any genes of the form *RPS11*, which includes, e.g.,
# MRPS11 and RPS11P1
g <- grep("RPS11", pl$`Gene Symbol`)
View(pl[g,])

# For this platform, the `Gene Symbol` column either contains a single
# gene name, or multiple gene names separated by " /// ". We
# therefore need to consider 4 possibilities (spaces are important):

# "^RPS11$" - exact match to single gene RPS11
# " RPS11 " - matches gene in middle of list, e.g., "X /// RPS11 /// Z"
# "^RPS11 " - matches gene that is first in the list, e.g., "RPS11 /// X"
# " RPS11$" - matches gene at the end of the list, e.g., "X /// RPS11"

g <- grep("^RPS11$| RPS11 |^RPS11 | RPS11$", pl$`Gene Symbol`)

# verify results
View(pl[g,])

# get probes
pl$ID[g]

##############################################################
# How can I find information about specific genes?
# A good reference for genes is the website Gene Cards
#  (http://www.genecards.org/)
##############################################################




