#########################################################################
# In this script we consider 2 normalization methods for
# RNA-seq data. However, these should not be used in practice
# to identify differentially expressed genes (even though they
# continue to be).

# RPKM - Reads per kilobase of transcript, per million mapped reads
# FPKM - Fragments per kilobase of transcript, per million mapped reads
# TPM - Transcripts per million

# RPKM/FPKM - for each gene, the number of reads/fragments you expect
#   to see for every 1000 bases in the gene, for however many million
#   fragments that were sequenced

# TPM - for each gene, the number of fragments you would expect to see
#   assuming that you sequenced 1 million full length genes
#
# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#########################################################################

library(dplyr)
library(edgeR) # this likely needs to be installed from Bioconductor

########################################################################
# Before we look at RNA-seq data, we need to cover the "sweep" 
# function. 

# sweep(m, MARGIN, STATS, FUN = "-")
#    This "sweeps" out a statistic on each row (when MARGIN = 1) or 
#    column (when MARGIN = 2). More specifically, the function (or
#    operation) is called with the row/column as the first argument
#    and the corresponding element of STATS as the second. The function
#    can be an operation ('-', '+', '/', '*').
########################################################################

# In the code below, we subtract the median from each row of m.
m <- matrix(1:20,ncol=4)
medians <- apply(m, 1, median)
s <- sweep(m, 1, medians, '-')
s

# normalize each row so that its values are between 0 and 1
maxes <- apply(m, 1, max)
s <- sweep(m, 1, maxes, '/')
s

############################################################
# Let's look at simulated RNA Seq data. We create a data 
# frame of read counts, assuming that the gene lengths are 
# specified in the gene_lengths vector below.
# Note the following:
#   - genes 1 and 2 have the same relative expression, but 
#   - gene 2 has twice as many reads as gene 1 (because 
#     gene 2 is twice as long)
#   - the total number of reads for the 2nd sample is roughly 
#     1.5x that of the first sample
#   - gene3 is differentially expressed (by roughly
#     a factor of 2)
############################################################

reads <- data.frame(s1 = c(10, 20, 10, 1)) 
rownames(reads) <- paste0("gene", 1:4)
reads <- reads %>% mutate(s2 = s1*1.5)
reads$s2[3] <- reads$s2[3]*2
reads$s2[4] <- 0
reads

####################################################################
# What measure can be used to compare genes.
# We will calculate RPKM values following the steps from
# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
####################################################################

# RPKM is the reads per kilobase of transcript, per million mapped reads
# 1. Count up the total reads in a sample and divide that number by 1,000,000 – 
#     this is our “per million” scaling factor.
# 2. Divide the read counts by the “per million” scaling factor. 
#    This normalizes for sequencing depth (total reads per sample), giving you 
#    reads per million (RPM)
# 3. Divide the RPM values by the length of the gene, in kilobases. 
#    This gives you RPKM.


# get total reads for each sample, per  million
N <- colSums(reads) / 1e6

# The length of genes 1 - 3 are specified below
gene_lengths <- c(100,200, 150, 125) 

# gene lengths in kilobases
L <- gene_lengths / 1000

# calculate reads per million mapped reads for each sample
RPM <- sweep(reads, 2, N, '/') 

# calculate reads per kilobase per transcript (gene), per million mapped reads
RPKM <- sweep(RPM, 1, L, '/') 
RPKM

# Conclusions
#   - Genes are now comparable WITHIN a sample:
#       - genes 1 and 2 have the same expression (within a sample)
#   - Genes are not comparable ACROSS or BETWEEN samples:
#       - the RPKM value for gene 1 is higher in sample 1 than in  
#         sample 2, but we know that this gene is not differentially expressed

# RPKM should NOT be used to compare genes across samples

# RPKM vs FPKM: RPKM is used when genome sequencing uses "single-end" reads. 
# When "paired-end" reads are used, sequencing is done in both directions, 
#  and a single molecule (the same one sequenced in both directions) 
#  can map to the gene. For paired-end reads, sequence assembly and 
#  quantification controls for this possible double counting; instead of 'reads', 
#  the term 'fragments' is used, and we call the normalization FPKM
#  (though the calculation is the same, just based on 'fragments', 
#  rather than 'reads')

#  Illustration of single vs paired-end reads:
#  https://www.yourgenome.org/sites/default/files/images/illustrations/bioinformatics_single-end_pair-end_reads_yourgenome.png

################################################################################
# TPM - transcripts per million
#   One reason whey RPKM/FPKM values are not comparable between samples
#   is because the total RPKM/FPKM values in a sample can differ. For 
#   an example, see RPKM-normalized counts table at:
#      https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#   This means that two RPKM values that are the same across samples does 
#   not mean that the gene is similarly expressed across samples

#   The TPM (transcripts per million) measure is an attempt to deal with
#   this. The TPM for a gene is its RPKM (or FPKM) value divided by the sum of all 
#   RPKM (or FPKM) values in the sample, times 1 million. 
################################################################################

N <- colSums(RPKM)
TPM <- sweep(RPKM, 2, N, '/') * 1e6

# Now the columns are normalized (they sum to the same value), but 
# we still are not able to compare across samples. Remember that only
# gene3 is differentially expressed, but we cannot see that from our
# data

# While FPKM data is readily available, and often used, more advanced
# methods should be used for processing count data when the goal is
# to find differentially expressed genes.

# We can also use the edgeR package to calculate RPKM/FPKM:
library(edgeR)
rpkm(reads, gene.length = gene_lengths)

