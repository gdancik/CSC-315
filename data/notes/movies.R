###########################################################
# Amazon movie reviews. The data used by this script is a
# processed version of data available from here:
# http://snap.stanford.edu/data/web-Movies.html
###########################################################

# read in file; file must be in your working directory or you must
# add the path to the file 
movies <- read.csv("movies_processed.txt")

# how many ratings are there?
nrow(movies)

# how many users are included?
length(levels(movies$UserID))

# how many movies are included?
length(levels(movies$productID))

###################################################################
# Let's summarize the ratings #
###################################################################
t = table(movies$Rating)
barplot(t / sum(t), main = "Relative frequency of ratings",
        xlab = "Rating", ylab = "Relative frequency", col = 1:5,
        ylim = c(0,0.60))
abline(h = 0)

############################################
# Get ratings for each product
############################################
s = split(movies$Rating, movies$productID)

#####################################
# Look at distribution of number of 
# ratings / product
#####################################
lengths = sapply(s, length)
hist(lengths, main = "# Ratings / Movie", ylab = "# ratings")

# what was rated the most frequently, and how many ratings? #
which.max(lengths)
max(lengths)

#####################################
# What are the top rated movies?
#####################################
avgs = sapply(s, mean)
sort(avgs, decreasing = TRUE)[1:10]

#####################################
# What are the top rated products?
# Require at least 200 ratings
#####################################
avg.200 = avgs[lengths >=200]
sort(avg.200, decreasing = TRUE)[1:10]

#############################################
### Where the $#%* is Shawshank Redemption? #
#############################################
m = match("B000P0J0EW", names(avgs))
s[m]    # all ratings are high
avgs[m] # average rating 

###############################################
### Customers who enjoyed Shawshank Redemption 
### also enjoyed.....
###############################################

keep = movies$productID%in%"B000P0J0EW"
users = movies$UserID[keep]   # users that rated Shawshank
keep = movies$UserID %in% users # movies rated by these users
ss.movies = movies[keep,] # movies rated by users that rated Shawshank

s = split(ss.movies$Rating, (ss.movies$productID))
lengths = lapply(s, length)
## compare: "0767803434" "0767804724" "0767811100"
boxplot(s[lengths > 1][1:5], col = "red")
