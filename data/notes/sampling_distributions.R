##################################################################
# sampling distributions
##################################################################

x = c(21,19,21,21,18)

## calculate population mean ##
mean(x)   ## population mean is 20

## draw histogram ##
xlim = c(min(x), max(x))
breaks = seq(17.5, 21.5, by= 1)
par.orig = par(mfrow = c(3,1), mar = c(3,3,2,1))
hist(x, prob = TRUE, xlim = xlim, 
     main = "Population distribution of x", breaks = breaks)
abline(v = mean(x), col = "red", lwd = 3)

###################################################################
## look at all combinations (order does not matter) of selecting 3 
## students
## set = FALSE because x contains duplicate values
## we do not repeat values; once a student is selected, he/she will
## not be selected again
###################################################################
library(gtools)
S = combinations(length(x), 3, x, set = FALSE, repeats.allowed = FALSE)
sample.means = rowMeans(S)
hist(sample.means, xlim = xlim, prob = TRUE, breaks = breaks,
     main = "Probability distribution of xbar when n = 3")
expected.value = mean(sample.means) ## 20, the same as the population mean!
abline(v = expected.value, col = "red", lwd = 3)

### repeat for n = 4 ###
library(gtools)
S = combinations(length(x), 4, x, set = FALSE, repeats.allowed = FALSE)
sample.means = rowMeans(S)
hist(sample.means, xlim = xlim, prob = TRUE, breaks = breaks,
     main = "Probability distribution of xbar when n = 4")
expected.value = mean(sample.means) ## 20, the same as the population mean!
abline(v = expected.value, col = "red", lwd = 3)

par(par.orig)


###################################################################
## constructs histogram and adds vertical red line for mean ##
###################################################################
plot.hist <-function(x, ...) {
  hist(x, ...)
  abline(v = mean(x), lwd = 3, col = "red")
}


###################################################################
## Central Limit Theorem when X ~ normally distributed
###################################################################

###################################################################
## returns the sample mean from 'n' randomly generated observations 
## from the standard normal distribution
###################################################################
get.sample.mean <-function(n) mean(rnorm(n))

x.population = rnorm(1000)
x.10 = replicate(5000, get.sample.mean(10))
x.30 = replicate(5000, get.sample.mean(30))
x.100 = replicate(5000, get.sample.mean(100))


xlim = c(min(x.population), max(x.population))
par.orig = par(mfrow = c(2,2), mar = c(3,3,2,1))
plot.hist(x.population, main = "Distribution of X", xlim = xlim)
plot.hist(x.10, main = "Distribution of sample mean (n = 10)", xlim = xlim)
plot.hist(x.30, main = "Distribution of sample mean (n = 30)", xlim = xlim)
plot.hist(x.100, main = "Distribution of sample mean (n = 100)", xlim = xlim)
par(par.orig)


###################################################################
## Central Limit Theorem when X ~ NOT normally distributed!
## we will give X the exponential distribution, which is skewed
###################################################################

## returns the sample mean from 'n' randomly generated observations from the
## exponential distribution
get.sample.mean <-function(n) mean(rexp(n))

x.population = rexp(1000)
x.10 = replicate(5000, get.sample.mean(10))
x.30 = replicate(5000, get.sample.mean(30))
x.100 = replicate(5000, get.sample.mean(100))


xlim = c(min(x.population), max(x.population))
par.orig = par(mfrow = c(2,2), mar = c(3,3,2,1))
plot.hist(x.population, main = "Distribution of X", xlim = xlim)
plot.hist(x.10, main = "Distribution of sample mean (n = 10)", xlim = xlim)
plot.hist(x.30, main = "Distribution of sample mean (n = 30)", xlim = xlim)
plot.hist(x.100, main = "Distribution of sample mean (n = 100)", xlim = xlim)
par(par.orig)


###################################################################
## Central Limit Theorem in practice:
## Suppose that heights have the distribution X ~ N(68, 1.7).
## We will calculuate probabilities involving the sample mean 
## from a sample of 20 individuals
###################################################################


###################################################################
## Visualization
###################################################################

x = seq(62, 74, by = .2)
y1 = dnorm(x, mean = 68, sd = 1.7) ## density of x
y2 = dnorm(x, mean = 68, sd = 1.7 / sqrt(20)) ## density of sample mean when n = 20
plot(x, y1, main = "Distribution of heights", type = "l", ylim = c(0, max(y1,y2)))
lines(x, y2, col = "red")
legend("topleft", c("X", "X.bar, n = 20"), col = c("black", "red"), lty = 1)

# Find P(X < 69)
pnorm(69, mean = 68, sd = 1.7)

# Find P(X > 67)
1 - pnorm(67, mean = 68, sd = 1.7)

# In a random sample of 20 individuals, find the probability
# that the sample mean is less than 69:
pnorm(69, mean = 68, sd = 1.7 / sqrt(20))

# In a random sample of 20 individuals, find the probability
# that the sample mean is greater than 67
1 - pnorm(67, mean = 68, sd = 1.7 / sqrt(20))

# In a random sample of 20 individuals, find the probability
# that the sample mean is between 67 and 67.5.
pnorm(67.5, mean = 68, sd = 1.7 / sqrt(20)) - pnorm(67, mean = 68, sd = 1.7 / sqrt(20))
