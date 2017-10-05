##################################################################
# The normal distribution
##################################################################

library(ggplot2)
library(dplyr)
library(cowplot)


##################################################################
# create a vector of heights
##################################################################
heights <- c(66.1759117740307, 66.1118375828377, 69.3398551802815, 70.4991167238977, 
    70.7354713429583, 66.8777424499321, 66.7103751135044, 67.998337409316, 
    69.5210869064879, 68.8470644101842, 70.0596931565498, 71.0154419501861, 
    65.7861428292448, 64.7948285607431, 66.9861649481436, 67.5283004226946, 
    66.6906302756055, 70.2184193623913, 70.6098030390015, 63.3569739129227, 
    68.534978119765, 67.3597759440788, 69.305240284645, 67.8014644375226, 
    65.8455616472892, 73.2715698159093, 69.0255738087185, 66.5503992464326, 
    68.2018131042349, 63.2807987441221, 67.0816561689806, 66.7636378663116, 
    67.8580308828024, 68.1618434342104, 66.6704530072032, 68.8866603679841, 
    66.9959315571359, 67.2598262605613, 66.7716288824807, 71.3268449901713, 
    65.7684215531473, 68.9392236496006, 68.8075214177273, 66.9828251741788, 
    67.0585081777219, 63.2197251962368, 67.8334338291738, 69.1256626969062, 
    67.098802016255, 71.8537548116061, 68.5077939511076, 67.4631190657698, 
    68.7682657975886, 68.4389671592818, 65.1083655300538, 70.6188289457821, 
    66.2051076389411, 68.6718172959333, 68.8217079591491, 68.7220495841377, 
    66.2941801179521, 67.203845719609, 66.6360550720183, 67.5497892468022, 
    63.7628983953615, 65.576280689069, 64.0067204772203, 67.9575607594972, 
    66.1883614584798, 70.3170715975757, 69.6617964114646, 63.828494096075, 
    66.9334761277201, 70.4796907579624, 67.2801882801361, 62.4496489533135, 
    68.0330782726793, 65.5903445506025, 67.4039768511255, 66.2852782333007, 
    67.2000392075503, 71.4179254017927, 66.8235099465514, 72.0666556398044, 
    63.7242656432218, 69.683115643616, 68.5226822988891, 69.8651235405071, 
    70.3909690825223, 68.8063397849254, 69.5386113193661, 64.9229208483115, 
    67.9294713565318, 72.1968849311748, 65.3545159580504, 64.3705367578319, 
    65.1015079618489, 66.4784003830516, 66.2087925112091, 67.7854148025578
  )

##################################################################
# This function constructs a histogram of X while optionally shading
# all values corresponding to P(X <= obs), and optionally adding
# a density curve if density = TRUE
# ... are additional arguments to hist
##################################################################
shaded.hist <- function(X, hist = TRUE, density = FALSE, obs = NULL, ...) {
  h = hist(X, plot = FALSE)
  col = rep("white", length(h$breaks))
  ylim = NULL
  if (density) {
    m = max( max(h$density), dnorm(mean(X), mean = mean(X), sd = sd(X)))
    m = m +0.1*m  
    ylim = c(0,m)
  }
  
  if (hist & !is.null(obs)) {
     col[h$breaks<obs] = "red"
  }
  lty = 1
  if (!hist) lty = 0
  hist(X, prob = TRUE, ylim = ylim, col = col, lty = lty, ...)
  
  if (density) {
    curve(dnorm(x, mean = mean(X), sd = sd(X)), col = "blue", lty = 2, lwd = 2, add = TRUE)
  }
}


#################################################
## If an individual is randomly selected, what 
## is the probability their height is less than
## or equal to 67 inches?
#################################################

hist(heights, prob = TRUE, main = "Histogram of Heights", xlab = "height")
p = sum(heights <=67) / length(heights)

#####################################################
## We can visualize this on a histogram. Let
## X be the height of a randomly selected individual
## P (X <= 67) is the cumulative density of
## the red bars, while P (X > 67) is the 
## cumulative density of the white bars, or
## 1 - P(X <=67)
#####################################################
shaded.hist(heights, obs = 67, main = "Histogram of Heights", xlab = "height")


#####################################################
## What is the probability P (67 < X < 70)?
## We can take:  P(X <=70) - P(X <= 67)
#####################################################

par.orig = par(mfrow = c(2,1), mar = c(3,3,2,1))
shaded.hist(heights, obs=70, main = "Histogram of Heights", xlab = "height")
shaded.hist(heights, obs=67, main = "Histogram of Heights", xlab = "height")
par(par.orig)

###################################################################
# In practice, probability distributions are visualized by curves
# which (approximately) reflect the distribution 
# (shape of the data) and density across the sample space.
# The cumulative density is equivalent to the area under the curve
###################################################################

shaded.hist(heights, density = TRUE, obs = 70, main = "Histogram of Heights With Normal Approximation", xlab = "height")

###################################################################
# Several distributions are given below. However we will focus
# on the most widely used distribution, the normal distribution
# (or bell curve)
###################################################################


# set up data frames with (x,y) values for various distributions
# for given 'x' values, use 'mutate' to add the corresponding 'y' value (density)
data.unif <- data.frame(x = seq(0,1,by=.1)) %>% mutate(y = dunif(x))
data.exp <- data.frame(x = seq(0,10,by=.1)) %>% mutate(y = dexp(x))
data.chisq <- data.frame(x = seq(0,10,by=.1)) %>% mutate(y = dchisq(x, df = 3))
data.norm <- data.frame(x = seq(0,10,by=.1)) %>% mutate(y = dnorm(x, mean = 5, sd = 1))

# generate each plot
plot.unif <- ggplot(data.unif, aes(x, y)) + geom_line(color = "blue", lty = 2) +
  ggtitle("Uniform Distribution") + labs(x = "x", y = "density") + theme_classic()
plot.exp <- ggplot(data.exp, aes(x, y)) + geom_line(color = "blue", lty = 2) +
  ggtitle("Exponential Distribution") + labs(x = "x", y = "density") + theme_classic()
plot.chisq <- ggplot(data.chisq, aes(x, y)) + geom_line(color = "blue", lty = 2) +
  ggtitle("Chi-Square Distribution") + labs(x = "x", y = "density") + theme_classic()
plot.norm <- ggplot(data.norm, aes(x, y)) + geom_line(color = "blue", lty = 2) +
  ggtitle("Normal Distribution") + labs(x = "x", y = "density") + theme_classic()

# display multiple plots on a grid
plot_grid(plot.unif, plot.exp, plot.chisq, plot.norm)



###################################################################
# The normal distribution is completely specified by its mean
# and standard deviation
###################################################################

# normal distribution with mean = 60 and sd = 72
d1 <- data.frame(x = seq(60,72,by=.1)) %>% mutate(y=dnorm(x, mean = 66, sd = 2))
p1 <- ggplot(data = d1, aes(x, y)) + geom_line() +
  ggtitle("Normal Distribution\n(mu = 66, sigma = 2)")  + ylab("density")

# normal distribution with mean = 0 and sd = 1 (i.e., standard normal distribution)
d2 <- data.frame(x = seq(-3,3,by=.1)) %>% mutate(y=dnorm(x, mean = 0, sd = 1))
p2 <- ggplot(data = d2, aes(x, y)) + geom_line() +
  ggtitle("Standard Normal Distribution\n(mu = 0, sigma = 1)")  + ylab("density")

# plot both distributions
plot_grid(p1, p2, nrow = 2)

########################################################################
# Suppose heights are normally distributed with mean = 68 and 
# sd = 1.7 inches. Find the probability that a randomly selected 
# person is less than (or equal to) 70 inches tall 
#
# pnorm(a, mean = mu, sd = sd) returns P(X <= a) when X ~ N(mu, sd)
#
########################################################################

pnorm(70, mean = 68, sd = 1.7)


########################################################################
# Shades the area under the normal curve between a and b 
# This area corresponds to P (a <= x <= b) when X ~ N (mean, sd)
# For P (x <= b) set a = -Inf
# For P (x >= a) set b = Inf
########################################################################
shade.norm <- function(a,b, mean = 0, sd = 1,  ...) {
  m1 = mean-4*sd
  if (!is.infinite(a)) m1 = min(a,m1)
  m2 = mean+4*sd
  if (!is.infinite(b)) m2 = max(b,m2)
  
  x = seq(m1,m2, length.out = 100)
  plot(x, dnorm(x, mean = mean, sd = sd), type = "l", ylab = "normal density", ...)
  if (is.infinite(a)) a = m1
  if (is.infinite(b)) b = m2
  r = seq(a,b,length.out=100)
  cord.x <- c(a,r,b)
  cord.y <- c(0,dnorm(r, mean = mean, sd = sd),0)
  polygon(cord.x,cord.y,col='skyblue')
  abline(h = 0)

  p = pnorm(b, mean = mean, sd = sd) - pnorm(a, mean = mean, sd = sd)
  legend("topleft", legend = paste("p = ", round(p,3), sep = ""))
  
}


########################################################################
# Assume that X ~ N(68, 1.7).
# 1. Find P(X <= 68)
# 2. Find P(X <= 65)
# 3. Find P(X > 65)
# 4. Find P(64.6 <= X <= 71.4)
########################################################################

#1. P(X <= 68)
shade.norm(-Inf, 68, mean = 68, sd = 1.7, 
           main = "X ~ N(68,1.7)\nP(X <= 68)")  # or
pnorm(68, mean = 68, sd = 1.7)

#2. P(X <= 65)
shade.norm(-Inf, 65, mean = 68, sd = 1.7,
           main = "X ~ N(68, 1.7)\nP(X<=65)")  # or
pnorm(65, mean = 68, sd = 1.7)

#3. P(X > 65)
shade.norm(65, Inf, mean = 68, sd = 1.7,
           main = "X ~ N(68, 1.7)\nP(X>65)")  # or
1 - pnorm(65, mean = 68, sd = 1.7)

#4. P(64.6 <= X <= 71.4)
shade.norm(64.6, 71.4, mean = 68, sd = 1.7,
           main = "X ~ N(68, 1.7)\nP(64.6 <= X <= 71.4)")  # or
pnorm(71.4, mean = 68, sd = 1.7) - pnorm(64.6, mean = 68, sd = 1.7)


########################################################################
# Percentiles (quantiles) from the normal distribution
# Assume still that X ~ N(68, 1.7)
# What is the 75th percentile?
# What is the 90th percentile?
########################################################################

qnorm(.75, mean = 68, sd = 1.7)
qnorm(.90, mean = 68, sd = 1.7)

########################################################################
# Verify 90th percentile
########################################################################
q = qnorm(.90, mean = 68, sd =1.7)
shade.norm(-Inf, q, mean = 68, sd = 1.7, 
           main = "90th percentile of X ~ N(68, 1.7)")


########################################################################
# Suppose that X is normally distributed. Find the probabilty
# that a randomly selected value from X is more than 2 standard
# deviations above the mean.
########################################################################
