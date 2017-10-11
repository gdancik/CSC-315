####################################################
# Intro to hypothesis testing:
####################################################

library(ggplot2)
library(dplyr)
library(cowplot)

# H0: p = 1/7
# H1: p != 1/7

# Study finds p.hat = 15/54  (so n = 54) 

# What is the distribution of p.hat (the sample proportion) 
# under the null hypothesis???

p.hat <- 15/54   # the sample proportion
p <- 1/7         # the true proportion under the null hypothesis
n <-  54         # the sample size

mu <- p          # the expected value of the sample proportion
sd <- sqrt(p*(1-p)/n) # the sd of the sample proportion

df.phat <- data.frame(x = seq(mu-4*sd, mu+4*sd, length.out = 100)) %>%
           mutate(y = dnorm(x, mean = mu, sd = sd))
                      
plot.phat <- ggplot(df.phat) + geom_line(aes(x, y)) +
  theme_classic() + ggtitle("p.hat ~ N(1/7, sqrt((1/7)(6/7)/54)") +
  labs(x = "p.hat", y = "density") +
  geom_vline(xintercept = p.hat) + 
  geom_text(aes(x=p.hat-.04, y=6), label = "observed p.hat")
plot.phat

###############################################################
# The test statistic is, in general equal to
# Z = (observed.value - mean) / standard deviation
# For the sample proportion, its test statistic is
# Z = (p.hat - mu) / sd, where sd = sqrt(p0(1-p0)/n),
# and is Z ~ N(0,1), under H0
###############################################################
Z = (p.hat - mu) / (sd)
x = seq(-4,4, length.out = 100)

df.z <- data.frame(x = seq(-4, 4, length.out = 100)) %>%
  mutate(y = dnorm(x))

plot.z <- ggplot(df.z) + geom_line(aes(x, y)) +
  theme_classic() + ggtitle("Z ~ N(0,1)") +
  labs(x = "Z", y = "density") +
  geom_vline(xintercept = Z) +
  geom_text(aes(x=Z-.8, y=.25), label = "observed Z")
plot.z

# compare distributions of p.hat with observed proportion and the 
# standard normal distribution with the observed Z statistic
plot_grid(plot.phat, plot.z, nrow = 2)


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

# get boundary of the tails for the p-value; boundaries
# are symmetric about the mean
diff = abs(p.hat - p) 

par.orig = par(mfrow = c(2,1), mar = c(4,4,2,1.3))
shade.norm(p - diff, p + diff, mean = mu, sd = sd, main = "p.hat")
abline(v = c(mu+diff,mu-diff), col = "red", lwd = 2) 
shade.norm(-abs(Z), abs(Z), main = "Z")
abline(v = c(Z,-Z), col = "red", lwd = 2)

###########################################
# calculate p-value from test statistic Z
# the p-value is the area in the tails
###########################################
p.value = 2*pnorm(-abs(Z))

###########################################
# use prop.test
###########################################
p2 = prop.test(15,54, p = 1/7, correct = FALSE)
p2$statistic # prop.test calculates X-squared statistic = Z^2
sqrt(p2$statistic) # confirm that our statistics match

###########################################
# more accurate with continuity correction
###########################################
p3 = prop.test(15,54, p = 1/7, correct = TRUE)
p3$p.value

# What is the conclusion regarding the null hypothesis, in the context
# of this problem?
