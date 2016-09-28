#1
tt = cbind("Not" = c(26,117, 172), "Pretty Happy" 
           = c(233,473,383), "Very Happy" = c(164,293,132))

rownames(tt) = c("Above average", "Average", "Below average")
tt=tt / rowSums(tt)
tt



#9 - 12
cereal = read.delim("http://pastebin.com/raw/0G6DrHyC")

plot(cereal$sugars, cereal$calories, xlab = "sugar content (g)", 
     ylab = "calories", pch = 19, 
     main = "Caloric and sugar content in various cereals")
cor(cereal$sugars, cereal$calories, use = "complete")
fit = lm(calories ~ sugars, data = cereal)
abline(fit, col = "blue")

# find and interpret y-intercept and slope

# predict calories for 10 and 25 grams of sugar
predict(fit, data.frame(sugars = 10))

# find residual for observation of 5 grams of sugar and 70 calories
# residual = observed - predicted
residual = 70 - predict(fit, data.frame(sugars = 5))
residual
