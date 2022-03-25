#'## 
#'## **Introduction**
#' 
#' 
#' This project is with my own data set on field work from May-August 2016. We transplanted little neck clams
#' onto clam garden (walled beach), and reference (non-walled beaches). We took sediment samples
#' to analyze the sediment for % carbonates, % organics and 8 different grain sizes, and measured growth
#' and survival after 16 weeks. 
#' 
#' Going to fit with a linear mixed effects model where beach type is the fixed
#' effect, and site is a random factor with plot nested in site.
#' Need the followint packages:
install.packages("nlme")
install.packages("multcomp")
install.packages("multcompView")
install.packages("lsmeans")
install.packages("lme4")
install.packages("lmerTest")
install.packages("TukeyC")


#' 
#'  
rm(list= ls()) 
#' Here is my data set:
#'
clam <- read.csv("clam.csv")
#'
#### Nested Anova with mixed effects model (nlme) ####
#' 
#' Since type is read in as a integer variable, convert it to factor
#' 
#clam$Type = as.factor(clam$Type)
clam$beach = as.factor(clam$beach)
#

#library(nlme)
#
#model = lme(surv ~ org, random = ~ 1|beach, data = clam)
#
#summary(model)
#anova.lme(model, type = "sequential", adjustSigma = FALSE)
##############################################
# without beach as random factor
#
#
# Use lm() or aov()
#
model1 <- aov(grow ~ carb, data = clam)
summary(model1)
#
plot(clam$grow,clam$carb)
#
modelm <- lm(grow ~ carb, data = clam)
summary(modelm)
# 
model2 <- aov(grow ~ org, data = clam)
summary(model2)
#
model2m <- lm(grow ~ org, data = clam)
summary(model2m)
#
model3<- lm(grow ~ rocks, data = clam)
summary(model3)
#
model4 <- lm(grow ~ srocks, data = clam)
summary(model4)
#
model5 <- lm(grow ~ vcsand, data = clam)
summary(model5)
#
model6 <- lm(grow ~ csand, data = clam)
summary(model6)
#
model7 <- lm(grow ~ sand, data = clam)
summary(model7)
#
model8 <- lm(grow ~ fsand, data = clam)
summary(model8)
#
model9 <- lm(grow ~ vfsand, data = clam)
summary(model9)
#
model10 <- lm(grow ~ silt, data = clam)
summary(model10)
#
plot(clam$beach,clam$sand)
plot(clam$beach,clam$carb)
plot(clam$beach,clam$surv)
#
