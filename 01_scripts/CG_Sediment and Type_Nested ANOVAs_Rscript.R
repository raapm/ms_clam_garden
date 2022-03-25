#**Analyzing Sediment with a Nested ANOVA**
#*
#install.packages("rmarkdown", lib = "C:/Program Files/R/R-4.0.2/library")
#cat("options(pkgType = 'binary')", file = "~/.Rprofile", sep = "\n", append = TRUE)

#
rm(list= ls()) 
# Here is my data set
#
cgsediment <- read.csv("cgsediment.csv")
#
#
# **Look for significant differences of different sediment characteristics
# and types**
# The following commands will install these packages if they are not already installed:
#if(!require(nlme)){install.packages("nlme")}
#if(!require(multcomp)){install.packages("multcomp")}
#if(!require(multcompView)){install.packages("multcompView")}
#if(!require(lsmeans)){install.packages("lsmeans")}
#if(!require(lme4)){install.packages("lme4")}
#if(!require(lmerTest)){install.packages("lmerTest")}
#if(!require(TukeyC)){install.packages("TukeyC")}
#
# **Nested anova with mixed effects model(nlme)**
library(nlme)
library(lme4)
#
cgsediment$Type = as.factor(cgsediment$Type)
#

model = lme(silt ~ Type, random = ~ 1|beach, data = cgsediment)
#
summary(model)
anova.lme(model)
#
# Analysis of random effect (beach)
#
model.fixed = gls(silt ~ Type, data = cgsediment)
anova(model, model.fixed)
#
# Post-hoc comparison of least-square means
#
library(multcomp)
#
posthoc = glht(model, linfct = mcp(Type = "Tukey"))
mcs = summary(posthoc,test = adjusted("single-step"))
mcs
cld(mcs, level = 0.05, decreasing = TRUE)
#
# Checking assumptions of model
hist(residuals(model))
rug(residuals(model))
plot(fitted(model), residuals(model))
plot(model)
#
# Mixed effects model with lmer
# 
library(lmerTest)
#
cgsediment$beach = as.factor(cgsediment$beach)
model = lmer(silt ~ Type + (1|beach), data = cgsediment, REML = TRUE)
anova(model)
#
#
# ** Nested anova with the aov function**
fit = aov(silt ~ Type + Error(beach), data = cgsediment)
summary(fit)
#
# Using Means Sq and Df values to get p-value for H = Type and Error = beach
pf(q= 5.556/6.056,df1=1,df2=4, lower.tail = FALSE)
#
# Shapiro-Wilk test for normality
#
shapiro.test(residuals(model))
#
# Bartlett's test for homogeneity of variance
#
bartlett.test(silt ~ interaction(Type,beach), data = cgsediment)
#
plot(cgsediment$Type, cgsediment$silt)
plot(cgsediment$beach, cgsediment$silt)
plot(cgsediment$surv, cgsediment$silt)
#

