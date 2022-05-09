# Physical data and survival data
# Initialized to GitHub 2022-03-25

# clear workspace
#rm(list=ls())

## Install packages and load libraries (PCA)
# install.packages("ggrepel")
# install.packages("tidyverse")
# install.packages("broom")
# install.packages("devtools")
# install.packages("ggbiplot") # not available with current version of R (v.4.1.2), requires install from github
library("ggrepel")
library("tidyverse")
library("broom")
library("devtools")
# install_github("vqv/ggbiplot")
library("ggbiplot")

## Install packages and load libraries (statistical analysis)
# install.packages("nlme")
# install.packages("multcomp")
# install.packages("multcompView")
# install.packages("lsmeans")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("TukeyC")
library("nlme")

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

# Set input filenames
#input.FN <- "02_input_data/cgsedimentPCA.csv" # the original input filename
input.FN <- "02_input_data/cgsedimentPCA_from_Monique_2022-04-27.csv" # new, with initial weight/ height data
input_AOV.FN <- "02_input_data/clam.csv"      # the original input filename

# Load data
beaches <- read.csv(file = input.FN)
head(beaches)
str(beaches)

row.names(beaches) <- beaches[,"Beach"] # provide row names
head(beaches)


#### 01. PCA based on physical characteristics ####
# PCA based on growth and physical attributes using prcomp
pca_fit <- beaches %>%
  select(where(is.numeric)) %>%  # only selects numeric columns
  scale() %>%                    # scales the numeric columns
  prcomp()

summary(pca_fit)

# Plot biplot
str(beaches$Group) 

# Plot PCA (ggplot)
pdf(file = "03_pheno_results/per_plot_abiotic_PCA.pdf", width = 6, height = 6)
par(mfrow = c(1,1))
p <- ggbiplot(pcobj = pca_fit, labels = row.names(beaches)
              , color = 'black', varname.adjust = 1.1, varname.size = 2.9
              , groups = beaches$Group, ellipse = TRUE
              ) + expand_limits(x = c(-2, 2), y = c(-2.5,2))

p # uncomment if want to keep grey grid

# optional to remove grey grid
# p <- p + theme(panel.background = element_rect(fill = "white", colour = NA) )
# p <- p + theme(panel.grid = element_blank(), legend.position = "none")
# p
dev.off()


#### 02. Statistical analysis, explanatory variable models ####
# Input data (TODO: should use the same file for both PCA and AOV)
clam <- read.csv(file = input_AOV.FN)
head(clam) # this eventually will be replaced by the same input file as above
#str(clam)

# Remove unneeded cols
drop_cols <- c("cage", "beach_label", "orig")
clam <- clam[, !(colnames(clam) %in% drop_cols)]
rm(drop_cols) # clean enviro
head(clam) # this eventually will be replaced by the same input file as above


#### 02.1 Abiotic variable effect on growth and survival ####
# Linear models of all variables on survival and growth

# Which columns are to be focused on as putative explanatory variables? 
non_explan_cols <- c("beach", "grow", "surv")

explan_vars <- colnames(
                         clam[, !(colnames(clam) %in% non_explan_cols)]
                        )
rm(non_explan_cols) # clean enviro

# Set nulls, Loop
abiotic_fx.list <- list(); voi <- NULL

for(i in 1:length(explan_vars)){
  
  voi <- explan_vars[i]
  
  # Linear model per explan variable, survival
  abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]] <- lm(clam[, "surv"] ~ clam[, voi])
  abiotic_fx.list[[paste0("surv_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0(voi, "surv.mod")]])
  
  # Linear model per explan variable, growth
  abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]] <- lm(clam[, "grow"] ~ clam[, voi])
  abiotic_fx.list[[paste0("grow_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0(voi, "grow.mod")]])
  
}

# Then write out the output
capture.output(abiotic_fx.list, file = "03_pheno_results/abiotic_variables_on_grow_surv_models.txt")


#### Original analysis of Growth by Sediment ####
# Fit a linear mixed effects model 
# Fixed effect: beach type; Random factor: site, with plot nested in site

# Required packages
# install.packages("nlme")
# install.packages("multcomp")
# install.packages("multcompView")
# install.packages("lsmeans")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("TukeyC")

library("nlme")
library("multcomp")
library("multcompView")
library("lsmeans")
library("lme4")
library("lmerTest")
library("TukeyC")

#### Nested Anova with mixed effects model (nlme) ####
#  Since type is read in as a integer [ no: character ] variable, convert it to factor
 
#clam$Type = as.factor(clam$Type)
#clam$beach = as.factor(clam$beach) ### TEST LATER TO SEE IF THIS AFFECTS ANYTHING  (BJGS)

# THIS APPEARS TO BE THE NESTED WITH RANDOM VARIABLE, CAME COMMENTED
#model = lme(surv ~ org, random = ~ 1|beach, data = clam)
#summary(model)
#anova.lme(model, type = "sequential", adjustSigma = FALSE)

# Without beach as random factor # THIS WAS THE UNCOMMENTED OPTION (BJGS)
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

#### Are sediment variables affected by Clam Gardens? ####
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

#### /END/ Are sediment variables affected by Clam Gardens? ####



#### Correlation of variables ####
#install.packages("corrplot")
library("corrplot")
colnames(clam)
to_cor_phenos.vec <- c("grow", "surv", "carb", "org", "rocks", "srocks", "vcsand", "csand", "sand"
                       , "fsand", "vfsand", "silt")
cor.set <- cor(clam[, to_cor_phenos.vec]
               , use = "pairwise.complete.obs" # cor bw ea pair var is done with all pairs of obs on those var 
               # note: only works with "pearson" method
               )

pdf(file = "03_pheno_results/surv_grow_abiotic_correlations.pdf", width = 4, height = 4)
par(mfrow = c(1,1))
corrplot(cor.set
         , method = "circle"
         #, order="hclust"
         , addshade = "all"
         #, type = "lower"
         #, outline = T
         , tl.col = "black"
         , tl.cex = 0.9
         ) #plot matrix
dev.off()

#### Boxplot of survival ####
pdf(file = "03_pheno_results/surv_grow_by_beach.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
#boxplot(clam$surv ~ clam$Type)
boxplot(clam$surv ~ clam$beach, las = 1
        , xlab = "Beach", ylab = "% Survival"
        )

#boxplot(clam$grow ~ clam$Type)
boxplot(clam$grow ~ clam$beach, las = 1
        , xlab = "Beach", ylab = "Growth"
)
dev.off()

# Consider survival and carbonates
pdf(file = "03_pheno_results/survival_by_carbonate.pdf", width = 8, height = 5)
par(mfrow=c(1,1))
plot(x = clam$carb, y = clam$surv
     , ylab = "Survival (%)"
     , xlab = "Carbonate (%)"
)

mod <- lm(surv ~ carb, data = clam)
summary(mod)

results <- summary(mod)
text(x = 13, y = 85, labels = paste0("adj. Rsq. = ", round(results$adj.r.squared, digits = 2)))
text(x = 13, y = 75, labels = paste0("p-value: ", round(results$coefficients["carb","Pr(>|t|)"], digits = 5)))
dev.off()

##### New information ####
plot(x = clam$inwt, y = clam$surv)
mod.inwt.surv <- lm(formula = clam$inwt ~ clam$surv)
mod.surv.inwt <- lm(formula = clam$surv ~ clam$inwt)
summary(mod.inwt.surv)
summary(mod.surv.inwt)

# Carbonate? 
plot(x = clam$inwt, y = clam$carb)
plot(y = clam$inwt, x = clam$carb)
mod1 <- lm(formula = clam$carb ~ clam$inwt)
summary(mod1)

mod.surv.carb <- lm(formula = clam$surv ~ clam$carb)
summary(mod.surv.carb)

boxplot(clam$surv ~ clam$beach)
boxplot(clam$inwt ~ clam$beach)

mod.surv.carb.inwt <- lm(formula = clam$surv ~ clam$carb + clam$inwt)
summary(mod.surv.carb.inwt)

par(mfrow=c(2,2))
plot(x = clam$carb, y = clam$surv)
plot(x = clam$inwt, y = clam$surv)
plot(x = clam$inwt, y = clam$carb)
plot(x = clam$carb, y = clam$grow)
mod.grow.by.carb <- lm(formula = clam$grow ~ clam$carb)
summary(mod.grow.by.carb)


plot(x = clam$inwt, y = clam$grow)

mod.grow.by.inwt <- lm(formula = clam$grow ~ clam$inwt)
summary(mod.grow.by.inwt)


plot(x = clam$grow, y = clam$surv)
#plot(y = clam$grow, x = clam$surv)
mod.grow.by.surv <- lm(formula = clam$grow ~ clam$surv)
summary(mod.grow.by.surv)

#### Split beaches by weight group ####
par(mfrow = c(1,2))
plot(x = clam$carb[clam$beach=="C" | clam$beach=="E" | clam$beach=="F"]
      , y = clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F"]
      , xlab = "carb (low weight beaches)"
       , ylab = "survivorship")

mod <- lm(formula = clam$carb[clam$beach=="C" | clam$beach=="E" | clam$beach=="F"]
          ~ clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F"]
          )
summary(mod)

plot(x = clam$carb[clam$beach=="A" | clam$beach=="B" | clam$beach=="D"]
     , y = clam$surv[clam$beach=="A" | clam$beach=="B" | clam$beach=="D"]
     , xlab = "carb (high weight beaches)"
     , ylab = "survivorship")

mod <- lm(formula = clam$carb[clam$beach=="A" | clam$beach=="B" | clam$beach=="D"]
          ~ clam$surv[clam$beach=="A" | clam$beach=="B" | clam$beach=="D"]
)
summary(mod)

par(mfrow = c(1,1))
plot(y = clam$carb, x = clam$surv)

