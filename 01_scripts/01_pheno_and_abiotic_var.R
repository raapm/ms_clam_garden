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


#### 02. Statistical analysis ####
# Input data (to update: should use the same file for both PCA and AOV)
clam <- read.csv(file = input_AOV.FN)
head(clam) # this eventually will be replaced by the same input file as above
str(clam)

#### Survival by abiotic variables ####
# Add all variables here # (#TODO#, need to update Table S4)
#### TODO: do linear models for both survival and growth below ####
# Effect of carbon on survival, linear model
model1 <- lm(surv ~ carb, data = clam)
summary(model1)

# Effect of organics on survival, linear model
model2 <- lm(surv ~ org, data = clam)
summary(model2)

# Effect of srocks on survival
model3<- lm(surv ~ srocks, data = clam)
summary(model3)

# Effect of silt on survival
model4 <- lm(surv ~ silt, data = clam)
summary(model4)

# Effect of sand on survival
model5 <- lm(surv ~ sand, data = clam)
summary(model5)


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

