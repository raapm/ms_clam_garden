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
input.FN <- "02_input_data/cgsedimentPCA.csv" # the original input filename
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

# # Plot PCA
# P <- pca_fit %>%
#   augment(beaches) %>%                              # ?
#   rename_at(vars(starts_with(".fitted")),           # renames the PCA .fitted variables as .fittedPC1 -> PC1
#             list(~str_replace(.,".fitted",""))) %>% # ?
#   ggplot(aes(x=PC1, 
#              y=PC2, 
#              ))+
#   geom_point()
# P <- P + theme(panel.background = element_rect(fill = "white", colour = "black"))
# P <- P + theme(panel.grid = element_blank())
# P

# Plot biplot
str(beaches$Group) 

### Not clear if this is required ###
#beaches.groups <- c("A1","A2","A3","B1","B2","B3"),c("C1", "C2","C3","D1","D2","D3")
 #                   ,c("E1","E2","E3","F1","F2","F3")
### /END/ Not clear if this is required ###

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

#### Remaining questions
# Q1: what is 'X' column? It is highly correlated with growth
# Q2: what are the 'Groups' ? Is there anything that defines this other than the PCA itself? 
# Q3: can we use the same data input as the other scripts? only one data input from raw data


#### 02. Statistical analysis ####
# Input data (to update: should use the same file for both PCA and AOV)
clam <- read.csv(file = input_AOV.FN)
head(clam) # this eventually will be replaced by the same input file as above
str(clam)

#### Nested ANOVA with mixed-effects model (nlme) ####
# Since type (CG or Ref) is read in as a character variable, convert it to factor
#clam$Type = as.factor(clam$Type)

# Since beach () is read in as a character variable, convert it to factor
clam$beach = as.factor(clam$beach)

#model = lme(surv ~ org, random = ~ 1|beach, data = clam)
#
#summary(model)
#anova.lme(model, type = "sequential", adjustSigma = FALSE)
# without beach as random factor

# Effect of carbon on survival
model1 <- aov(surv ~ carb, data = clam)
summary(model1)

# Effect of carbon on survival, linear model
modelm <- lm(surv ~ carb, data = clam)
summary(modelm)

# Effect of organics on survival, AOV
model2 <- aov(surv ~ org, data = clam)
summary(model2)

# Effect of organics on survival, linear model
model2m <- lm(surv ~ org, data = clam)
summary(model2m)

# Effect of 'srocks' on survival
model3<- lm(surv ~ srocks, data = clam)
summary(model3)

# Effect of silt on survival
model4 <- lm(surv ~ silt, data = clam)
summary(model4)

# Effect of sand on survival
model5 <- lm(surv ~ sand, data = clam)
summary(model5)


#### Questions. models ####
# Q5. How to determine the most significantly associated variable? Stepwise process of removing variables

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
