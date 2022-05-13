# Physical data and survival data
# Initialized to GitHub 2022-03-25

# clear workspace
#rm(list=ls())

#### Front Matter ####
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
#install.packages("corrplot")

## Install packages and load libraries (statistical analysis)
# install.packages("nlme")
# install.packages("multcomp")
# install.packages("multcompView")
# install.packages("lsmeans")
# install.packages("lme4")
# install.packages("lmerTest")
# install.packages("TukeyC")
library("nlme")
library("lme4")
library("multcomp")
#library("lmerTest") # has to be loaded later to load over lme4 when needed
library("corrplot")
library("dplyr")

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)


#### 00. Load data ####
# Set input filenames
input.FN <- "02_input_data/cg_sediment_phenos_2022-05-12.csv" # new, all-in-one file

# TODO: delete these lines, as we now use one single input file #
#input.FN <- "02_input_data/cgsedimentPCA_from_Monique_2022-04-27.csv" # new, with initial weight/ height data
#input.FN <- "02_input_data/cgsedimentPCA.csv" # the original input filename
#input_AOV.FN <- "02_input_data/clam.csv"      # the original input filename
# Input sediment data (#TODO: should be replaced) [ THIS WAS THE INPUT PROVIDED WITH THE SCRIPT ]
#cgsediment <- read.csv("02_input_data/cgsediment.csv")
#cgsediment <- read.csv("02_input_data/cg_sediment_data_2022-03-25.csv")
# /END/ TODO: delete these lines, as we now use one single input file #

# Load data
sed_pheno.df <- read.csv(file = input.FN)
head(sed_pheno.df)
str(sed_pheno.df)

# Set variables as character if not including in the PCA
sed_pheno.df$plot <- as.character(sed_pheno.df$plot)
sed_pheno.df$day <- as.character(sed_pheno.df$day)
str(sed_pheno.df)

row.names(sed_pheno.df) <- sed_pheno.df[,"beach_id"] # provide row names
head(sed_pheno.df)

# Add grouping vector
sed_pheno.df$Group <- "NA"
sed_pheno.df[grep(pattern = "A|B", x = sed_pheno.df$beach), "Group"] <- "A"
sed_pheno.df[grep(pattern = "C|D", x = sed_pheno.df$beach), "Group"] <- "B"
sed_pheno.df[grep(pattern = "E|F", x = sed_pheno.df$beach), "Group"] <- "C"


#### 01. PCA based on physical characteristics ####
# PCA based on growth and physical attributes using prcomp
pca.df <- select_if(sed_pheno.df, is.numeric) # only selects numeric columns

pca_fit <- pca.df %>% 
  scale() %>%                    # scales the numeric columns
  prcomp()

summary(pca_fit)

# Plot biplot
str(sed_pheno.df$Group) 

# Plot PCA (ggplot)
pdf(file = "03_pheno_results/per_plot_abiotic_PCA.pdf", width = 6, height = 6)
par(mfrow = c(1,1))
p <- ggbiplot(pcobj = pca_fit, labels = row.names(sed_pheno.df)
              , color = 'black', varname.adjust = 1.1, varname.size = 2.9
              , groups = sed_pheno.df$Group, ellipse = TRUE
              ) + expand_limits(x = c(-2, 2), y = c(-2.5,2))

p # uncomment if want to keep grey grid

# optional to remove grey grid
# p <- p + theme(panel.background = element_rect(fill = "white", colour = NA) )
# p <- p + theme(panel.grid = element_blank(), legend.position = "none")
# p
dev.off()


#### 02. Statistical analysis, linear models ####
# Input data (TODO: should use the same file for both PCA and AOV)
clam <- read.csv(file = input_AOV.FN)
head(clam) # this eventually will be replaced by the same input file as above
#str(clam)

# Remove unneeded cols
drop_cols <- c("cage", "beach_label", "orig")
clam <- clam[, !(colnames(clam) %in% drop_cols)]
rm(drop_cols) # clean enviro
head(clam) # this eventually will be replaced by the same input file as above


#### 02.1 Effect of abiotic variable on growth and survival ####
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
  abiotic_fx.list[[paste0("surv_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]])
  
  # Linear model per explan variable, growth
  abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]] <- lm(clam[, "grow"] ~ clam[, voi])
  abiotic_fx.list[[paste0("grow_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]])
  
}

# Then write out the output
capture.output(abiotic_fx.list, file = "03_pheno_results/abiotic_variables_on_grow_surv_models.txt")


#### 02.2 Effect of clam garden status on sediment variables ####
# Fit a linear mixed effects model 
# Fixed effect: beach type; Random factor: site, with plot nested in site

# Required packages
# install.packages("nlme")
# install.packages("multcompView")
# install.packages("lsmeans")
# install.packages("TukeyC")
# library("multcompView")
# library("lsmeans")
# library("TukeyC")

head(clam)

# Nested ANOVA with mixed effects model (nlme)
clam$Type <- as.factor(clam$Type)
clam$beach <- as.factor(clam$beach)

# THIS APPEARS TO BE THE NESTED WITH RANDOM VARIABLE, CAME COMMENTED OUT
#model = lme(surv ~ org, random = ~ 1|beach, data = clam)
#summary(model)
#anova.lme(model, type = "sequential", adjustSigma = FALSE)

# Plot abiotic variable levels by beach (exploratory only)
par(mfrow = c(2,2))
boxplot(clam$surv ~ clam$beach)
boxplot(clam$sand ~ clam$beach) # note descending effect on sand
boxplot(clam$carb ~ clam$beach) 
boxplot(clam$inwt ~ clam$beach) 

#### Complex statistical analyses starts, with only one variable shown ####
head(sed_pheno.df)
head(cgsediment)

# Linear mixed-effects model, with nested random effect
cgsediment$Type <- as.factor(cgsediment$Type) # Make CG/ Ref into factor
model <- lme(silt ~ Type, random = ~ 1|beach, data = cgsediment) # example with silt
summary(model)
anova.lme(model)

# Linear model without nested random effect (for reference only)
mod.no.nest <- lm(silt ~ Type, data = cgsediment)
summary(mod.no.nest)

# # Analysis of random effect (beach) - fit linear model using generalized least squares
# model.fixed <- gls(silt ~ Type, data = cgsediment)
# anova(model, model.fixed)

# # Post-hoc comparison of least-square means (general linear hypotheses)
# posthoc <- glht(model, linfct = mcp(Type = "Tukey"))
# mcs <- summary(posthoc,test = adjusted("single-step"))
# mcs

# # Set up a compact letter display of all pair-wise comparisons
# cld(mcs, level = 0.05, decreasing = TRUE)

# Checking assumptions of model
hist(residuals(model))
rug(residuals(model))
plot(fitted(model), residuals(model))
plot(model)

# # Mixed effects model with lmer
# library(lmerTest) # Keep this here, to overload lmer from lme4
# cgsediment$beach <- as.factor(cgsediment$beach)
# model <- lmer(silt ~ Type + (1|beach), data = cgsediment, REML = TRUE)
# anova(model)

# # Nested ANOVA with the aov function
# fit <- aov(silt ~ Type + Error(beach), data = cgsediment)
# summary(fit)

# # Using Means Sq and Df values to get p-value for H = Type and Error = beach
# pf(q= 5.556/6.056, df1=1, df2=4, lower.tail = FALSE) # the F distribution (probability distribution function)
# # TODO: discuss with Monique

# Shapiro-Wilk test for normality
shapiro.test(residuals(model)) # p < 0.05 is a significant deviation from normality

# Bartlett's test for homogeneity of variance
bartlett.test(silt ~ interaction(Type, beach), data = cgsediment) # p < 0.05 is significant deviation from homogeneity

# par(mfrow = c(2,2))
# plot(cgsediment$Type, cgsediment$silt)
# plot(cgsediment$beach, cgsediment$silt)
# plot(cgsediment$surv, cgsediment$silt)

#### TODO: automate all variables, using only the required functions ####
# /END section / #


#### 03. Correlation of variables ####
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

#### 04. Extra plotting and exploration ####
# Survival and growth by beach
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

#### 05. Effect of in-weight exploration ####
# Is there a trend of higher inwt and higher survival? 
plot(x = clam$inwt, y = clam$surv) 
mod.surv.inwt <- lm(clam$surv ~ clam$inwt) # response by terms
summary(mod.surv.inwt) # yes, p = 0.03

# Is carbonate also correlated with inwt? 
plot(x = clam$inwt, y = clam$carb) # higher inwt, lower carbonate
mod1 <- lm(clam$carb ~ clam$inwt)
summary(mod1)          # yes, p = 0.009

# Is survival associated with carbonate? 
plot(x = clam$surv, y = clam$carb)
mod.surv.carb <- lm(clam$surv ~ clam$carb)
summary(mod.surv.carb) # yes, p = 0.0011

## Why does this remove the effect of inwt? Not clear, do not use 
# mod.surv.carb.inwt <- lm(clam$surv ~ clam$carb + clam$inwt)
# summary(mod.surv.carb.inwt)

# Is growth associated with inwt? 
plot(x = clam$inwt, y = clam$grow)
mod.grow.by.inwt <- lm(clam$grow ~ clam$inwt)
summary(mod.grow.by.inwt) # no, p = 0.37

# Is growth associated with carbonate? 
mod.grow.by.carb <- lm(clam$grow ~ clam$carb)
summary(mod.grow.by.carb) # yes, p = 0.017

# Plotting these variables to view these trends
par(mfrow=c(2,3))
plot(x = clam$carb, y = clam$surv)
plot(x = clam$inwt, y = clam$surv)
plot(x = clam$inwt, y = clam$carb)
plot(x = clam$carb, y = clam$grow)
plot(x = clam$inwt, y = clam$grow)
plot(x = clam$grow, y = clam$surv)

# Is growth correlated to survival? 
mod.grow.by.surv <- lm(clam$grow ~ clam$surv)
summary(mod.grow.by.surv) # yes, p = 0.02


## Further test: split beaches by weight group ##
# Observe again the differences in inwts between the beaches
par(mfrow = c(1,1))
boxplot(clam$inwt ~ clam$beach) 
#boxplot(clam$inwt ~ clam$beach)  ####  CHECK WITH INITIAL HEIGHT
# Beaches A, C, E, F are the lowest in weight beaches

# IS THERE A SIGNIFICANT DIFFERENCE IN INWT IN BEACH A, C, E, F

# Is there a significant effect of carbonate on survivorship in low in-weight beaches? 
par(mfrow = c(2,2))
plot(x = clam$carb[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
      , y = clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
      , xlab = "carb (low weight beaches)"
       , ylab = "survivorship")

mod <- lm(clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"] ~ 
                      clam$carb[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
          )
                    
summary(mod) # p = 0.026

# Is there a significant effect of carbonate on survivorship in high in-weight beaches? 
plot(x = clam$carb[clam$beach=="B" | clam$beach=="D"]
     , y = clam$surv[clam$beach=="B" | clam$beach=="D"]
     , xlab = "carb (high weight beaches)"
     , ylab = "survivorship")

mod <- lm(clam$surv[clam$beach=="B" | clam$beach=="D"] ~ clam$carb[clam$beach=="B" | clam$beach=="D"])
summary(mod) # no, p = 0.2 (but note the lower n)

# Then the complementary analysis: 
# Is there an effect of in.wt on survivorship in low in-weight beaches?  
plot(x = clam$inwt[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
     , y = clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
     , xlab = "inwt (low weight beaches)"
     , ylab = "survivorship")

mod <- lm(clam$surv[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"] ~ 
            clam$inwt[clam$beach=="C" | clam$beach=="E" | clam$beach=="F" | clam$beach=="A"]
)
summary(mod) # no, p = 0.5796

# Is there an effect of in.wt on survivorship in high weight beaches? 
plot(x = clam$inwt[clam$beach=="B" | clam$beach=="D"]
     , y = clam$surv[clam$beach=="B" | clam$beach=="D"]
     , xlab = "inwt (high weight beaches)"
     , ylab = "survivorship")

mod <- lm(clam$surv[clam$beach=="B" | clam$beach=="D"] ~ clam$inwt[clam$beach=="B" | clam$beach=="D"])
summary(mod) # yes,  p = 0.004 (but be aware of low n))

# In conclusion: 
# When separating beaches into low or high inwt beaches, 
# we see a significant effect of carbonate on survivorship in the low inwt beaches, but not in the high inwt
# we see a significant effect of inwt on survivorship in the high inwt beaches, but not in the low inwt beaches


#### END ####

