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
# Set input filename
pheno.FN <- "02_input_data/cg_sediment_phenos_2022-06-24.csv"

# Load data
sed_pheno.df <- read.csv(file = pheno.FN)
head(sed_pheno.df)

# User select variable #
# Growth - only use one measurement, either height based on weight or height, as selected below (choose the one to drop)
growth_to_drop <- "grow.wt.surv" # Keep height
#growth_to_drop <- "grow.ht.all" # Keep weight

# Remove the unselected growth estimate
print(paste0("You have chosen to drop ", growth_to_drop, " and therefore keep the alternate growth measure"))
print("Dropping now")
sed_pheno.df <- sed_pheno.df[, -(which(colnames(sed_pheno.df)==growth_to_drop))] # remove the unselected
head(sed_pheno.df)
colnames(sed_pheno.df)[grep(pattern = "grow", x = colnames(sed_pheno.df))] <- "grow" # rename the retained
head(sed_pheno.df)

# Set variables that are not to be included in the PCA as characters
sed_pheno.df$plot <- as.character(sed_pheno.df$plot)
sed_pheno.df$day <- as.character(sed_pheno.df$day)
str(sed_pheno.df)

row.names(sed_pheno.df) <- sed_pheno.df[,"beach_id"] # provide row names
head(sed_pheno.df)

# Add grouping vector (manually determined)
sed_pheno.df$Group <- "NA"
sed_pheno.df[grep(pattern = "A|B", x = sed_pheno.df$beach), "Group"] <- "A"
sed_pheno.df[grep(pattern = "C|D", x = sed_pheno.df$beach), "Group"] <- "B"
sed_pheno.df[grep(pattern = "E|F", x = sed_pheno.df$beach), "Group"] <- "C"

# # Attempts with data transformations (exploratory)
# bartlett.test(sed_pheno.df$inwt ~ interaction(Type, beach), data = sed_pheno.df)
# bartlett.test(log2(sed_pheno.df$inwt) ~ interaction(Type, beach), data = sed_pheno.df)
# bartlett.test(log10(sed_pheno.df$inwt) ~ interaction(Type, beach), data = sed_pheno.df)
# boxplot(sed_pheno.df$inwt ~ sed_pheno.df$Type * sed_pheno.df$beach)
# boxplot(log2(sed_pheno.df$inwt) ~ sed_pheno.df$Type * sed_pheno.df$beach)
# bartlett.test(sed_pheno.df$inwt ~ Type, data = sed_pheno.df)
# boxplot(sed_pheno.df$inwt ~ sed_pheno.df$Type)


#### 01. PCA based on physical characteristics ####
# PCA based on growth and physical attributes using prcomp
pca.df <- select_if(sed_pheno.df, is.numeric) # only selects numeric columns

colnames(pca.df)

# Rename column headers for figure
colnames(pca.df)[colnames(pca.df)=="inwt"] <- "initial weight"
colnames(pca.df)[colnames(pca.df)=="inht"] <- "initial height"
colnames(pca.df)[colnames(pca.df)=="grow"] <- "growth"
colnames(pca.df)[colnames(pca.df)=="surv"] <- "survival"
colnames(pca.df)[colnames(pca.df)=="carb"] <- "carbonate"
colnames(pca.df)[colnames(pca.df)=="org"] <- "organics"
colnames(pca.df)[colnames(pca.df)=="srocks"] <- "small rocks"
colnames(pca.df)[colnames(pca.df)=="vcsand"] <- "very course sand"
colnames(pca.df)[colnames(pca.df)=="csand"] <- "course sand"
colnames(pca.df)[colnames(pca.df)=="fsand"] <- "fine sand"
colnames(pca.df)[colnames(pca.df)=="vfsand"] <- "very fine sand"

# Scale and run PCA
pca_fit <- pca.df %>% 
  scale() %>%                    # scales the numeric columns
  prcomp()

summary(pca_fit)

# Plot PCA (ggplot)
pdf(file = "03_pheno_results/per_plot_abiotic_PCA.pdf", width = 6, height = 6)
par(mfrow = c(1,1))
p <- ggbiplot(pcobj = pca_fit, labels = row.names(sed_pheno.df)
              , color = 'black', varname.adjust = 1.3, varname.size = 2.9
              , groups = sed_pheno.df$Group, ellipse = TRUE
              ) + expand_limits(x = c(-2, 2), y = c(-2.5,2))

summary_details <- summary(pca_fit)

custom_x_axis  <-  paste0("Standardized PC1 ("
                       , round(
                            summary_details[["importance"]][2,1] * 100
                            , digits = 1
                            )
                       , "%)"
                      )

custom_y_axis  <-  paste0("Standardized PC2 ("
                          , round(
                            summary_details[["importance"]][2,2] * 100
                            , digits = 1
                          )
                          , "%)"
)

p + xlab(custom_x_axis) + ylab(custom_y_axis)
  

#p # uncomment for default xlab/ylab, if want to keep grey grid

# optional to remove grey grid
# p <- p + theme(panel.background = element_rect(fill = "white", colour = NA) )
# p <- p + theme(panel.grid = element_blank(), legend.position = "none")
# p
dev.off()

# For confirmation of custom values, should look the same
pdf(file = "03_pheno_results/per_plot_abiotic_PCA_default_axes.pdf", width = 6, height = 6)

p <- ggbiplot(pcobj = pca_fit, labels = row.names(sed_pheno.df)
              , color = 'black', varname.adjust = 1.3, varname.size = 2.9
              , groups = sed_pheno.df$Group, ellipse = TRUE
              ) + expand_limits(x = c(-2, 2), y = c(-2.5,2))

p
dev.off()

#### 02. Statistical analysis, linear models ####
head(sed_pheno.df)
str(sed_pheno.df)

# Remove unneeded cols
drop_cols <- c("plot", "beach_id", "orig", "Group")
sed_pheno.df <- sed_pheno.df[, !(colnames(sed_pheno.df) %in% drop_cols)]
rm(drop_cols) # clean enviro
head(sed_pheno.df)

#### 02.0 Evaluate inwt and inht per beach

pdf(file = "03_pheno_results/input_weight_height_by_beach.pdf", width = 7, height = 5)
par(mfrow=c(1,2))

# Weight
boxplot(sed_pheno.df$inwt ~ sed_pheno.df$beach
        , las = 1
        , xlab = "Beach"
        , ylab = "input weight (g)"
        )

# Height
boxplot(sed_pheno.df$inht ~ sed_pheno.df$beach
        , las = 1
        , xlab = "Beach"
        , ylab = "input height (cm)"
)

dev.off()

# Significance testing
#Weight
mod.wt <- aov(formula = sed_pheno.df$inwt ~ sed_pheno.df$beach)
summary(mod.wt)
TukeyHSD(mod.wt)

#Height
mod.ht <- aov(formula = sed_pheno.df$inht ~ sed_pheno.df$beach)
summary(mod.ht)
TukeyHSD(mod.ht)


#### 02.1 Effect of abiotic variables on growth and survival ####
# Linear models of all variables on survival and growth

# Which columns are to be considered response variables? 
response_variables <- c(
                     # "beach", 
                      "grow", "surv"
                      )

explan_vars <- colnames(
                         sed_pheno.df[, !(colnames(sed_pheno.df) %in% response_variables)]
                        )
rm(response_variables) # clean enviro

print(paste0("These are your explanatory variables, tested based on their effect on growth and survival:"))
print(explan_vars)

# Set nulls, Loop to run lm per explan variable
abiotic_fx.list <- list(); voi <- NULL

for(i in 1:length(explan_vars)){
  
  voi <- explan_vars[i]
  
  # Linear model per explan variable, survival
  abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]] <- lm(sed_pheno.df[, "surv"] ~ sed_pheno.df[, voi])
  abiotic_fx.list[[paste0("surv_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]])
  
  # Linear model per explan variable, growth
  abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]] <- lm(sed_pheno.df[, "grow"] ~ sed_pheno.df[, voi])
  abiotic_fx.list[[paste0("grow_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]])
  
}

# Then write out the output
capture.output(abiotic_fx.list, file = "03_pheno_results/abiotic_variables_on_grow_surv_models.txt")

# Extract the important values for the table
explan_vars

# Select only continuous variables
cont_explan_vars <- explan_vars[grep(pattern = "Type|beach|day", x = explan_vars, invert = T)]

data.mat <- matrix(nrow = length(cont_explan_vars), ncol = 5)
colnames(data.mat) <- c("explan_var", "surv_pval", "surv_rsq", "grow_pval", "grow_rsq")
data.df <- as.data.frame(data.mat)
data.df

for(i in 1:length(cont_explan_vars)){
  
  voi <- cont_explan_vars[i]
  
  data.df[i ,"explan_var"] <- voi
  
  # Survival
  out.obj <- summary(abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]])
  rsq <- out.obj$adj.r.squared
  #pval <- out.obj$coefficients["sed_pheno.df[, voi]" ,"Pr(>|t|)"] # does not work
  pval <- out.obj$coefficients[
                                grep(pattern = "sed_pheno.df", x = rownames(out.obj$coefficients))
                                ,"Pr(>|t|)"
                                ]
                              
  data.df[i, "surv_pval"] <- pval
  data.df[i, "surv_rsq"] <- rsq
  
  # Growth
  out.obj <- summary(abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]])
  rsq <- out.obj$adj.r.squared
  #pval <- out.obj$coefficients["sed_pheno.df[, voi]" ,"Pr(>|t|)"] # does not work
  pval <- out.obj$coefficients[
    grep(pattern = "sed_pheno.df", x = rownames(out.obj$coefficients))
    ,"Pr(>|t|)"
  ]
  
  data.df[i, "grow_pval"] <- pval
  data.df[i, "grow_rsq"] <- rsq
  
  
}

data.df

# Save out results
write.csv(x = data.df, file = "03_pheno_results/abiotic_variables_on_grow_surv_table.csv", quote = F)

# Second run, without Beach D (the significantly bigger at planting beach)
sed_pheno.df.bck <- sed_pheno.df # backup the original
sed_pheno.df <- sed_pheno.df[sed_pheno.df$beach!="D", ]

# Re-run, as above
# Set nulls, Loop to run lm per explan variable
abiotic_fx.list <- list(); voi <- NULL

for(i in 1:length(explan_vars)){
  
  voi <- explan_vars[i]
  
  # Linear model per explan variable, survival
  abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]] <- lm(sed_pheno.df[, "surv"] ~ sed_pheno.df[, voi])
  abiotic_fx.list[[paste0("surv_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]])
  
  # Linear model per explan variable, growth
  abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]] <- lm(sed_pheno.df[, "grow"] ~ sed_pheno.df[, voi])
  abiotic_fx.list[[paste0("grow_by_", voi, ".summary")]] <- summary(abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]])
  
}

# Then write out the output
capture.output(abiotic_fx.list, file = "03_pheno_results/abiotic_variables_on_grow_surv_models_no_beach_D.txt")

# Extract the important values for the table
explan_vars

# Select only continuous variables
cont_explan_vars <- explan_vars[grep(pattern = "Type|beach|day", x = explan_vars, invert = T)]

data.mat <- matrix(nrow = length(cont_explan_vars), ncol = 5)
colnames(data.mat) <- c("explan_var", "surv_pval", "surv_rsq", "grow_pval", "grow_rsq")
data.df <- as.data.frame(data.mat)
data.df

for(i in 1:length(cont_explan_vars)){
  
  voi <- cont_explan_vars[i]
  
  data.df[i ,"explan_var"] <- voi
  
  # Survival
  out.obj <- summary(abiotic_fx.list[[paste0("surv_by_", voi, ".mod")]])
  rsq <- out.obj$adj.r.squared
  #pval <- out.obj$coefficients["sed_pheno.df[, voi]" ,"Pr(>|t|)"] # does not work
  pval <- out.obj$coefficients[
    grep(pattern = "sed_pheno.df", x = rownames(out.obj$coefficients))
    ,"Pr(>|t|)"
  ]
  
  data.df[i, "surv_pval"] <- pval
  data.df[i, "surv_rsq"] <- rsq
  
  # Growth
  out.obj <- summary(abiotic_fx.list[[paste0("grow_by_", voi, ".mod")]])
  rsq <- out.obj$adj.r.squared
  #pval <- out.obj$coefficients["sed_pheno.df[, voi]" ,"Pr(>|t|)"] # does not work
  pval <- out.obj$coefficients[
    grep(pattern = "sed_pheno.df", x = rownames(out.obj$coefficients))
    ,"Pr(>|t|)"
  ]
  
  data.df[i, "grow_pval"] <- pval
  data.df[i, "grow_rsq"] <- rsq
  
  
}

data.df

write.csv(x = data.df, file = "03_pheno_results/abiotic_variables_on_grow_surv_table_no_beach_D.csv", quote = F)

# Clean space 
rm(explan_vars)
sed_pheno.df <- sed_pheno.df.bck # Revert back to the full dataset

#### 02.2 Effect of clam garden status on sediment variables ####
# Fit a linear mixed effects model 
## where fixed effect is beach type, and random factor is beach (plot nested in beach)
head(sed_pheno.df)
str(sed_pheno.df)

# Update character to factor
sed_pheno.df$Type <- as.factor(sed_pheno.df$Type)     # Make CG/ Ref into factor
sed_pheno.df$beach <- as.factor(sed_pheno.df$beach)   # Make beach into factor

# Nested ANOVA with mixed effects model (nlme)
# Define variables that are NOT response variables
non_resp_vars <- c("Type", "beach", "day")

explan_vars <- colnames(sed_pheno.df)

resp_vars <- explan_vars[grep(pattern = paste0(non_resp_vars, collapse = "|")
                            , x = explan_vars
                            , invert = T
                            )]
rm(non_resp_vars)
print(paste0("The following are response variables to be considered: "))
print(resp_vars)

# Set nulls
cg_fx.list <- list(); voi <- NULL; temp.df <- NULL; mod <- NULL

for(i in 1:length(resp_vars)){
  
  voi <- resp_vars[i]
  
  # Need to separate out the section of interest, including only the voi, the type, and beach (due to potential limitations of lme syntax)
  temp.df <- sed_pheno.df[,c(voi, "Type", "beach")]
  colnames(x = temp.df)[1] <- "select_voi"
  
  # Generate summary statistics
  print(paste0("***Summary statistics for ", voi, "***"))
  cg_fx.list[[paste0(voi, "_summary_statistics")]] <- tapply(temp.df$select_voi, temp.df$Type, summary)
  cg_fx.list[[paste0(voi, "_sd")]] <- tapply(temp.df$select_voi, temp.df$Type, sd)
  
  # Linear mixed-effects model, with nested random effect; response variable as a function of Type
  mod <- lme(select_voi ~ Type, random = ~ 1|beach, data = temp.df)
  
  # Store results in list
  cg_fx.list[[paste0(voi, "_by_CG_type.mod")]] <- mod
  cg_fx.list[[paste0(voi, "_by_CG_type.summary")]] <- summary(mod)
  cg_fx.list[[paste0(voi, "_by_CG_type.anova.lme.summary")]] <- anova.lme(mod) # redundant
  cg_fx.list
  
  ## Assumption tests
  # Shapiro-Wilk test for normality
  cg_fx.list[[paste0(voi, "_by_CG_type_mod_residuals_shapiro")]] <- shapiro.test(residuals(mod))
  # p < 0.05 is a significant deviation from normality
  
  # Bartlett's test for homogeneity of variance
  cg_fx.list[[paste0(voi, "_by_CG_type_mod_bartlett")]] <- bartlett.test(select_voi ~ interaction(Type, beach), data = temp.df)
  # p < 0.05 is significant deviation from homogeneity
  
  # histogram of residuals
  pdf(file = paste0("03_pheno_results/model_assumption_graphs/hist_residuals_nested_mod_", voi, ".pdf"), width = 5, height = 5)
  hist(residuals(mod), main = voi, xlab = "Residuals of nested model")
  rug(residuals(mod))
  dev.off()
  
  # residuals by fitted plot
  pdf(file = paste0("03_pheno_results/model_assumption_graphs/standard_residuals_by_fitted_val_", voi, ".pdf"), width = 5, height = 5)
  print(plot(mod)) # print is required here, or else a plot will not be produced
  dev.off()
  
  # Plot the data
  pdf(file = paste0("03_pheno_results/resp_var_by_CG_type_and_beach_", voi, ".pdf")
      , width = 9.5, height = 6)
  par(mfrow = c(2,2))
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type, ylab = voi, xlab = "Type")
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$beach, ylab = voi, xlab = "beach")
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type *  sed_pheno.df$beach
          , ylab = voi, xlab = "Type x beach")
  dev.off()
}

# Then write out the output
capture.output(cg_fx.list, file = "03_pheno_results/clam_garden_on_response_variables_models_nested.txt")


# Also collect data without nesting for reference only
# Set nulls
cg_fx_no_nest.list <- list(); voi <- NULL; temp.df <- NULL

for(i in 1:length(resp_vars)){
  
  voi <- resp_vars[i]
  
  # Need to separate out the section of interest, including only the voi, the type, and beach (due to potential limitations of lme syntax)
  temp.df <- sed_pheno.df[,c(voi, "Type", "beach")]
  colnames(x = temp.df)[1] <- "select_voi"
  
  # Linear model, without nesting; response variable as a function of Type
  mod <- lm(select_voi ~ Type, data = temp.df)
  
  # Store results in list
  cg_fx_no_nest.list[[paste0(voi, "_by_CG_type.mod")]] <- mod
  cg_fx_no_nest.list[[paste0(voi, "_by_CG_type.summary")]] <- summary(mod)
  cg_fx_no_nest.list[[paste0(voi, "_by_CG_type.anova.summary")]] <- anova(mod)
  cg_fx_no_nest.list
  
  ## Assumption tests
  # Shapiro-Wilk test for normality
  cg_fx_no_nest.list[[paste0(voi, "_by_CG_type_mod_residuals_shapiro")]] <- shapiro.test(residuals(mod))
  # p < 0.05 is a significant deviation from normality
  
  # Bartlett's test for homogeneity of variance
  cg_fx_no_nest.list[[paste0(voi, "_by_CG_type_mod_bartlett")]] <- bartlett.test(select_voi ~ Type, data = temp.df)
  # p < 0.05 is significant deviation from homogeneity
  
  # histogram of residuals
  pdf(file = paste0("03_pheno_results/model_assumption_graphs/hist_residuals_not_nested_mod_", voi, ".pdf"), width = 5, height = 5)
  hist(residuals(mod), main = voi, xlab = "Residuals of nested model")
  rug(residuals(mod))
  dev.off()
  
  # residuals by fitted plot
  pdf(file = paste0("03_pheno_results/model_assumption_graphs/standard_residuals_not_nested_by_fitted_val_", voi, ".pdf"), width = 5, height = 5)
  print(plot(mod))
  dev.off()

  # Plot the data
  pdf(file = paste0("03_pheno_results/resp_var_by_CG_type_and_beach_not_nested_", voi, ".pdf")
      , width = 9.5, height = 6)
  par(mfrow = c(2,2))
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type, ylab = voi, xlab = "Type")
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$beach, ylab = voi, xlab = "beach")
  boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type *  sed_pheno.df$beach
          , ylab = voi, xlab = "Type x beach")
  dev.off()
  
  
}

# Then write out the output
capture.output(cg_fx_no_nest.list, file = "03_pheno_results/clam_garden_on_response_variables_models_not_nested.txt")


#### 02.3 Effect of beach location on all variables ####
# Identify numeric cols
numeric_vars <- colnames(select_if(sed_pheno.df, is.numeric)) # only selects numeric columns

sed_pheno.df

# Also collect data without nesting for reference only
# Set nulls
beach_fx.list <- list(); voi <- NULL; temp.df <- NULL

for(i in 1:length(numeric_vars)){
  
  voi <- numeric_vars[i]
  
  # Need to separate out the section of interest, including only the voi, the type, and beach (due to potential limitations of lme syntax)
  temp.df <- sed_pheno.df[,c(voi, "beach")]
  colnames(x = temp.df)[1] <- "select_voi"
  
  # One-way ANOVA
  mod <- lm(select_voi ~ beach, data = temp.df)
  
  # Store results in list
  beach_fx.list[[paste0(voi, "_by_beach.mod")]] <- mod
  beach_fx.list[[paste0(voi, "_by_beach.summary")]] <- summary(mod)
  beach_fx.list[[paste0(voi, "_by_beach.anova.summary")]] <- anova(mod)
  beach_fx.list
  
  # Generate summary statistics
  print(paste0("***Summary statistics for ", voi, "***"))
  beach_fx.list[[paste0(voi, "_means")]] <- tapply(temp.df$select_voi, temp.df$beach, mean)
  beach_fx.list[[paste0(voi, "_sd")]] <- tapply(temp.df$select_voi, temp.df$beach, sd)
  
  ## Assumption tests
  # Shapiro-Wilk test for normality
  beach_fx.list[[paste0(voi, "_by_beach_mod_residuals_shapiro")]] <- shapiro.test(residuals(mod))
  # p < 0.05 is a significant deviation from normality
  
  # Bartlett's test for homogeneity of variance
  beach_fx.list[[paste0(voi, "_by_beach_mod_bartlett")]] <- bartlett.test(select_voi ~ interaction(beach), data = temp.df)
  # p < 0.05 is significant deviation from homogeneity
  
  # # histogram of residuals
  # pdf(file = paste0("03_pheno_results/model_assumption_graphs/hist_residuals_not_nested_mod_", voi, ".pdf"), width = 5, height = 5)
  # hist(residuals(mod), main = voi, xlab = "Residuals of nested model")
  # rug(residuals(mod))
  # dev.off()
  # 
  # # residuals by fitted plot
  # pdf(file = paste0("03_pheno_results/model_assumption_graphs/standard_residuals_not_nested_by_fitted_val_", voi, ".pdf"), width = 5, height = 5)
  # print(plot(mod))
  # dev.off()
  # 
  # # Plot the data
  # pdf(file = paste0("03_pheno_results/resp_var_by_CG_type_and_beach_not_nested_", voi, ".pdf")
  #     , width = 9.5, height = 6)
  # par(mfrow = c(2,2))
  # boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type, ylab = voi, xlab = "Type")
  # boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$beach, ylab = voi, xlab = "beach")
  # boxplot(sed_pheno.df[,voi] ~ sed_pheno.df$Type *  sed_pheno.df$beach
  #         , ylab = voi, xlab = "Type x beach")
  # dev.off()
  
  
}

# Then write out the output
capture.output(beach_fx.list, file = "03_pheno_results/by_beach_numeric_vars_ANOVA.txt")


#### 03. Correlation of variables ####
colnames(sed_pheno.df)
to_cor_phenos.vec <- c("grow", "surv", "inwt", "inht", "carb", "org", "rocks", "srocks", "vcsand", "csand", "sand"
                       , "fsand", "vfsand", "silt")

# Exclude outlier beach? Toggle below
exclude_beach <- FALSE
#exclude_beach <- TRUE

if(exclude_beach==TRUE){
  
  sed_pheno_for_corr.df <- sed_pheno.df[sed_pheno.df$beach!="D", ]
  special_tag <- "_excl_beach_d"
  
}else{

  sed_pheno_for_corr.df <- sed_pheno.df  
  special_tag <- ""
  
}


# Correlate
cor.set <- cor(sed_pheno_for_corr.df[, to_cor_phenos.vec]
               , use = "pairwise.complete.obs" # cor bw ea pair var is done with all pairs of obs on those var 
               # note: only works with "pearson" method
               )

# Rename variables
colnames(cor.set) <- gsub(pattern = "grow", replacement = "growth", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "grow", replacement = "growth", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "surv", replacement = "survival", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "surv", replacement = "survival", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "inwt", replacement = "input weight ", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "inwt", replacement = "input weight", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "inht", replacement = "input height ", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "inht", replacement = "input height", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "carb", replacement = "carbonate", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "carb", replacement = "carbonate", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "org", replacement = "organic", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "org", replacement = "organic", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "srocks", replacement = "small rocks", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "srocks", replacement = "small rocks", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "vcsand", replacement = "very course sand", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "vcsand", replacement = "very course sand", x = rownames(cor.set))

colnames(cor.set) <- gsub(pattern = "csand", replacement = "course sand", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "csand", replacement = "course sand", x = rownames(cor.set))

# Must be done before 
colnames(cor.set) <- gsub(pattern = "vfsand", replacement = "very fine sand", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "vfsand", replacement = "very fine sand", x = rownames(cor.set))

# Must be done after
colnames(cor.set) <- gsub(pattern = "fsand", replacement = "fine sand", x = colnames(cor.set))
rownames(cor.set) <- gsub(pattern = "fsand", replacement = "fine sand", x = rownames(cor.set))


cor.set

# Plot
pdf(file = paste0("03_pheno_results/surv_grow_abiotic_correlations", special_tag, ".pdf"), width = 7, height = 7)
par(mfrow = c(1,1))
corrplot(cor.set
         , method = "circle"
         #, order="hclust"
         , addshade = "all"
         #, type = "lower"
         #, outline = T
         , tl.col = "black"
         , tl.cex = 1
         ) #plot matrix
dev.off()


#### 04. Extra plotting and exploration ####
# Survival and growth by beach
pdf(file = "03_pheno_results/surv_grow_by_beach.pdf", width = 8, height = 5)
par(mfrow = c(1,2))
boxplot(sed_pheno.df$surv ~ sed_pheno.df$beach, las = 1
        , xlab = "Beach", ylab = "Survival (%)"
        , col = c("grey", "white")
        )
mtext(text = "A"
      #, text = expression(paste(bold("A")))
      , side = 3, line = 1.5
      , at = -0.5
      #, at=par("usr")[1]+0.05*diff(par("usr")[1:2])
      )
#legend("topright", legend = c("CG", "Ref"), fill = c("grey", "white"))


boxplot(sed_pheno.df$grow ~ sed_pheno.df$beach, las = 1
        , xlab = "Beach", ylab = "Growth (%)"
        , col = c("grey", "white")
)

mtext(text = "B"
      #, text = expression(paste(bold("A")))
      , side = 3, line = 1.5
      , at = -0.5
      #, at=par("usr")[1]+0.05*diff(par("usr")[1:2])
      )

legend("topright", legend = c("CG", "Ref"), fill = c("grey", "white"))

dev.off()

# Consider survival and carbonates
pdf(file = "03_pheno_results/survival_by_carbonate.pdf", width = 8, height = 5)
par(mfrow=c(1,1))
plot(x = sed_pheno.df$carb, y = sed_pheno.df$surv
     , ylab = "Survival (%)"
     , xlab = "Carbonate (%)"
)

mod <- lm(surv ~ carb, data = sed_pheno.df)
summary(mod)

results <- summary(mod)
text(x = 13, y = 85, labels = bquote("Adjusted R"^2*"="*.(round(results$adj.r.squared, digits = 2))))
text(x = 13, y = 75
     , labels = bquote(italic(p)*"-value: "*.(round(results$coefficients["carb","Pr(>|t|)"], digits = 5))))

dev.off()


#### 05. Effect of in-weight exploration ####
# Is there a trend of higher inwt and higher survival? 
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$surv) 
mod.surv.inwt <- lm(sed_pheno.df$surv ~ sed_pheno.df$inwt) # response by terms
summary(mod.surv.inwt) # yes, p = 0.03

# Is carbonate also correlated with inwt? 
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$carb) # higher inwt, lower carbonate
mod1 <- lm(sed_pheno.df$carb ~ sed_pheno.df$inwt)
summary(mod1)          # yes, p = 0.005

# Is survival associated with carbonate? 
plot(x = sed_pheno.df$surv, y = sed_pheno.df$carb)
mod.surv.carb <- lm(sed_pheno.df$surv ~ sed_pheno.df$carb)
summary(mod.surv.carb) # yes, p = 0.0011

# Is growth associated with inwt? 
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$grow)
mod.grow.by.inwt <- lm(sed_pheno.df$grow ~ sed_pheno.df$inwt)
summary(mod.grow.by.inwt) 
# growth by weight: no, p = 0.37
# growth by height: no (close), p = 0.05

# Is growth associated with carbonate? 
mod.grow.by.carb <- lm(sed_pheno.df$grow ~ sed_pheno.df$carb)
summary(mod.grow.by.carb) 
# growth by weight: yes, p = 0.017
# growth by height: yes, p = 0.003

# Plotting these variables to view these trends
par(mfrow=c(2,3))
plot(x = sed_pheno.df$carb, y = sed_pheno.df$surv)
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$surv)
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$carb)
plot(x = sed_pheno.df$carb, y = sed_pheno.df$grow)
plot(x = sed_pheno.df$inwt, y = sed_pheno.df$grow)
plot(x = sed_pheno.df$grow, y = sed_pheno.df$surv)

# Is growth correlated to survival? 
mod.grow.by.surv <- lm(sed_pheno.df$grow ~ sed_pheno.df$surv)
summary(mod.grow.by.surv) 
# growth by weight: yes, p = 0.02
# growth by height: yes, p = 0.0004

## Further test: split beaches by weight group ##
# Observe again the differences in inwts between the beaches
par(mfrow = c(1,2))
boxplot(sed_pheno.df$inwt ~ sed_pheno.df$beach)
boxplot(sed_pheno.df$inht ~ sed_pheno.df$beach) 
# Beaches A, C, E, F are the lowest in weight beaches; arguably these have no significant difference between them
## and it is only beach B and D that are elevated relative to the others

# (#TODO) formally test: IS THERE A SIGNIFICANT DIFFERENCE IN INWT IN BEACH A, C, E, F? 

### Carbonate on survivorship, stratified by weight class beaches
# Is there a significant effect of carbonate on survivorship in low in-weight beaches? 
pdf(file = "03_pheno_results/weight_stratified_beaches_and_carbonate_inwt.pdf", width = 6, height = 6)
par(mfrow = c(2,2))
plot(x = sed_pheno.df$carb[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
      , y = sed_pheno.df$surv[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
      , xlab = "carb (low weight beaches)"
       , ylab = "survivorship")

mod <- lm(sed_pheno.df$surv[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"] ~ 
                      sed_pheno.df$carb[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
          )
                    
summary(mod) # yes; p = 0.026

# Is there a significant effect of carbonate on survivorship in high in-weight beaches? 
plot(x = sed_pheno.df$carb[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"]
     , y = sed_pheno.df$surv[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"]
     , xlab = "carb (high weight beaches)"
     , ylab = "survivorship")

mod <- lm(sed_pheno.df$surv[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"] ~ sed_pheno.df$carb[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"])
summary(mod) # no, p = 0.2 (but note the lower n)

### Inwt on survivorship, stratified by weight class beaches (complementary analysis)
# Is there a significant effect of in.wt on survivorship in low in-weight beaches?  
plot(x = sed_pheno.df$inwt[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
     , y = sed_pheno.df$surv[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
     , xlab = "inwt (low weight beaches)"
     , ylab = "survivorship")

mod <- lm(sed_pheno.df$surv[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"] ~ 
            sed_pheno.df$inwt[sed_pheno.df$beach=="C" | sed_pheno.df$beach=="E" | sed_pheno.df$beach=="F" | sed_pheno.df$beach=="A"]
)
summary(mod) # no, p = 0.53

# Is there an effect of in.wt on survivorship in high weight beaches? 
plot(x = sed_pheno.df$inwt[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"]
     , y = sed_pheno.df$surv[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"]
     , xlab = "inwt (high weight beaches)"
     , ylab = "survivorship")

mod <- lm(sed_pheno.df$surv[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"] ~ sed_pheno.df$inwt[sed_pheno.df$beach=="B" | sed_pheno.df$beach=="D"])
summary(mod) # yes,  p = 0.004 (be aware of low n))

dev.off()

# In conclusion for this section: 
# When separating beaches into low or high inwt beaches, 
# we see a significant effect of carbonate on survivorship in the low inwt beaches, but not in the high inwt
# we see a significant effect of inwt on survivorship in the high inwt beaches, but not in the low inwt beaches
# Conclude: the significant effect of carbonate on survivorship is independent of inwt for the equal and low inwt beaches

#### END ####

