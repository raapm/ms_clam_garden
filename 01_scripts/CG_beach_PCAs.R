# Physical data and survival data
# Initialized to GitHub 2022-03-25

# clear workspace
#rm(list=ls())

# Install packages and load libraries
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



#
library(ggbiplot)

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

# Set input filenames
input.FN <- "02_input_data/cgsedimentPCA.csv" # the original input filename

# Load data
beaches <- read.csv(file = input.FN)
head(beaches)
str(beaches)

row.names(beaches) <- beaches[,"Beach"] # provide row names
head(beaches)

### May remove ###
#beach.pca <- prcomp(beaches[,c(1:)], center = TRUE, scale. = TRUE)
#summary(beach.pca)
#
#prcomp(beaches)

#library(palmerpenguins)

### /end/ May remove ###

#### 01. PCA based on physical characteristics ####
# PCA based on growth and physical attributes using prcomp
pca_fit <- beaches %>%
  select(where(is.numeric)) %>%  # only selects numeric columns
  scale() %>%                    # scales the numeric columns
  prcomp()

summary(pca_fit)

# Plot PCA
P <- pca_fit %>%
  augment(beaches) %>%                              # ?
  rename_at(vars(starts_with(".fitted")),           # renames the PCA .fitted variables as .fittedPC1 -> PC1
            list(~str_replace(.,".fitted",""))) %>% # ?
  ggplot(aes(x=PC1, 
             y=PC2, 
             ))+
  geom_point()
P <- P + theme(panel.background = element_rect(fill = "white", colour = "black"))
P <- P + theme(panel.grid = element_blank())
P

# Plot biplot
str(beaches$Group) 

### Not clear if this is required ###
#beaches.groups <- c("A1","A2","A3","B1","B2","B3"),c("C1", "C2","C3","D1","D2","D3")
 #                   ,c("E1","E2","E3","F1","F2","F3")
### /END/ Not clear if this is required ###

p <- ggbiplot(pcobj = pca_fit, labels = row.names(beaches)
              , color = 'black', varname.adjust = 1.1, varname.size = 2.9
              , groups = beaches$Group, ellipse = TRUE
              ) + expand_limits(x = c(-2, 2), y = c(-2.5,2))

p
p <- p + theme(panel.background = element_rect(fill = "white", colour = NA) )
p <- p + theme(panel.grid = element_blank(), legend.position = "none")
p 

#### Remaining questions
# Q1: what is 'X' column? It is highly correlated with growth
# Q2: what are the 'Groups' ? Is there anything that defines this other than the PCA itself? 
# Q3: can we use the same data input as the other scripts? only one data input from raw data
