#
rm(list= ls()) 
# Here is my data set
#
# update packages
#
#update.packages(ask=FALSE)
#
# Load data
#
beaches <- read.csv("cgsedimentPCA.csv")
#
str(beaches)
#

#beach.pca <- prcomp(beaches[,c(1:)], center = TRUE, scale. = TRUE)
#summary(beach.pca)
#
#prcomp(beaches)

#
#install.packages("tidyverse")
#install.packages("broom")

library(tidyverse)
library(broom)
#library(palmerpenguins)
#
#
row.names(beaches)<-beaches[,1]
#
pca_fit <- beaches %>%
  select(where(is.numeric)) %>%
  scale() %>%
  prcomp()
#
summary(pca_fit)
#
#
P <- pca_fit %>%
  augment(beaches) %>%
  rename_at(vars(starts_with(".fitted")),
            list(~str_replace(.,".fitted",""))) %>%
  ggplot(aes(x=PC1, 
             y=PC2, 
             ))+
  geom_point()
P <- P + theme(panel.background = element_rect(fill = "white", colour = "black"))
P <- P + theme(panel.grid = element_blank())
P
#
#
install.packages("devtools")
library(devtools)
#install_github("vqv/ggbiplot")
#
#install.packages("ggrepel")
library(ggrepel)

#
library(ggbiplot)
#
str(beaches$Group)
#
#beaches.groups <- c("A1","A2","A3","B1","B2","B3"),c("C1", "C2","C3","D1","D2","D3")
 #                   ,c("E1","E2","E3","F1","F2","F3")

p <- ggbiplot(pca_fit,labels = row.names(beaches),color = 'black', varname.adjust = 1.1, varname.size = 2.9,
              groups = beaches$Group, ellipse = TRUE)+expand_limits(x = c(-2, 2), y = c(-2.5,2))

p
p <- p + theme(panel.background = element_rect(fill = "white", colour = NA) )
p <- p + theme(panel.grid = element_blank(), legend.position = "none")
p 


