#Barplot of Functional clusters between High/Medium/Low survival groups
# Present clear visual for paper
# Monique_May 17/2021
# R version 4.0.4
# R packages: 

getwd()
#
rm(list = ls())
ls()
#
library(rio)
#
Funclust <- import("cgHML.csv")
#
str(Funclust)
#
library(ggplot2)
#
Funclust$Comparison <- factor(Funclust$Comparison, levels = c("High vs Low survival", "High vs Medium survival", "Low vs High survival"))
Funclust$Tissue <- factor(Funclust$Tissue, levels = c("Gill", "DG"))
Funclust$Term <- as.character(Funclust$Term)
Funclust$Genes <- as.integer(Funclust$Genes)
#
str(Funclust)
#
#
#plot(Funclust$Term, Funclust$Genes)

?stat_count

p <- ggplot(data = Funclust, aes(x = Term,y = Genes, fill = Tissue)) + geom_col( )+ scale_fill_grey() 
p <- p + theme(panel.background = element_rect(fill = "white", colour = "black"))
p <- p + theme(panel.grid = element_blank())
p <- p + facet_wrap(~Comparison)+ coord_flip() 
p

p <- ggplot(data = Funclust, aes(x = Term,y = Genes, fill = Tissue)) + geom_col() 
p <- p + facet_wrap(~Comparison)+ theme_set(theme_bw())+ coord_flip()
p


?facet_wrap
p1 + theme(panel.background = element_rect(fill = "white", colour = "black"))


