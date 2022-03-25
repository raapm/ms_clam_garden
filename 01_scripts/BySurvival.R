
# To identify sediment characteristics that may have a positive impact on clam survival and growth
# 
# by Monique Raap
# Started March 22, 2021
# R version 4.0.4
# R packageS: rio, ggplot2, ggeasy, patchwork. 
############################################

library(nlme)
library(lme4)


getwd()

Sediment <- import("Data\\cgsediment.csv")