# CG by beach
# To identify sediment characteristics that may have a positive impact on clam survival and growth
# 
# by Monique Raap
# Started March 22, 2021
# R version 4.0.4
# R packageS: rio, ggplot2, ggeasy, patchwork. 
############################################

library(rio)
library(ggplot2)
library(ggeasy)
library(patchwork)

getwd()

Beach <- import("Data\\clambybeach.csv")

# Looking at survival and all the sediment factors

p1 <- ggplot(data = Beach, aes( x = carb, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p2 <- ggplot(data = Beach, aes( x = org, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p3 <- ggplot(data = Beach, aes( x = rocks, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p4 <- ggplot(data = Beach, aes( x = srocks, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p5 <- ggplot(data = Beach, aes( x = vcsand, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p6 <- ggplot(data = Beach, aes( x = csand, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p7 <- ggplot(data = Beach, aes( x = sand, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p8 <- ggplot(data = Beach, aes( x = fsand, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p9 <- ggplot(data = Beach, aes( x = vfsand, y = surv, color = beach)) + geom_point()+ easy_remove_legend()

p10 <- ggplot(data = Beach, aes( x = silt, y = surv, color = beach)) + geom_point()+guides(col = guide_legend(label.hjust = 0.5))

plotAll <- p1 + p2 + p3 + p4+ p5 + p6 + p7 + p8 + p9 + p10

plotAll

ggsave(plot = plotAll, file = "Plots//CG_Surv vs sediment plots.pdf") 

# Looking at growth and all the sediment factors

pg1 <- ggplot(data = Beach, aes( x = carb, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg2 <- ggplot(data = Beach, aes( x = org, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg3 <- ggplot(data = Beach, aes( x = rocks, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg4 <- ggplot(data = Beach, aes( x = srocks, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg5 <- ggplot(data = Beach, aes( x = vcsand, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg6 <- ggplot(data = Beach, aes( x = csand, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg7 <- ggplot(data = Beach, aes( x = sand, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg8 <- ggplot(data = Beach, aes( x = fsand, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg9 <- ggplot(data = Beach, aes( x = vfsand, y = grow, color = beach)) + geom_point()+ easy_remove_legend()

pg10 <- ggplot(data = Beach, aes( x = silt, y = grow, color = beach)) + geom_point()+ guides(col = guide_legend(label.hjust = 0.5))

plotgAll <- pg1 + pg2 + pg3 + pg4+ pg5 + pg6 + pg7 + pg8 + pg9 + pg10

plotgAll

ggsave(plot = plotgAll, file = "Plots//CG_Grow vs sediment plots.pdf") 



