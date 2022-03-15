#
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.3")
#


## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

#

# clear workspace
rm(list=ls())

citation("edgeR")
library("limma")
library("edgeR")

# ******All libraries First**********************

cgrawdata <- read.csv("cgrnaseqBTv1.csv", header = TRUE)
#Create A DGEList object for easy manipulation
y <- DGEList(counts = cgrawdata[ ,2:33], genes = cgrawdata[ ,1])
#
#
# Both
colnames(cgrawdata)[1] <- "jong.id"
#
# sort by jong.id
cgrawdata.sorted <- cgrawdata[order(cgrawdata$jong.id),]
str(cgrawdata.sorted)
#
# Merge the identifiers from the complete list right in the beginning

uniprot.map <- read.table(file="project155.uniprot_blastp.txt"
                          , header = F)
colnames(uniprot.map) <- c("jong.id" , "uniprot")
str(uniprot.map)

# sort by jong.id
uniprot.map.sorted <- uniprot.map[order(uniprot.map$jong.id),]


data <- merge(x = cgrawdata.sorted, y = uniprot.map.sorted, by = "jong.id")

str(data)

# Move ID column to last position
data <- data[,c(1:35,37,36)]

str(data)


# Create A DGEList object for easy manipulation
y <- DGEList(counts = data[ ,2:33], genes = data[c(1,34,35,36,37)])
str(y)
y$genes




###### Subset samples for tissue #####
# 

y.gill <- y[,1:14]
colnames(y.gill)

y.dig <- y[,15:28]
colnames(y.dig)

#y.gCGRef <- y[,29:30]
#colnames(y.gCGRef)

y.gillCG <- y[,29]
colnames(y.gillCG)

y.gillRef <- y[,30]
colnames(y.gillRef)

y.dig.CG <- y[,31]
colnames(y.dig.CG)

y.dig.Ref <- y[,32]
colnames(y.dig.Ref)

y.all <- y[,1:28]
colnames(y.all)

#y <- y.gill
#y <- y.all
y <- y.dig 
#y <- y.gillCG
#y <- y.gillRef



#### Library Sizes ####

str(y[[2]])
y[[2]]$lib.size # lib sizes


# Google search: subsetting lists


### TMM  normalization  ####
# Normalization is  applied  to  this  dataset  to  account  for  compositional  difference
# between the libraries

y <- calcNormFactors(y, method="TMM")
y$samples




# Data Exploration:
#### BCV/PCA ####

# examine  the  samples  for  outliers  and  for  other relationships.
# The  function plotMDS produces  a  plot  in  which  distances  between  samples
# correspond to leading biological coefficient of variation (BCV) between those samples:

par(mar = c(4,4,4,4))
par(mfrow= c(1,1))

# Figure X: Multidimensional scaling of Pacific littleneck clam gill tissue libraries, logFC; gene expression fold change. The first letter (C: clam garden, or N: non-walled beach) signifies the beach type, and the second letter (A – F) signifies the beach, and the following integer (1 – 3) is the plot number. 
# Label groups
targets <- read.csv("cgtargetsv2.csv", header = TRUE)
targets <- targets[1:14,]
targets

# Combine all the experimental factors into one combined factor
Group <- factor(paste(targets$type, targets$beach, sep = "."))
#Group <- factor(targets$type)
#cbind(targets, Group = Group)
Group

#targets$beach = as.factor(targets$beach)
#targets$type = as.factor(targets$type)
#
Beach <- factor(targets$beach, levels = c("A","C", "E", "B", "D", "F"))
Type <- factor(targets$type, levels = c("CG", "UM"))
Plot <- factor(targets$plot, levels = c("1", "2", "3"))
data.frame(Type, Beach, Plot)
#Group <- factor(targets$type)
#cbind(targets, Group = Group)



# Data Exploration:
#### BCV/PCA ####

# examine  the  samples  for  outliers  and  for  other relationships.
# The  function plotMDS produces  a  plot  in  which  distances  between  samples
# correspond to leading biological coefficient of variation (BCV) between those samples:
#
plotMDS(y)

# Plot with colors, and a legend
points <- c(0,1,2,15,16,17)
#
plotMDS(y, pch = points[Beach], ylab = "Leading logFC dim 2",
        xlab = "Leading logFC dim 1", main = "")
legend("bottomright", legend = levels(Group), pch = points, ncol =2)
# 
#
# MDS plot "gene.selection" option
plotMDS(y, gene.selection = "common", pch = points[Beach], ylab = "Leading logFC dim 2",
        xlab = "Leading logFC dim 1", main = "")
legend("topright", legend = levels(Group), pch = points, ncol =2)
# 