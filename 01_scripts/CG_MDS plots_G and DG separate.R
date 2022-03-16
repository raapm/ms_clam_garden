# Data exploration script
# Initialized to GitHub 2022-03-15

# clear workspace
#rm(list=ls())

# Install packages and load libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("edgeR")
# BiocManager::install("limma")

library("limma")
library("edgeR")

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

# Set input filenames
input.FN <- "02_input_data/cgrnaseqBTv1.csv" # this is the input filename for "CG_MDS plots_G and DG separate.R"
annot.FN <- "02_input_data/project155.uniprot_blastp.txt.gz"

#### Read in all libraries data ####
cgrawdata <- read.csv(file = input.FN, header = TRUE)
head(cgrawdata)

# All columns
colnames(cgrawdata)

# Selected for DGEList
colnames(cgrawdata)[2:33]

# Rename unigene to 'jong.id'
colnames(cgrawdata)[1] <- "jong.id"

# Sort by unigene
# cgrawdata.sorted <- cgrawdata[order(cgrawdata$jong.id),]
# str(cgrawdata.sorted)
## to confirm: I think this is unecessary? 

#### Read in annotation information ####
uniprot.map <- read.table(file = gzfile(annot.FN)
                          , sep = "\t"
                          , header = F
                          )

head(uniprot.map)
colnames(uniprot.map) <- c("jong.id" , "uniprot")
uniprot.map$jong.id <- as.character(uniprot.map$jong.id)

head(uniprot.map)
str(uniprot.map)
dim(uniprot.map)


#### Combine counts and annotation
data <- merge(x = cgrawdata, y = uniprot.map, by = "jong.id")
dim(cgrawdata)     # 52000
dim(uniprot.map)   # 42708
dim(data) # limited by uniprot.map

colnames(data)

#### Question: why are there fewer entries in the uniprot.map than the counts? ####

#### Prepare DGEList ####
y <- DGEList(counts = data[, which(!(colnames(data) %in% c("Evalue", "bp", "ID", "uniprot", "jong.id")))] # select only count cols
             , genes = data[, c("jong.id", "Evalue", "bp", "ID", "uniprot")] # select only annotation cols
             ) 

str(y)
head(y$genes)

# Create a DGEList object
### OTHER METHOD APPLIED: 
# y <- DGEList(counts = cgrawdata[, 2:33], genes = cgrawdata[ ,1])


###### Tissue-separated #####
# Based on numeric selection, it appears that the gill tissues are the following samples: 
colnames(y)[1:14]  # Gill
colnames(y)[15:28] # Digestive gland

# Gill samples are those without '.1' at the end; these have a total of 4 characters in the string
y.gill <- y[, which(nchar(colnames(y))==4)]
colnames(y.gill)

# Digestive gland samples are those with '.1' at the end; these have a total of 6 characters in the string
y.dig <- y[, which(nchar(colnames(y))==6)]
colnames(y.dig)

# y.gillCG <- y[,29]
# colnames(y.gillCG)
# 
# y.gillRef <- y[,30]
# colnames(y.gillRef)
# 
# y.dig.CG <- y[,31]
# colnames(y.dig.CG)
# 
# y.dig.Ref <- y[,32]
# colnames(y.dig.Ref)
# 
# y.all <- y[,1:28]
# colnames(y.all)

##### Manually choose which unit to be analyzed 
#y <- y.gill
#y <- y.all
y <- y.dig 
#y <- y.gillCG
#y <- y.gillRef


#### Library Sizes ####
str(y[[2]])
y[[2]]$lib.size # lib sizes


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