## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

#

# clear workspace
rm(list=ls())
#
if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

BiocManager::install("limma")

#
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

BiocManager::install("edgeR")
#
getwd()
#install.packages("edgeR")
#install.packages("limma")

citation("edgeR")
library(limma)
library(edgeR)

cgrawdata <- read.csv("cgrnaseqHML.csv", header = TRUE)
# Create A DGEList object for easy manipulation
# y <- DGEList(counts = cgrawdata[ ,2:33], genes = cgrawdata[ ,1])



colnames(cgrawdata)[1] <- "jong.id"


# sort by jong.id
cgrawdata.sorted <- cgrawdata[order(cgrawdata$jong.id),]
str(cgrawdata.sorted)


# head(cgrawdata)

# ?DGEList
# ncol(cgrawdata)
# dim(cgrawdata)
#names(cgrawdata)
#ncol(cgrawdata)
#summary(cgrawdata)
#?data

# 

# 

# Create A DGEList object for easy manipulation
# y <- DGEList(counts = cgrawdata[ ,2:33], genes = cgrawdata[ ,1])

#colnames(y)

#### Merging Identifiers ####

# Merge the identifier from differential gene lists with uniprot information

# setwd("~/")
# deg <- read.table(file = " name of DE list"
#                 , header = T)

# str(deg)
# colnames(deg)[1] <- "jong.id"

# sort(deg)

# sort by jong.id
# deg.sorted <- deg[order(deg$jong.id),]
# str(deg.sorted)

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

#### Annotation ####
# install.packages("org.HS.eg.db") not available for R version 3.3.3
# library(org.HS.eg.db)




###### Subset samples for tissue #####
# Select only Digestive gland individuals

y.gill <- y[,1:14]
colnames(y.gill)

y.dig <- y[,15:28]
colnames(y.dig)

y.gCGRef <- y[,29:30]
colnames(y.gCGRef)

y.dCGRef <- y[,31:32]
colnames(y.dCGRef)

y.all <- y[,1:32]
colnames(y.all)

y <- y.gill
#y <- y.all
#y <- y.dig 
#y <- y.gCGRef
#y <- y.dCGRef




##### BEGIN ####

#### Annotation ####
#install.packages("org.Hs.eg.db")
#require(org.Hs.eg.db)
#idfound <- y$genes$RefSeqID %in% mappedRkeys(org.Hs.egREFSEQ)
#y <- y[idfound,]
#dim(y)


# view DGEList
# str(y)
# dim(y)

#### Library Sizes ####

str(y[[2]])
y[[2]]$lib.size # lib sizes

# str(y[[7]]) # subscript out of bounds error message

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

par(mar = rep(2, 4))
plotMDS(y)
plotMDS(y, gene.selection = "common" ) #  PCA 2 plot 
summary(cpm(y))
# Dimension one separates the 2 different tissue types but the Exposed 
# and Control samples don't appear to separate out in dimension 2. 

##########################################33
# # Annotation

##### Filtering #####
# Filtering out genes with very low counts
dim(cpm(y))
cpm(y)[1:10,1:5]
(y)[1:10, 1:5] # An object of class "DGEList"


keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
dim(y)
colnames(y)

summary(cpm(y))
str(y)



#### Export background gene list ####

write.table(y$genes, file = " CG_DG_background_gene_list_OCT24_2021.txt")

##############
# Clustering, heatmaps etc..

# logcpm <- cpm(y, prior.count = 2, log = TRUE)
# plotMDS(logcpm)


#### Estimation Dispersion ####
# Gives a common dispersion
y <- estimateCommonDisp(y)

# estimate tagwise and common in one run (recommended)
#y <- estimateDisp(y)

# estimate Tagwise dispersions
y <- estimateTagwiseDisp(y)
y <- estimateGLMTagwiseDisp(y)

##################################################################

##### Differential Expression ####



## Differential Expression with all combinations of multiple factors

############################################################

# Defining each treatment combination as a group

# The data frame targets describes the treatment conditions applied to each sample

targets <- read.csv("cgtargetsv1.csv", header = TRUE)
targets <- targets[1:14,]

# Combine all the experimental factors into one combined factor
Group <- factor(paste(targets$type, targets$beach, sep = "."))
cbind(targets, Group = Group)
#targets.full <- cbind(targets, Group = Group)

# Set Control as the reference level
#targets.full$Treat <- relevel(targets$Treat, ref = "Control")
#targets$Treat <- relevel(targets$Treat, ref = "Control")

# Each treatment time for each treatment is a group

#################################################################
####  Form design matrix  ####
################################################################


# Try with linear models: 
# Want to look at interactions between time and treatment
design <- model.matrix(~ 0 + Group)
#design <- model.matrix(~targets.full$Treat * targets.full$Time)
#design <- model.matrix(~targets$Treat + targets$Time + targets$Treat:targets$Time)
#fit <- glmQLFit(y, design)
colnames(design)

#colnames(design) <- levels(targets.full$Group)

#### Likelihood ratio test ####
# fit <- glmFit(y, design) 


# qlf <- glmQLFTest(fit, (fit)

#### Quasi-Likelihood F test ####
fit <- glmQLFit(y, design)

# Estimating the dispersion
# y <- estimate.QLFDisp(y, design, robust = TRUE)
y$common.dispersion
plotBCV(y)
# y <- estimateQLFTagwiseDisp(y)
# y$tagwise.dispersion

##########################################################

# Comparisons we might wish to make:

########################################################33
my.contrasts <- makeContrasts(
  ABvsCD = (GroupCG.A + GroupRef.B) - (GroupCG.C + GroupRef.D),
  ABvsEF = (GroupCG.A + GroupRef.B) - (GroupCG.E + GroupRef.F),
  CDvsEF = (GroupCG.C + GroupRef.D) - (GroupCG.E + GroupRef.F), levels = design)
#
fit <- glmQLFit(y, design)

anov <- glmQLFTest(fit, contrast = my.contrasts)
topTags(anov, n= 23200)

o <- order(anov$table$PValue)
cpm(y)[o[1:23200],]
tab <- cpm(y)[o[1:23200],]
write.table(tab, file = "CG_G_AB_CD_EF_ANOVA_cpm_top genes_Oct24.2021.txt")


tab <- topTags(anov, n = 23200)
write.table(tab, file = "CG_G_AB_CD_EF_Anova_Oct24.2021.txt")
#
tab <- cpm(y)[o[1:150],]
#
#cn <- c("DR18","DC15","DC4","DC5","DC6","DR10", "DR11","DR12","DC1","DC2","DC3", "DR7", "DR8", "DR9") 
#

heatmap(tab)



#my.contrasts <- makeContrasts( 
#CG.1vsRef.4 = GroupCG.1 - GroupRef.4,
#CG.2vsRef.3 = GroupCG.2 - GroupRef.3,
#CG.5vsRef.6 = GroupCG.5 - GroupRef.6,
#CGvsRef = (GroupCG.1 + GroupCG.2 + GroupCG.5) - (GroupRef.4 + GroupRef.3 + GroupRef.6), levels = design)
#
my.contrasts <- makeContrasts(
  CGvRef = (GroupCG.1 + GroupCG.2 + GroupCG.5) - (GroupRef.3 = GroupRef.4 + GroupRef.6),
  levels = design)


# Find genes responding to Mp after 3h
# qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Mp.3hvs0h"])
# topTags(qlf)

# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 

# Find genes responding to Mp after 15 days
# qlf <- glmQLFTest(fit, contrast = my.contrasts[,"Mp.15dvs0h"])
# topTags(qlf)

# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 

# Find genes in the Control group that respond after 3 h
# qlf <- glmQLFTest(fit, contrast = my.contrasts[, "C.3hvs0h"])
# topTags(qlf)

# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 


# Find genes in the control group that respond after 15days
# qlf <- glmQLFTest(fit, contrast = my.contrasts[, "C.15dvs0h"])
# topTags(qlf)

# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 



# Find genes that are different between KB25CG and RefSWC
qlf <- glmQLFTest(fit, contrast = my.contrasts[, "CG.1vsRef.4"])
topTags(qlf, n = 1000 )

# Export DE genes list 
tab <- topTags(qlf, n = 1000)
write.table(tab, file = "CG_DE_KB25_SWC.txt")

# tr <- glmTreat(fit, contrast = my.contrasts[, "MpvsC.0h"], lfc = log2(1.2))
# topTags(tr)


# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 

# Find genes that are DE between CGBB and KB08 REf
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"CG.2vsRef.3"])
topTags(qlf, n = 1000)

o <- order(qlf$table$PValue)
cpm(y)[o[1:500],]
tab <- cpm(y)[o[1:500],]
write.table(tab, file = "CG_DE_BBvsKB08_cpm_top genes.txt")

# Export DE genes list 
tab <- topTags(qlf, n = 1000)
write.table(tab, file = "CG_DE_BBvsKB08.txt")

# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down. 

# Find genes that are DE between KB11vsREfKB11
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"CG.5vsRef.6"])
topTags(qlf, n= 1000)

o <- order(qlf$table$PValue)
cpm(y)[o[1:500],]

tab <- cpm(y)[o[1:500]]
write.table(tab, file = "CG_DE_KB11s_cpm_top_genes.txt")

# Export DE genes list 
tab <- topTags(qlf, n = 1000)
write.table(tab, file = "CG_DE_KB11.txt")

# Find genes that are DE between CGs and Refs
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"CGvsRef"])
topTags(qlf, n= 1000)

o <- order(qlf$table$PValue)
cpm(y)[o[1:500],]

tab <- cpm(y)[o[1:500]]
write.table(tab, file = "CG_DE_CGsvsRefs_cpm_top_genes.txt")

# Export DE genes list 
tab <- topTags(qlf, n = 1000)
write.table(tab, file = "CG_DE_CGsvsRefs.txt")


####################################################

# Use glmTreat to narrow down the list of DE genes and focus on genes that are more
# biological meaningful - DE above a fold-change threshold

# tr <- glmTreat(fit, contrast = my.contrasts[,"MpvsC.15d"], lfc = log2(1.2))
# topTags(tr)
#
# plot all the logFC's against average count size, highlighting the DE genes
# plotMD(qlf)
# abline(h=c(-1,1), col = "blue")
# The blue lines indicate 2 fold up or down.

######################################################


########################################################

#### ANOVA-Like testing ####

############################################################

# anov <- glmQLFTest(fit, contrast = my.contrasts)
# topTags(anov)

# only gave first 4 contrasts

my.contrasts <- makeContrasts( 
  CG.1vsRef.4 = GroupCG.1 - GroupRef.4,
  CG.2vsRef.3 = GroupCG.2 - GroupRef.3,
  CG.5vsRef.6 = GroupCG.5 - GroupRef.6,
  CGvsRef = (GroupCG.1 + GroupCG.2 + GroupCG.5) - (GroupRef.4 + GroupRef.3 + GroupRef.6), levels = design)



anov <- glmQLFTest(fit, contrast = my.contrasts)
topTags(anov, n= 20000)

o <- order(qlf$table$PValue)
cpm(y)[o[1:20000],]
tab <- cpm(y)[o[1:20000],]
write.table(tab, file = "CG_ANOVA_cpm_top genes.txt")


tab <- topTags(anov, n = 20000)
write.table(tab, file = "CG_DE_Anova.txt")

?heatmap


heatmap(y)

# my.contrasts <- makeContrasts( 
# Mp.15dvs0h = Mp.15d - Mp.0h,
# C.15dvs0h = Control.15d - Control.0h,
#MpvsC.0h = Mp.0h - Control.0h,
#MpvsC.3h = (Mp.3h - Mp.0h) - (Control.3h - Control.0h),
#MpvsC.15d = (Mp.15d - Mp.0h) - (Control.15d - Control.0h), levels = design)

# anov <- glmQLFTest(fit, contrast = my.contrasts)
# topTags(anov, n = 40)


################################################################

################################################################
### Differential Expression

#### The design matrix

Treat <- factor(c("C", "E", "E", "C", "C", "E", "E", "C", "E", "C", "C", "E", "E", "C", "C", "E", "C", "E", "C", "E", "E", "C", "C", "E" ))
Time <- factor(c("0", "0", "0", "0", "0", "0","0", "0", "3h", "3h", "3h", "3h", "3h", "3h","3h", "3h", "15d", "15d", "15d", "15d", "15d", "15d", "15d", "15d"))

data.frame(Sample=colnames(y), Treat, Time)
design <- model.matrix(~Treat+ Time)
rownames(design) <- colnames(y)
design

### Estimating the dispersion

# y <- estimateDisp(y)


y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

plotBCV(y)
### Differential Expression
fit <- glmFit(y, design)

## Conduct likelihood ratio tests for control vs exposed time differences and show top genes

lrt <- glmLRT(fit)
topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col ="blue")

##########################################################



-------------------------------------------------------------------
  
  
  
  
  