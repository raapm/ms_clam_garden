# Script to analyze transcriptomic data
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

options(scipen = 9999999)


#### Method 1. cpm data ####
# Set input filenames
input.FN <- "02_input_data/out.matrix.csv"
pheno.FN <- "02_input_data/cg_sediment_data_2022-03-25.csv" # the original input filename

#### Read in phenotype data ####
phenos.df <- read.csv(file = pheno.FN)
phenos.df$sample.id <- paste0("P", phenos.df$plot)

phenos.df

#### Read in raw count data ####
# Import counts file (i.e. the final output from Simple_reads_to_counts repo)
my.counts <- read.csv(file = "02_input_data/out.matrix.csv")
my.counts[1:5,1:5]

# Set first column as rownames
rownames(my.counts) <- my.counts[,1]
my.counts[1:5, 1:5]

# Round all expression levels to no decimal place, remove transcript.id column (now rownames)
my.counts.round <- round(my.counts[,-1])
str(my.counts.round) # should be numeric
my.counts.round[1:5, 1:5]

# Prepare individual tissues datasets, and combined tissues (for MDS)
# Create dataframe for all samples
my.counts.round.all <- my.counts.round

# Create dataframe for gill samples
my.counts.round.gill <- my.counts.round[, grep(pattern = "\\.CG\\_G", x = colnames(my.counts.round))]
dim(my.counts.round.gill)
colnames(my.counts.round.gill)

# Create dataframe for digestive gland (dig) samples
my.counts.round.dig <- my.counts.round[, grep(pattern = "\\.CG\\_G", x = colnames(my.counts.round)
                                              , invert = T)]
dim(my.counts.round.dig)
colnames(my.counts.round.dig)

# Create vector for all three datatypes
datatypes <- c("all", "gill", "dig")

# Create list with all three datatypes
input_dataframes.list <- list()
input_dataframes.list[["all"]] <- my.counts.round.all
input_dataframes.list[["gill"]] <- my.counts.round.gill
input_dataframes.list[["dig"]] <- my.counts.round.dig


#### Per datatype analysis, filter ####
# User-set variables
min.reads.mapping.per.transcript <- 10 # Variable to find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.ind <- 5 # choose the minimum number of individuals that need to pass the threshold

# Use the list above and the three different datatypes to produce three different analyses, depending on the tissue type(s) included
datatypes

doi <- NULL; doi.summary <- list(); cpm.filt <- NULL; keep <- NULL; doi.DGEList.filt <- list()
for(i in 1:length(datatypes)){
  
  # Select the datatypes
  doi <- input_dataframes.list[[datatypes[i]]]
  
  # create DGElist
  doi.DGEList <- DGEList(counts = doi)
  
  # Statistics on library sizes
  doi.summary[[paste0(datatypes[[i]], "_align_summary")]]    <- summary(doi.DGEList$samples$lib.size)
  doi.summary[[paste0(datatypes[[i]], "_align_summary_sd")]] <- sd(doi.DGEList$samples$lib.size)
  
  # Determine cpm.filt for the library with the fewest transcripts to be the most stringent on low expression cutoff
  cpm.filt <- min.reads.mapping.per.transcript / min(doi.DGEList$samples$lib.size) * 1000000
  print(paste0("The cpm filter threshold for ", datatypes[i], " is determined as: ", cpm.filt)) # min cpm filt
  doi.summary[[paste0(datatypes[[i]], "_cpm.filt")]] <- cpm.filt
  
  # Identify tags passing filter
  keep <- rowSums(cpm(doi.DGEList)>cpm.filt) >= min.ind # Find which transcripts pass the filter
  table(keep) # gives number passing, number failing
  doi.summary[[paste0(datatypes[[i]], "_post-filter_retained")]] <- table(keep)
  
  # Subset DGEList to only those passing expression filter
  doi.DGEList <- doi.DGEList[keep, , keep.lib.sizes=FALSE] #keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
  dim(doi.DGEList)
  
  # Save out the filtered DGEList
  doi.DGEList.filt[[datatypes[i]]] <- doi.DGEList
  
  print("Your data is saved as a slot in doi.DGEList.filt")
  print("Summaries of raw data and filtering steps are in doi.summary")
  
}


#### 3. Normalization and Data Visualization ####
# Use the list above (doi.DGEList.filt) and the three different datatypes to produce three different analyses, depending on the tissue type(s) included
datatypes

doi <- NULL
for(i in 1:length(datatypes)){
  
  # Select the datatypes
  doi <- doi.DGEList.filt[[datatypes[i]]]
  
  # Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
  doi <- calcNormFactors(doi, method = c("TMM"))
  doi$samples
  
  # Statistics on library sizes
  doi.summary[[paste0(datatypes[[i]], "_norm_factors")]]    <- doi$samples
  
  # Plot normalization factors
  pdf(file = paste0("04_txomic_results/", datatypes[i], "_norm_factors.pdf"), width = 6, height = 6)
  plot(doi$samples$norm.factors ~ doi$samples$lib.size
       , xlab = "Library Size", ylab = paste0("Normalization factors (", datatypes[i], ")")
       , las = 1)
  dev.off()
  
  # Estimate dispersions (measure inter-library variation per tag)
  doi <- estimateDisp(doi) # note that this can use a design matrix when provided 
  doi$prior.df # est. overall var. across genome for dataset
  doi.summary[[paste0(datatypes[[i]], "_overall_var")]] <- doi$prior.df
  
  print("The coefficient of variation for biological variation for this dataset is: ")
  sqrt(doi$common.disp) #coeff of var, for biol. var
  doi.summary[[paste0(datatypes[[i]], "_coeff_var_biol")]] <- sqrt(doi$common.disp)
  
  # Plot Biological coefficient of variation by mean log cpm
  pdf(file = paste0("04_txomic_results/", datatypes[i], "_BCV_by_mean_log_cpm.pdf"), width = 6, height = 6)
  plotBCV(doi, ylab = paste0("Biological coefficient of variation (", datatypes[i] ,")"))
  dev.off()
  
  # MDS plot
  pdf(file = paste0("04_txomic_results/", datatypes[i], "_mds_dim_1_2.pdf"), width = 9, height = 5)
  plotMDS(x = doi, cex= 0.8) # note that this is supposed to be run on whatever you wrote calcNormFact() to
  dev.off()
  
  # MDS plotting options to work with:
  # labels = character vector of sample names or labels. Defaults to colnames(x).
  # pch = plotting symbol or symbols. See points for possible values. Ignored if labels is non-NULL.
  # TODO: this will be adjusted once the phenotypes are brought in
  
  }


#### Exploratory MDS plots ####
datatypes
#doi <- doi.DGEList.filt[["all"]] # choose from
#doi <- doi.DGEList.filt[["gill"]] # choose from
#doi <- doi.DGEList.filt[["dig"]] # choose from

# What are the samples, and in what order?
sample_order.df <- rownames(doi$samples)
sample_order.df <- as.data.frame(sample_order.df)
colnames(sample_order.df)[1] <- "RNAseq.id"
sample_order.df

#sample_order.df$match.id <- sample_order.df$RNAseq.id

sample_order.df <- separate(data = sample_order.df, col = "RNAseq.id"
                      , sep = "_P", into = c("drop", "plot_and_suffix")
                       #, remove = FALSE
                      )

sample_order.df <- separate(data = sample_order.df, col = "plot_and_suffix"
                            , sep = "\\.eff", into = c("plot", "suffix")
                            #, remove = FALSE
)

sample_order.df
sample_order.df$sample.id <- paste0("P", sample_order.df$plot)

# Put the phenotypes dataframe into the same order as the samples in the DGEList
phenos_for_present_samples.df <- merge(x = sample_order.df, y = phenos.df, by = "sample.id")

ordered_phenos.df <- phenos_for_present_samples.df[order(match(phenos_for_present_samples.df[,"sample.id"], sample_order.df[,"sample.id"])),]
head(ordered_phenos.df)

ordered_phenos.df$sample.id
sample_order.df$sample.id

# Selecting custom labels for MDS Plot 
# What do we have to choose from? 
colnames(ordered_phenos.df)

# Type, Beach, sample.id
custom_labels <- paste0(ordered_phenos.df$Type, "_"
                        , ordered_phenos.df$beach, "_"
                        , ordered_phenos.df$sample.id 
       )

# Survival, Beach, sample.id
custom_labels <- paste0(ordered_phenos.df$surv, "_"
                        , ordered_phenos.df$beach, "_"
                        , ordered_phenos.df$sample.id 
)

# Sand, silt, sample.id
custom_labels <- paste0(ordered_phenos.df$sand, "_"
                        , ordered_phenos.df$silt, "_"
                        , ordered_phenos.df$sample.id 
)

# Day, Carbon, sample.id
custom_labels <- paste0(ordered_phenos.df$day, "_"
                        , ordered_phenos.df$carb, "_"
                        , ordered_phenos.df$sample.id 
)


# Plot
plotMDS(x = doi, cex= 0.8, labels = custom_labels)
# save out as 5 x 8


#### END NEW ####


#### Method 2. Original Files ####
# Goes from input data through normalization to MDS plot
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
