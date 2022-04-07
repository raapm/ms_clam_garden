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

#### 00. Front Matter ####

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

## User set filenames:
uniprot.FN <- "02_input_data/project155.uniprot_blastp.txt.gz"
contig_annot.FN <- "02_input_data/cgrnaseqBTv1.csv"
input.FN <- "02_input_data/out.matrix.csv"
pheno.FN <- "02_input_data/cg_sediment_data_2022-03-25.csv"

## User-set variables
min.reads.mapping.per.transcript <- 10 # Variable to find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.ind <- 5 # choose the minimum number of individuals that need to pass the threshold

## Plotting options
options(scipen = 9999999)

#### 01a. Import Annotation ####
# contig ID and uniprot ID
id_and_uniprot_id.df <- read.table(file = gzfile(uniprot.FN), sep = "\t", header = F)
colnames(id_and_uniprot_id.df) <- c("contig.id", "uniprot.id")
head(id_and_uniprot_id.df)
dim(id_and_uniprot_id.df)

# contig ID and general annotation
RPKM_and_annot.df <- read.table(file = contig_annot.FN, header = T, sep = ",")
dim(RPKM_and_annot.df)
head(RPKM_and_annot.df)

RPKM_and_annot.df <- RPKM_and_annot.df[, which(colnames(RPKM_and_annot.df) %in% c("unigene", "ID", "Evalue", "bp"))] # which columns have useful annotation? 
head(RPKM_and_annot.df)
dim(RPKM_and_annot.df)

# Combine contig ID and uniprot ID
annot.df <- merge(x = id_and_uniprot_id.df, RPKM_and_annot.df, by.x = "contig.id", by.y = "unigene", all.y = T)
annot.df <- annot.df[,c("contig.id", "uniprot.id", "bp", "ID", "Evalue")] # reorder cols
head(annot.df)
dim(annot.df)

#### 01b. Import Phenotypes ####
phenos.df <- read.csv(file = pheno.FN)
phenos.df$sample.id <- paste0("P", phenos.df$plot)

head(phenos.df)


#### 02. Import count data ####
# Import counts file (i.e. the final output from Simple_reads_to_counts repo)
my.counts <- read.csv(file = "02_input_data/out.matrix.csv")
my.counts[1:5,1:5] # note that here the identifier is referred to as 'transcript.id', but is the same as 'contig.id' above

# Set first column as rownames
rownames(my.counts) <- my.counts[,1]
my.counts[1:5, 1:5]
dim(my.counts)
my.counts <- round(x = my.counts, digits = 0) # Round all effective counts to the nearest whole number

# Confirm this did not impact the transcript.id 
table(my.counts$transcript.id==rownames(my.counts)) # All rows should show TRUE here

# Add annotation information count data
my.counts <- merge(x = my.counts, y = annot.df, by.x = "transcript.id", by.y = "contig.id", all.x = T) # use all.x = T in order to ensure contigs that are not in the annotation are retained
dim(my.counts)
head(my.counts)

# Note: we lose rownames here, so add them again
rownames(my.counts) <- my.counts[,"transcript.id"]
my.counts[1:5,1:5]

# Finalize dataframe formats
my.counts.round <- my.counts # Use variable as per rest of script
head(my.counts.round)
str(my.counts.round) # Sample info should be numeric (NOTE: #TODO: This was due to using 'round' on the full dataset above)
my.counts.round$transcript.id <- as.character(my.counts.round$transcript.id)
head(my.counts.round)
my.counts.round[1:5,1:5]


#### 03. Prepare datasets (all, gill, dig) ####
# Create dataframe for all samples
my.counts.round.all <- my.counts.round

# Create dataframe for gill samples
gill.cols <- grep(pattern = "\\.CG\\_G", x = colnames(my.counts.round))
annot.cols <- which(colnames(my.counts.round) %in% c("transcript.id", "uniprot.id", "bp", "ID", "Evalue"))
retain.cols <- c(gill.cols, annot.cols)

my.counts.round.gill <- my.counts.round[, retain.cols]
dim(my.counts.round.gill)
colnames(my.counts.round.gill)

# Create dataframe for digestive gland (dig) samples
dig.cols <- grep(pattern = "\\.CG\\_DG", x = colnames(my.counts.round))
annot.cols <- which(colnames(my.counts.round) %in% c("transcript.id", "uniprot.id", "bp", "ID", "Evalue"))
retain.cols <- c(dig.cols, annot.cols)

my.counts.round.dig <- my.counts.round[, retain.cols]
dim(my.counts.round.dig)
colnames(my.counts.round.dig)

# Create vector for all three datatypes
datatypes <- c("all", "gill", "dig")

# Create list with all three datatypes
input_dataframes.list <- list()
input_dataframes.list[["all"]] <- my.counts.round.all
input_dataframes.list[["gill"]] <- my.counts.round.gill
input_dataframes.list[["dig"]] <- my.counts.round.dig


#### 04. Create DGELists and Filter  ####
# Use the list above and the three different datatypes to produce three different analyses, depending on the tissue type(s) included
datatypes

doi <- NULL; doi.summary <- list(); cpm.filt <- NULL; keep <- NULL; doi.DGEList.filt <- list()
for(i in 1:length(datatypes)){
  
  # Select the datatypes
  doi <- input_dataframes.list[[datatypes[i]]]
  
  # create DGElist, with both counts and genes
  doi.DGEList <- DGEList(counts = doi[, grep(pattern = "eff.counts", x = colnames(doi))]
                         , genes = doi[, grep(pattern = "eff.counts", x = colnames(doi), invert = T)]
                         )
  
  # How many columns of annotation info? 
  dim(doi.DGEList$genes)
  
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


#### 5. Normalization and Data Visualization ####
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


# #### 06. Exploratory MDS plots, remains early draft (bjgs) ####
# datatypes
# #doi <- doi.DGEList.filt[["all"]] # choose from
# #doi <- doi.DGEList.filt[["gill"]] # choose from
# #doi <- doi.DGEList.filt[["dig"]] # choose from
# 
# # What are the samples, and in what order?
# sample_order.df <- rownames(doi$samples)
# sample_order.df <- as.data.frame(sample_order.df)
# colnames(sample_order.df)[1] <- "RNAseq.id"
# sample_order.df
# 
# #sample_order.df$match.id <- sample_order.df$RNAseq.id
# 
# sample_order.df <- separate(data = sample_order.df, col = "RNAseq.id"
#                       , sep = "_P", into = c("drop", "plot_and_suffix")
#                        #, remove = FALSE
#                       )
# 
# sample_order.df <- separate(data = sample_order.df, col = "plot_and_suffix"
#                             , sep = "\\.eff", into = c("plot", "suffix")
#                             #, remove = FALSE
# )
# 
# sample_order.df
# sample_order.df$sample.id <- paste0("P", sample_order.df$plot)
# 
# # Put the phenotypes dataframe into the same order as the samples in the DGEList
# phenos_for_present_samples.df <- merge(x = sample_order.df, y = phenos.df, by = "sample.id")
# 
# ordered_phenos.df <- phenos_for_present_samples.df[order(match(phenos_for_present_samples.df[,"sample.id"], sample_order.df[,"sample.id"])),]
# head(ordered_phenos.df)
# 
# ordered_phenos.df$sample.id
# sample_order.df$sample.id
# 
# # Selecting custom labels for MDS Plot 
# # What do we have to choose from? 
# colnames(ordered_phenos.df)
# 
# # Type, Beach, sample.id
# custom_labels <- paste0(ordered_phenos.df$Type, "_"
#                         , ordered_phenos.df$beach, "_"
#                         , ordered_phenos.df$sample.id 
#        )
# 
# # Survival, Beach, sample.id
# custom_labels <- paste0(ordered_phenos.df$surv, "_"
#                         , ordered_phenos.df$beach, "_"
#                         , ordered_phenos.df$sample.id 
# )
# 
# # Sand, silt, sample.id
# custom_labels <- paste0(ordered_phenos.df$sand, "_"
#                         , ordered_phenos.df$silt, "_"
#                         , ordered_phenos.df$sample.id 
# )
# 
# # Day, Carbon, sample.id
# custom_labels <- paste0(ordered_phenos.df$day, "_"
#                         , ordered_phenos.df$carb, "_"
#                         , ordered_phenos.df$sample.id 
# )
# 
# 
# # Plot
# plotMDS(x = doi, cex= 0.8, labels = custom_labels)
# # save out as 5 x 8


#### 07. Tissue-Specific Expression ####
datatypes
gill.DGEList <- doi.DGEList.filt[["gill"]]
dig.DGEList  <- doi.DGEList.filt[["dig"]] 

expr_gill.vec <- rownames(gill.DGEList$counts)
expr_dig.vec  <- rownames(dig.DGEList$counts)
length(expr_gill.vec)
length(expr_dig.vec)

gill_specific_genes.vec <- setdiff(x = expr_gill.vec, y = expr_dig.vec)
dig_specific_genes.vec  <- setdiff(x = expr_dig.vec, y = expr_gill.vec)

length(gill_specific_genes.vec)
length(dig_specific_genes.vec)

# Data check: 
head(gill_specific_genes.vec)
"25341993" %in% rownames(gill.DGEList$counts) # to confirm the correct reading of setdiff


#### 08. Export background list (expressed genes)
write.table(x = gill.DGEList$genes, file = "04_txomic_results/background_gene_list_gill.txt", sep = "\t", quote = F)
write.table(x = gill.DGEList$genes, file = "04_txomic_results/background_gene_list_gill.txt", sep = "\t", quote = F)


#### 09. Differential expression ####
datatypes
dig.DGEList

#### HERE ####

# Estimate dispersions
dig.DGEList <- estimateCommonDisp(dig.DGEList)
dig.DGEList <- estimateTagwiseDisp(dig.DGEList)
dig.DGEList <- estimateGLMTagwiseDisp(dig.DGEList)
#TODO: confirm these do not write over each other

dig.DGEList

#gill.DGEList

# Phenos are here 
head(phenos.df)

# Use previously developed method to connect the dataframe based on the order of the samples in the DGElist
dig.DGEList$samples

# Then use this in a model matrix command, as per: 

##### v.0.1 code from CG_edgeR_AB_CD_EF_ANOVA.R #####
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
























#### OLD CODE ####
# Want to look at interactions between time and treatment
design <- model.matrix(~ 0 + Group)
#design <- model.matrix(~targets.full$Treat * targets.full$Time)
#design <- model.matrix(~targets$Treat + targets$Time + targets$Treat:targets$Time)
#fit <- glmQLFit(y, design)
colnames(design)
#colnames(design) <- levels(targets.full$Group)

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


heatmap(tab)
#### \END\ OLD CODE ####
# Many more examples provided
















