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

## Plotting options
options(scipen = 9999999)

#### Read in annotation information ####
# contig ID and uniprot ID
id_and_uniprot_id.df <- read.table(file = gzfile("02_input_data/project155.uniprot_blastp.txt.gz"), sep = "\t", header = F)
# TODO: add these filenames as user variables above
colnames(id_and_uniprot_id.df) <- c("contig.id", "uniprot.id")
head(id_and_uniprot_id.df)
dim(id_and_uniprot_id.df)

# contig ID and general annotation
RPKM_and_annot.df <- read.table(file = "02_input_data/cgrnaseqBTv1.csv", header = T, sep = ",")
# TODO: add these filenames as user variables above
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


#### Method 1. cpm data ####
# Set input filenames
input.FN <- "02_input_data/out.matrix.csv"
pheno.FN <- "02_input_data/cg_sediment_data_2022-03-25.csv" # the original input filename

#### Read in phenotype data ####
phenos.df <- read.csv(file = pheno.FN)
phenos.df$sample.id <- paste0("P", phenos.df$plot)

head(phenos.df)

#### Read in raw count data ####
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


#### Add annotation information ####
my.counts <- merge(x = my.counts, y = annot.df, by.x = "transcript.id", by.y = "contig.id", all.x = T) # use all.x = T in order to ensure contigs that are not in the annotation are retained
dim(my.counts)
head(my.counts)
# note: we lose rownames here, so add them again
rownames(my.counts) <- my.counts[,"transcript.id"]
my.counts[1:5,1:5]


# # Round all expression levels to no decimal place, remove transcript.id column (now rownames)
# my.counts <- my.counts[,-1]
# 
# # Which columns are sample data cols? 
# sample.cols <- grep(pattern = "eff.counts", colnames(x = my.counts))
# 
# # Which columns are annotation data? 
# annot.cols <- grep(pattern = "eff.counts", colnames(x = my.counts), invert = T)

my.counts.round <- my.counts # Use variable as per rest of script
head(my.counts.round)
str(my.counts.round) # Sample info should be numeric (NOTE: #TODO: This was due to using 'round' on the full dataset above)
my.counts.round$transcript.id <- as.character(my.counts.round$transcript.id)
head(my.counts.round)
my.counts.round[1:5,1:5]

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
  
  # create DGElist, with both counts and genes
  doi.DGEList <- DGEList(counts = doi[, grep(pattern = "eff.counts", x = colnames(doi))]
                         , genes = doi[, grep(pattern = "eff.counts", x = colnames(doi), invert = T)]
                         )
  
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


#### Exploratory MDS plots, remains early draft ####
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


#### 4. Tissue-specific expression ####
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


#### 5. Differential expression ####







