# Script to analyze transcriptomic data
# Initialized to GitHub 2022-03-15

# clear workspace
#rm(list=ls())

#### 00. Front Matter ####
# Install packages and load libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("edgeR")
# BiocManager::install("limma")

library("limma")
library("edgeR")
library("tidyr")

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)

## User set filenames:
uniprot.FN <- "02_input_data/project155.uniprot_blastp.txt.gz"
contig_annot.FN <- "02_input_data/cgrnaseqBTv1.csv"
input.FN <- "02_input_data/out.matrix_cg_2022-06-07.csv"
pheno.FN <- "02_input_data/cg_sediment_phenos_2022-06-24.csv"

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
my.counts <- read.csv(file = input.FN)
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
str(my.counts.round) # transcript.id needs to be made into a character (NOTE: it is numeric due to using 'round' on the full dataset above)
my.counts.round$transcript.id <- as.character(my.counts.round$transcript.id)
str(my.counts.round)
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
  
  # Select the datatype
  doi <- input_dataframes.list[[datatypes[i]]]
  
  # Create DGElist, with both counts and genes
  doi.DGEList <- DGEList(counts = doi[, grep(pattern = "eff.counts", x = colnames(doi))]
                         , genes = doi[, grep(pattern = "eff.counts", x = colnames(doi), invert = T)]
                         )
  
  # How many columns of annotation info? 
  dim(doi.DGEList$genes)
  colnames(doi.DGEList$genes)
  
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



#### 05. Normalization and Data Visualization ####
# Use the list above (doi.DGEList.filt) and the three different datatypes to produce three different analyses, depending on the tissue type(s) included
datatypes

doi <- NULL; doi.DGEList.filt.norm <- list()
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
  doi <- estimateDisp(doi) # note that this can use a design matrix when provided, note: here we are not looking at any specific variable, leave group all as '1'
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
  
  # Assign the current DGElist into a list
  doi.DGEList.filt.norm[[datatypes[i]]] <- doi
  
}

# Output is doi.DGEList.filt.norm
names(doi.DGEList.filt.norm)
doi.DGEList.filt.norm$all$samples
doi.DGEList.filt.norm$gill$samples
doi.DGEList.filt.norm$dig$samples
# norm factors have been added
# Here forward use doi.DGEList.filt.norm


#### 06. MDS plots ####
datatypes

# Plot all samples (both tissues) with concise labels
doi <- doi.DGEList.filt.norm[["all"]]

# MDS plot
# Note: run on whatever you wrote calcNormFact() to
pdf(file = paste0("04_txomic_results/all_samples_mds_dim_1_2.pdf"), width = 8, height = 6)
plotMDS(x = doi
        , labels = 
          gsub(pattern = "\\.eff.*", replacement = ""
               , x = gsub(pattern = ".*\\.CG\\_", replacement = "", x = colnames(doi$counts))
          )
        , cex= 1
        #, xlab = # This is where one would update the axis label, but I think it is fine as is now
        ) 
dev.off()


# Plot each tissue with specified variables of interest
#par(mfrow=c(1,1))

tissues.to.include <- c("gill", "dig")
phenos_for_present_samples.df <- NULL;  ordered_phenos.df <- NULL
for(i in 1:length(tissues.to.include)){
  
  pdf(file = paste0("04_txomic_results/" 
                    , tissues.to.include[i]
                    , "_samples_mds_dim_1_2.pdf"), width = 11.3, height = 5.2)
  
  
  # Reporting
  print(paste0("Plotting ", tissues.to.include[i]))
  
  # Choose the tissue DGEList
  doi <- doi.DGEList.filt.norm[[tissues.to.include[i]]] # choose from 
  
  # What are the samples, and in what order?
  sample_order.df <- rownames(doi$samples)
  sample_order.df <- as.data.frame(sample_order.df)
  colnames(sample_order.df)[1] <- "RNAseq.id"
  sample_order.df
  
  # Clean up sample IDs
  sample_order.df <- separate(data = sample_order.df, col = "RNAseq.id"
                              , sep = "_P", into = c("drop", "plot_and_suffix")
                              #, remove = FALSE
  )
  
  sample_order.df <- separate(data = sample_order.df, col = "plot_and_suffix"
                              , sep = "\\.eff", into = c("plot", "suffix")
                              #, remove = FALSE
  )
  
  sample_order.df
  
  # Add character to allow sample IDs to match the phenotype df
  sample_order.df$sample.id <- paste0("P", sample_order.df$plot)
  
  # Add order of the samples
  sample_order.df$true.order <- seq(1:nrow(sample_order.df))
  
  # Merge the samples (including order) with the phenotypes dataframe
  phenos_for_present_samples.df <- merge(x = sample_order.df, y = phenos.df, by = "sample.id")
  
  # Put back in the order that exists in the DGEList
  ordered_phenos.df <- phenos_for_present_samples.df[order(phenos_for_present_samples.df$true.order), ]
  
  # ordered_phenos.df <- phenos_for_present_samples.df[order(match(phenos_for_present_samples.df[,"sample.id"], sample_order.df[,"sample.id"])),]
  head(ordered_phenos.df)
  
  # Selecting custom labels for MDS Plot
  # What do we have to choose from?
  #colnames(ordered_phenos.df)
  
  # Create ID for plot: beach type, % survival, carbonate, [ sand ], sample ID
  
  # Type, survival percent, carbonate percent, Beach, sample.id
  custom_labels <- paste0(ordered_phenos.df$Type, "_"
                          , ordered_phenos.df$surv, "_"
                          , ordered_phenos.df$carb, "_"
                          , ordered_phenos.df$beach, "_"
                          , ordered_phenos.df$sample.id
  )
  
  # Plot
  plotMDS(x = doi, cex= 0.8, labels = custom_labels
          #, main = tissues.to.include[i]
          )
  dev.off()
  
}


#### 07. Tissue-Specific Expression ####
datatypes

n.genes.per.tissue.for.heatmap <- 200

# Simplify object calls
all.DGEList <- doi.DGEList.filt.norm[["all"]]
gill.DGEList <- doi.DGEList.filt.norm[["gill"]]
dig.DGEList  <- doi.DGEList.filt.norm[["dig"]] 

# Obtain the contig names that are expressed in each DGEList
expr_gill.vec <- rownames(gill.DGEList$counts)
expr_dig.vec  <- rownames(dig.DGEList$counts)
length(expr_gill.vec)
length(expr_dig.vec)

# Compare the DGEList expressed contigs to see tissue-specific contigs
gill_specific_genes.vec <- setdiff(x = expr_gill.vec, y = expr_dig.vec)
dig_specific_genes.vec  <- setdiff(x = expr_dig.vec, y = expr_gill.vec)
length(gill_specific_genes.vec)
length(dig_specific_genes.vec)

# Collect the annotation for tissue-specific genes
gill_specific_annot.df <- gill.DGEList$genes[gill.DGEList$genes$transcript.id %in% gill_specific_genes.vec, ]
dig_specific_annot.df <- dig.DGEList$genes[dig.DGEList$genes$transcript.id %in% dig_specific_genes.vec, ]

# # Data check: 
# head(gill_specific_genes.vec)
# "25341993" %in% rownames(gill.DGEList$counts) # to confirm the correct reading of setdiff

# Export tissue-specific genes and associated annotation
write.table(x = gill_specific_annot.df, file = "04_txomic_results/tissue-specific_genes_gill.txt", sep = "\t", quote = F
            , row.names = F)
write.table(x = dig_specific_annot.df, file = "04_txomic_results/tissue-specific_genes_dig.txt", sep = "\t", quote = F
            , row.names = F)

## Make plot for tissue-specific genes
# Need to subset first, as there are too many currently to make a proper heatmap

# Which genes show tissue-specific expression? 
genes_to_retain.vec <- c(gill_specific_genes.vec, dig_specific_genes.vec)

# Confirm sizes
length(gill_specific_genes.vec) + length(dig_specific_genes.vec) == length(genes_to_retain.vec)
length(genes_to_retain.vec)

# Extract linear cpm for selecting the top expressed
to_id_top_expressed <- cpm(y = all.DGEList, log=FALSE)
to_id_top_expressed.df <- as.data.frame(to_id_top_expressed)
to_id_top_expressed.df[1:5,1:5]
dim(to_id_top_expressed.df)

# Keep only tissue specific
to_id_top_expressed.df <- to_id_top_expressed.df[rownames(to_id_top_expressed.df) %in% genes_to_retain.vec, ]
dim(to_id_top_expressed.df)
length(genes_to_retain.vec) # a couple will be missing (diff cpm filter), can do a setdiff here as needed to see
length(intersect(x = genes_to_retain.vec, y = rownames(to_id_top_expressed.df)))

# Separate into tissue type
gill_linear.df  <- to_id_top_expressed.df[, grep(pattern = "\\.CG\\_G", x = colnames(to_id_top_expressed.df))]
dig_linear.df  <- to_id_top_expressed.df[, grep(pattern = "\\.CG\\_DG", x = colnames(to_id_top_expressed.df))]

# Add up total number of reads aligned to feature
gill_linear.df$sum.counts <- rowSums(gill_linear.df)
dig_linear.df$sum.counts <- rowSums(dig_linear.df)

gill_linear.df[1:2, ]
dig_linear.df[1:2, ]

# Sort by top expressors, use variable set above
gill_linear.df <- gill_linear.df[order(gill_linear.df$sum.counts, decreasing = T), ]
head(gill_linear.df)
top.expr.gill.sp <- rownames(gill_linear.df)[1:n.genes.per.tissue.for.heatmap]

dig_linear.df <- dig_linear.df[order(dig_linear.df$sum.counts, decreasing = T), ]
head(dig_linear.df)
top.expr.dig.sp <- rownames(dig_linear.df)[1:n.genes.per.tissue.for.heatmap]

# Technically I should be able to use these top.expr.dig.sp in the order they are in now
top.expr.gill.sp.df <- as.data.frame(x = top.expr.gill.sp)
colnames(top.expr.gill.sp.df)[1] <- "contig.id"
#top.expr.gill.sp.df$true.order <- seq(1:nrow(top.expr.gill.sp.df))
head(top.expr.gill.sp.df)

top.expr.dig.sp.df <- as.data.frame(x = top.expr.dig.sp)
colnames(top.expr.dig.sp.df)[1] <- "contig.id"
#top.expr.dig.sp.df$true.order <- seq(1:nrow(top.expr.dig.sp.df))
head(top.expr.dig.sp.df)

# Plotting order
plotting_order.df <- rbind(top.expr.gill.sp.df, top.expr.dig.sp.df)
plotting_order.df$true.order <- seq(1:nrow(plotting_order.df))
head(plotting_order.df)
tail(plotting_order.df)

#### Get the actual data for plotting
# Convert data to log2 CPM
logcounts <- cpm(y = all.DGEList, log=TRUE)

# Select only those that were specified to be retained
top_genes_to_retain.vec <- c(top.expr.gill.sp, top.expr.dig.sp)
length(top_genes_to_retain.vec)

logcounts_tissue_specific <- logcounts[rownames(logcounts) %in% top_genes_to_retain.vec, ] # all retained
dim(logcounts_tissue_specific)
logcounts_tissue_specific[1:5, 1:5]

logcounts_tissue_specific.df <- as.data.frame(logcounts_tissue_specific)
#head(logcounts_tissue_specific.df)
logcounts_tissue_specific.df$contig.id <- rownames(logcounts_tissue_specific.df)
head(logcounts_tissue_specific.df)

logcounts_tissue_specific.df <- merge(x = logcounts_tissue_specific.df, y = plotting_order.df, by = "contig.id")
head(logcounts_tissue_specific.df) # Lost rownames, add back from contig.id
rownames(logcounts_tissue_specific.df) <- logcounts_tissue_specific.df$contig.id
head(logcounts_tissue_specific.df) # Got them back

# Order by top expression
logcounts_tissue_specific.df <- logcounts_tissue_specific.df[order(logcounts_tissue_specific.df$true.order), ]
logcounts_tissue_specific.df$true.order
logcounts_tissue_specific.df<- logcounts_tissue_specific.df[, grep(pattern = "contig.id|true.order", x = colnames(logcounts_tissue_specific.df), invert = T)]

colnames(logcounts_tissue_specific.df) <- gsub(pattern = "\\.eff.*", replacement = ""
     , x = gsub(pattern = ".*\\.CG\\_", replacement = "", x = colnames(logcounts_tissue_specific.df))
          )

logcounts_tissue_specific.df <- as.matrix(x = logcounts_tissue_specific.df)

# test <- head(logcounts_tissue_specific.df, n = 20)
# test <- as.matrix(test)

pdf(file = "04_txomic_results/tissue_specific_heatmap.pdf", width = 6, height = 6)
heatmap(x = logcounts_tissue_specific.df
        , Rowv = NA
        , labRow = NA # Removes contig ids
        ) # Red is expression
dev.off()


#### 08. Export background list (expressed genes) ####
write.table(x = gill.DGEList$genes, file = "04_txomic_results/background_gene_list_gill.txt", sep = "\t", quote = F
            , row.names = F)
write.table(x = dig.DGEList$genes, file = "04_txomic_results/background_gene_list_dig.txt", sep = "\t", quote = F
            , row.names = F)

# # Output all dgelist (may not be necessary)
# write.table(x = all.DGEList$genes, file = "04_txomic_results/background_gene_list_gill.txt", sep = "\t", quote = F
#              , row.names = F)


#### 09. Differential expression ####
datatypes
#dig.DGEList
#gill.DGEList

lists_of_interest <- list()
lists_of_interest[["dig"]] <- dig.DGEList
lists_of_interest[["gill"]] <- gill.DGEList


# Assign the bins for carbonate
colnames(phenos.df)

# Carbonate
explan_variable <- "carb"
length(phenos.df[, explan_variable])  # How many samples? 
summary(phenos.df[, explan_variable]) # What are the characteristics of this variable? 

# How many samples in each bin? 
table(phenos.df[, explan_variable] <= 3.25)
table(phenos.df[, explan_variable] > 3.25 & phenos.df[, explan_variable] < 8)
table(phenos.df[, explan_variable] >= 8)

# Assign samples to bins
head(phenos.df)
phenos.df$carb.bins <- rep(NA, times = length(phenos.df[, explan_variable])) # Create vector

phenos.df[which(phenos.df[, explan_variable] <= 3.25), "carb.bins"] <- "low"
phenos.df[which(phenos.df[, explan_variable] >= 8), "carb.bins"] <- "high"
phenos.df[which(phenos.df[, explan_variable] > 3.25 & 
                              phenos.df[, explan_variable] < 8)
                      , "carb.bins"] <- "med"

phenos.df[,c("plot", "carb", "carb.bins")]

# Survival
explan_variable <- "surv"
length(phenos.df[, explan_variable])  # How many samples? 
summary(phenos.df[, explan_variable]) # What are the characteristics of this variable? 
# Replace this with the summary of the actual retained samples

# How many samples in each bin? (note: the values below are based on quartiles for only the retained beaches)
table(phenos.df[, explan_variable] <= 62.5)
table(phenos.df[, explan_variable] > 62.5 & phenos.df[, explan_variable] < 90)
table(phenos.df[, explan_variable] >= 90)

# Assign samples to bins
head(phenos.df)
phenos.df$surv.bins <- rep(NA, times = length(phenos.df[, explan_variable])) # Create vector

phenos.df[which(phenos.df[, explan_variable] <= 62.5), "surv.bins"] <- "low"
phenos.df[which(phenos.df[, explan_variable] >= 90), "surv.bins"] <- "high"
phenos.df[which(phenos.df[, explan_variable] > 62.5 & 
                  phenos.df[, explan_variable] < 90)
          , "surv.bins"] <- "med"

phenos.df[,c("plot", "surv", "surv.bins")]


### Loop to analyze data
names(lists_of_interest)

datatype_current <- NULL
for(i in 1:length(names(lists_of_interest))){
  
  target_list <- lists_of_interest[[i]]
  datatype_current <- names(lists_of_interest)[i]
  
  ### Defining group ###
  head(phenos.df) # Phenos are here 
  
  ## Provide group detail based on order of the samples in the DGEList
  sample_order.df <- as.data.frame(rownames(target_list$samples))
  colnames(sample_order.df) <- "sample"
  sample_order.df
  
  # Remove suffix
  sample_order.df$sample <- gsub(pattern = ".eff.counts", replacement = "", x = sample_order.df$sample)
  head(sample_order.df)
  
  # Isolate <tissue>_<plot> identifier
  sample_order.df <- separate(data = sample_order.df, col = "sample", into = c("prefix", "plot_and_suffix")
                              , sep = "CG_", remove = F)
  head(sample_order.df)
  
  # Separate <tissue>_<plot>
  sample_order.df <- separate(data = sample_order.df, col = "plot_and_suffix", into = c("tissue", "plot")
                              , sep = "_", remove = T
  )

  head(sample_order.df) # This is the order of the samples in the DGEList
  head(phenos.df)       # This is the pheno data (regardless of tissue)
                        # Do not assume these are in same order
  
  # Ensure both are character vectors
  str(sample_order.df)
  str(phenos.df$sample.id)
  
  # Document DGEList sample order using an index
  sample_order.df$DGEList.order <- seq(1:length(sample_order.df$sample))
  head(sample_order.df)
  
  # Combine the pheno data with the DGEList sample names
  samples_and_phenos.df <- merge(x = sample_order.df, y = phenos.df, by.x = "plot", by.y = "sample.id", all.x = T)
  samples_and_phenos.df # ensure there are no missing rows, this is NOT expected to be in order yet
  
  # Return to sample order defined by DGEList.order index
  samples_and_phenos.df <- samples_and_phenos.df[order(samples_and_phenos.df$DGEList.order), ]
  samples_and_phenos.df # Now it should be in order
  
  # Visually compare and confirm
  cbind(sample_order.df$plot, samples_and_phenos.df$plot) # These should be in the same order
  table(sample_order.df$plot == samples_and_phenos.df$plot)
  
  #### DE - phenotype of interest 1: Clam Garden ####
  colnames(samples_and_phenos.df)
  
  # Set the explanatory variable
  categorical_explan_variable <- "Type"
  
  # Assign this phenotype to the group vector in the DGEList
  grouping.vec <- as.factor(samples_and_phenos.df[, categorical_explan_variable])
  grouping.vec <- relevel(x = grouping.vec, ref = "Ref") # Assign the reference as the reference factor level
  target_list$samples$group <- grouping.vec
  
  # To see the assigned group, and the sample names
  target_list$samples$group
  rownames(target_list$samples)
  
  # Estimate dispersions for multiple group experiment
  target_list <- estimateCommonDisp(y = target_list) # Group is already assigned to the DGEList and should not be supplied here or will produce error. The function automatically uses the held group. 
  target_list <- estimateTagwiseDisp(target_list)
  
  # Build model design
  design <- model.matrix(~target_list$samples$group)
  fit <- glmQLFit(y = target_list, design = design)
  results <- glmQLFTest(glmfit = fit)
  all_results <- topTags(object = results
                         , n = length(target_list$genes$transcript.id)
                         , adjust.method = "none"
                         , p.value = 0.001
  )
  
  output.FN <- paste0("04_txomic_results/DE_results_CG_v_REF_", datatype_current, ".txt")
  
  write.table(x = all_results, file = output.FN, quote = F, sep = "\t", row.names = F)
  
  
  #### DE - phenotype of interest 2: carbonate ####
  # Set the explanatory variable
  categorical_explan_variable <- "carb.bins"
  
  # Assign this phenotype to the group vector in the DGEList
  grouping.vec <- as.factor(samples_and_phenos.df[, categorical_explan_variable])
  grouping.vec <- relevel(x = grouping.vec, ref = "low") # Assign the reference as the reference factor level
  target_list$samples$group <- grouping.vec
  
  # Estimate dispersions for multiple group experiment
  target_list <- estimateCommonDisp(y = target_list) # Group is already assigned to the DGEList and should not be supplied here or will produce error. The function automatically uses the held group. 
  target_list <- estimateTagwiseDisp(target_list)
  
  # Build model design
  design <- model.matrix(~0+target_list$samples$group) # Inclusion of 0 allows one to specify the reference when doing multiple contrasts
  fit <- glmQLFit(y = target_list, design = design)
  
  # Rename design columns with usable names (dynamic)
  colnames(design) <- gsub(pattern = "target\\_list\\$samples\\$", replacement = "", x = colnames(design))
  
  # Provide contrasts, since there are more than two groups
  my.contrasts <- makeContrasts(     groupmed-grouplow
                                   , grouphigh-grouplow
                                   , grouphigh-groupmed
                                   , levels=colnames(design)
  )
  
  # Run DE test
  results <- glmQLFTest(glmfit = fit, contrast = my.contrasts)
  all_results <- topTags(object = results
                         , n = length(target_list$genes$transcript.id)
                         , adjust.method = "none"
                         , p.value = 0.001
  )
  
  output.FN <- paste0("04_txomic_results/DE_results_carbonate_", datatype_current, ".txt")
  
  write.table(x = all_results, file = output.FN, quote = F, sep = "\t", row.names = F)
  
  #### DE - phenotype of interest 3: survival ####
  # Set the explanatory variable
  categorical_explan_variable <- "surv.bins"
  
  # Assign this phenotype to the group vector in the DGEList
  grouping.vec <- as.factor(samples_and_phenos.df[, categorical_explan_variable])
  grouping.vec <- relevel(x = grouping.vec, ref = "low") # Assign the reference as the reference factor level
  target_list$samples$group <- grouping.vec
  
  # Estimate dispersions for multiple group experiment
  target_list <- estimateCommonDisp(y = target_list) # Group is already assigned to the DGEList and should not be supplied here or will produce error. The function automatically uses the held group. 
  target_list <- estimateTagwiseDisp(target_list)
  
  # Build model design
  design <- model.matrix(~0+target_list$samples$group) # Inclusion of 0 allows one to specify the reference when doing multiple contrasts
  fit <- glmQLFit(y = target_list, design = design)
  
  # Rename design columns with usable names (dynamic)
  colnames(design) <- gsub(pattern = "target\\_list\\$samples\\$", replacement = "", x = colnames(design))
  
  # Provide contrasts, since there are more than two groups
  my.contrasts <- makeContrasts(     groupmed-grouplow
                                     , grouphigh-grouplow
                                     , grouphigh-groupmed
                                     , levels=colnames(design)
  )
  
  # Run DE test
  results <- glmQLFTest(glmfit = fit, contrast = my.contrasts)
  all_results <- topTags(object = results
                         , n = length(target_list$genes$transcript.id)
                         , adjust.method = "none"
                         , p.value = 0.001
  )
  
  output.FN <- paste0("04_txomic_results/DE_results_surv_", datatype_current, ".txt")
  
  write.table(x = all_results, file = output.FN, quote = F, sep = "\t", row.names = F)
  
  
}


#### 10. Gene Of Interest plotting ####
### Clam garden consistent DEGs
lists_of_interest[["gill"]]$samples # the filtered, normalized data used for DE analysis
names(lists_of_interest)

# Genes to plot: 
goi.list <- list(); gois_prepared.list <- list()
goi.list[["CG"]] <- c("25388346", "25389949")
goi.list[["dig_surv"]] <- c("25356128", "25357396", "25364008")
goi.list[["gill_surv"]] <- c("25371088", "25357391", "25385377")

cpm.data <- NULL; cg.gois <- NULL; surv.gois <- NULL; gois <- NULL; tissue.of.interest <- NULL; DGE.data <- NULL
for(i in 1:length(names(lists_of_interest))){
  
  # Tissue this round
  tissue.of.interest <- names(lists_of_interest)[i]
  
  # Obtain the DGEList of interest
  DGE.data <- lists_of_interest[[tissue.of.interest]]
  
  # obtain log2 cpm
  cpm.data <- cpm(y = DGE.data, normalized.lib.sizes = T, log = TRUE)
  
  # rows = transcripts; cols = samples
  dim(cpm.data)
  colnames(cpm.data) # samples
  
  # Transpose for plotting purposes
  cpm.data <- t(cpm.data)
  dim(cpm.data)
  rownames(cpm.data)
  colnames(cpm.data)[1:5]
  
  # Which GOIs for this tissue? 
  cg.gois <- goi.list[["CG"]]
  surv.gois <- goi.list[[ paste0(names(lists_of_interest)[i],"_surv") ]] # Survival genes
  
  # Obtain the cpm data for only the GOIs
  gois <- cpm.data[, c(cg.gois          # clam garden
                       , surv.gois      # survival 
                       )] 
  
  gois.df <- as.data.frame(gois)
  
  # Obtain matching character from the rownames for upcoming phenotype merge
  gois.df$sample <- rownames(gois.df)
  gois.df$sample <- gsub(pattern = ".eff.counts", replacement = "", x = gois.df$sample)
  gois.df <- separate(data = gois.df, col = "sample", into = c("prefix", "plot_and_suffix")
                      , sep = "CG_", remove = T)
  gois.df <- separate(data = gois.df, col = "plot_and_suffix", into = c("tissue", "plot")
                      , sep = "_", remove = T
  )
  
  
  head(gois.df)
  
  # Merge with phenotypes
  gois.df <- merge(x = gois.df, y = samples_and_phenos.df, by = "plot")
  head(gois.df)
  
  # TODO: although it doesn't technically matter, here we don't technically use tissue-specific pheno data, should be ok
  
  gois_prepared.list[[names(lists_of_interest)[i]]] <- gois.df
  
}


# Plotting GOIs
# Obtain data from: gois_prepared.list

gois_prepared.list[["dig"]]$Type <- factor(x = gois_prepared.list[["dig"]]$Type, levels = c("Ref", "CG"))
gois_prepared.list[["gill"]]$Type <- factor(x = gois_prepared.list[["gill"]]$Type, levels = c("Ref", "CG"))

# Clam garden plot, consistent GOIs: 
pdf(file = "04_txomic_results/GOI_clam_garden.pdf", width = 7, height = 5)
par(mfrow=c(2,2), mar= c(3,3,3,1) + 0.2, mgp = c(2,0.75,0))
boxplot(gois_prepared.list$dig$`25388346` ~ gois_prepared.list$dig$Type, las = 1
        , ylab = "log2 cpm(vwd), digestive gland", xlab = "Beach Type")

boxplot(gois_prepared.list$gill$`25388346` ~ gois_prepared.list$gill$Type, las = 1
        , ylab = "log2 cpm(vwd), gill", xlab = "Beach Type")

boxplot(gois_prepared.list$dig$`25389949` ~ gois_prepared.list$dig$Type, las = 1
        , ylab = "log2 cpm(pear1), digestive gland", xlab = "Beach Type")

boxplot(gois_prepared.list$gill$`25389949` ~ gois_prepared.list$gill$Type, las = 1
        , ylab = "log2 cpm(pear1), gill", xlab = "Beach Type")
dev.off()


# Survival plot, tissue-specific GOIs
## Also consider R-squared value derived from these normalized cpm values shown in the plot

# If need to highlight a single beach...
# Toggle the following to TRUE: 
plot_specific_beach <- FALSE
#plot_specific_beach <- TRUE

# Digestive gland plot cols
plot.col <- rep(x = "black", times = length(gois_prepared.list$dig$beach))
special_tag <- ""

if(plot_specific_beach==TRUE){

  plot.col[which(gois_prepared.list$dig$beach=="D")] <- "red"
  special_tag <- "_colour_by_specific_beach"
  
}

plot.col
special_tag

pdf(file = paste0("04_txomic_results/GOI_expr_survival", special_tag,".pdf"), width = 9, height = 5)
par(mfrow=c(2,3)
    , mar= c(3,3,3,1) + 0.2
    , mgp = c(2,0.75,0)
    )
## Dig
# 25356128
plot(x = gois_prepared.list$dig$surv, y = gois_prepared.list$dig$`25356128`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( hsp70-12a ), Dig."
     , col = plot.col
     )

mod <- lm(gois_prepared.list$dig$`25356128` ~ gois_prepared.list$dig$surv)
out.mod <- summary(mod)
out.mod
text(x = 60, y = -2, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))

# 25357396
plot(x = gois_prepared.list$dig$surv, y = gois_prepared.list$dig$`25357396`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( aminopeptidase ), Dig."
     , col = plot.col
)

mod <- lm(gois_prepared.list$dig$`25357396` ~ gois_prepared.list$dig$surv)
out.mod <- summary(mod)
out.mod
text(x = 85, y = 5, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))

# 25364008
plot(x = gois_prepared.list$dig$surv, y = gois_prepared.list$dig$`25364008`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( mucin-2 ), Dig."
     , col = plot.col
)

mod <- lm(gois_prepared.list$dig$`25364008` ~ gois_prepared.list$dig$surv)
out.mod <- summary(mod)
out.mod
text(x = 85, y = 4, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))


## Gill
plot.col <- rep(x = "black", times = length(gois_prepared.list$gill$beach))

if(plot_specific_beach==TRUE){
  
  plot.col[which(gois_prepared.list$gill$beach=="D")] <- "red"  
  
}

# 25371088
plot(x = gois_prepared.list$gill$surv, y = gois_prepared.list$gill$`25371088`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( hsp70-12a ), Gill"
     , col = plot.col
)

mod <- lm(gois_prepared.list$gill$`25371088` ~ gois_prepared.list$gill$surv)
out.mod <- summary(mod)
out.mod
text(x = 85, y = 2, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))

# 25357391
plot(x = gois_prepared.list$gill$surv, y = gois_prepared.list$gill$`25357391`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( hsp90 ), Gill"
     , col = plot.col
)

mod <- lm(gois_prepared.list$gill$`25357391` ~ gois_prepared.list$gill$surv)
out.mod <- summary(mod)
out.mod
text(x = 85, y = 5, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))

# 25385377
plot(x = gois_prepared.list$gill$surv, y = gois_prepared.list$gill$`25385377`, las = 1
     , xlab = "Survival (%)"
     , ylab = "log2 cpm( toll-like-rec. 1 ), Gill"
     , col = plot.col
)

mod <- lm(gois_prepared.list$gill$`25385377` ~ gois_prepared.list$gill$surv)
out.mod <- summary(mod)
out.mod
text(x = 85, y = 2.4, labels = bquote("Adjusted R"^2*"="*.(round(out.mod$adj.r.squared, digits = 2))))

dev.off()


#### Export doi.summary
doi.summary
capture.output(doi.summary, file = "04_txomic_results/DGEList_summaries.txt")
