# ms_clam_garden
A repository used to analyze data used in support of the manuscript 'Effects of ancient anthropogenic clam gardens on the growth, survival, and transcriptome of Pacific littleneck clams (Leukoma staminea)'.      

### Requirements       

See requirements within the R scripts for R libraries.      
All analysis is done within the main directory of the repository.     


### 01. Abiotic factor and growth/ survival analysis
Requires the following input data obtained from FigShare(#TODO):      
`02_input_data/cg_sediment_data_2022-03-25.csv`      
`02_input_data/cgsedimentPCA.csv`           

Use the script interactively:      
`01_scripts/01_pheno_and_abiotic_var.R`       

Prepare a PCA of sample locations based on site metadata and survival/ growth, and investigate statistical correlates of survival and growth by running interactively:        


### 02. Transcriptomic analysis
#### 02.a Raw reads to counts
Follow the methods using paired-end reads against _de novo_ transcriptome from [Simple_reads_to_counts](https://github.com/bensutherland/Simple_reads_to_counts)        

Take the output from the repo and put it in the current repo, as per:       
`02_input_data/out.matrix.csv` (gene expression data)        


#### 02.b Analysis of transcriptome data
Use the following interactive script:      
`01_scripts/02_transcriptomic_analysis.R`     

Note: this will also use `02_input_data/cg_sediment_data_2022-03-25.csv`      

This script will allow you to analyze all data together ('all'), or individually by tissue ('gill' or 'dig' for digestive gland.       

Set the variables:         
Read mapping cutoff: `min.reads.mapping.per.transcript` (default = 10)      
Minimum individuals passing threshold: `min.ind`       

Go to the section `Exploratory MDS plots` to check custom labels on the mds plot. This will require choosing a datatype (all, gill, dig) from a previously-generated list.       

This script will also let you find the gill- or digestive-gland-specific expressed genes by comparing those that passed filters in one tissue to the other tissue, and visa-versa.       



