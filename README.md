# ms_clam_garden
A repository used to analyze data used in support of the manuscript 'Effects of ancient anthropogenic clam gardens on the growth, survival, and transcriptome of Pacific littleneck clams (Leukoma staminea)'.      

### Requirements       


All analysis is done within the main directory of the repository.     


### 01. Abiotic factor and growth/ survival analysis
Requires the following input data:      
`02_input_data/cg_sediment_data_2022-03-25.csv`      
`02_input_data/cgsedimentPCA.csv`           

Prepare a PCA of sample locations based on site metadata and survival/ growth, and investigate statistical correlates of survival and growth by running interactively:        
`01_scripts/01_pheno_and_abiotic_var.R`       


### 02. Transcriptomic analysis
#### 02.a Raw reads to counts
To add.     

#### 02.b Analysis of transcriptome data
Use the following interactive script:      
`01_scripts/02_transcriptomic_analysis.R`     

Requires the following input data
`02_input_data/out.matrix.csv` (gene expression data)       
`02_input_data/cg_sediment_data_2022-03-25.csv` (phenotype data)      

This script will allow you to analyze all data together ('all'), or individually by tissue ('gill' or 'dig' for digestive gland.       

Set the variables:         
Read mapping cutoff: `min.reads.mapping.per.transcript` (default = 10)      
Minimum individuals passing threshold: `min.ind`       

Go to the section `Exploratory MDS plots` to check custom labels on the mds plot. This will require choosing a datatype (all, gill, dig) from a previously-generated list.       

Note: original (v.0.1) analysis remains in section: 'Method 2. Original Files'





