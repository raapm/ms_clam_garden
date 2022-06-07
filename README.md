# ms_clam_garden
Analyze data for the manuscript 'Effects of ancient anthropogenic clam gardens on the growth, survival, and transcriptome of Pacific littleneck clams (Leukoma staminea)'.      
All analysis is done within the main directory of the repository.     

### Requirements       
See requirements within [Simple_reads_to_counts](https://github.com/bensutherland/Simple_reads_to_counts) for paired-end reads against a _de novo_ transcriptome        
See requirements within the R scripts for R libraries.      

Requires input data obtained from FigShare: DOI: 10.6084/m9.figshare.19735753         

### 01. Abiotic factor and growth/ survival analysis
Requires the following input data (obtained from FigShare, see above):      
`02_input_data/cg_sediment_data_2022-03-25.csv`      

Use the script interactively:      
`01_scripts/01_pheno_and_abiotic_var.R`       

This will allow you to: 
- load data
- plot a sample-based PCA from abiotic variables, with beaches overlaid
- statistical analysis of the effects of abiotic variables on survival and growth
- statistical analysis of the effects of clam garden state on abiotic variables
- plot correlations between variables

Outputs will be in `03_pheno_results`.      


### 02. Transcriptomic analysis
#### 02.a Raw reads to counts
Requires raw read data from SRA (BioProject PRJNA818991, BioSamples SAMN26893882-SAMN26893909).      

Follow the methods using paired-end reads against _de novo_ transcriptome from [Simple_reads_to_counts](https://github.com/bensutherland/Simple_reads_to_counts)        

Take the output from the repo and put it in the current repo, as per:       
`02_input_data/out.matrix_cg_2022-06-07.csv` (gene expression data)        
note: this matrix is also available from FigShare (see DOI above).      

#### 02.b Analysis of transcriptome data
Requires the following input data (obtained from FigShare, see above):      
Sediment/ phenotypic data: `02_input_data/cg_sediment_phenos_2022-05-12.csv`      
Read counts from (02.a):   `02_input_data/out.matrix_cg_2022-06-07.csv`     
Uniprot IDs:               `02_input_data/project155.uniprot_blastp.txt.gz`      
Other annotation:          `02_input_data/cgrnaseqBTv1.csv`      

Use the script interactively:      
`01_scripts/02_transcriptomic_analysis.R`     
...and set the following variables:      
- `min.reads.mapping.per.transcript`, the read mapping cutoff (default = 10)      
- `min.ind`, the minimum individuals passing threshold (default = 5)        


This will allow you to: 
- load data and format
- prepare three datasets for analysis (i.e., gill, digestive gland, all)
- filter for low expression and normalize, then generate MDS plot
- evaluate tissue-specific expression and generate heatmap
- identify differentially expressed genes by Clam Garden or by survival (%)
- plot genes of interest

Outputs will be in `04_txomic_results`.      

**The manuscript is available here**: [preprint](https://www.biorxiv.org/search/ben%252Bj%252Bg%252Bsutherland)
