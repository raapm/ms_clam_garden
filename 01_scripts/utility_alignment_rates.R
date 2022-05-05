align.df <- read.table(file = "../Simple_reads_to_counts/10_log_files/aligned_trimmed_data_to_lsta_txome_log_overall_align_rates.txt")
align.df <- as.data.frame(align.df)

align.df
mean(align.df$V1)
sd(align.df$V1)
