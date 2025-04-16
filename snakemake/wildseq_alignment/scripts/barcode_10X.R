#!/usr/bin/env Rscript
library(dplyr)
library(argparse)


#parse data
parser <- ArgumentParser(description= 'barcode_10X')
parser$add_argument('--sample', help = 'sample_index')
parser$add_argument('--theoretical_bowtie', help= 'bamtofastq_R2_theoretical.bowtie')
parser$add_argument('--cb_ub_txt', help= 'WILDseq.CB.UB.txt')
parser$add_argument('--output_file', help= 'CBC_table_full.txt')
xargs<- parser$parse_args()


#debugging msgs
if (!file.exists(xargs$theoretical_bowtie)) stop("Theoretical bowtie file not found!")
if (!file.exists(xargs$cb_ub_txt)) stop("WILDseq.CB.UB.txt file not found!")
cat("Processing sample:", xargs$sample, "\n")
cat("Input files loaded successfully.\n")

#script
bowtie <- read.table(xargs$theoretical_bowtie, sep="\t")
colnames(bowtie) <- c("read","strand","barcode_name","offset","seq","qualities","X", "mismatches")

CB.UB <- read.table(xargs$cb_ub_txt)
colnames(CB.UB) <- c("read", "CBC", "UMI")

table <- bowtie %>%
  mutate(read = gsub (" 3:N:0:0", "", as.character(read))) %>%
  inner_join(CB.UB, by = "read") %>%
  select(CBC, UMI, barcode_name) %>%
  distinct()


write.table(table, xargs$output_file, sep="\t", row.names = F, quote = F)


# bowtie$read <- as.character(bowtie$read)
# bowtie$read <- gsub (" 3:N:0:0", "", bowtie$read)
# table <- merge(CB.UB, bowtie, by="read")
# table <- table[,c("CBC", "UMI", "barcode_name")]
# table <- distinct(table)