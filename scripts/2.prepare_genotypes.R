#!/usr/bin/Rscript
## R script to prepare the file of genotypes for GWAS  

library("tidyverse")
library("data.table")

## parameters
phenotype_file = "phenotypes_morphometry.csv"

## read and transform
writeLines(" - reading in the genotypic data ...")
unfiltered_genotypes <- readxl::read_xlsx("unfiltered.xlsx", sheet = 1, col_names = TRUE)

## get matrix of genotypes (0/1 coded)
M <- unfiltered_genotypes[,-c(1:12)]
print(paste(nrow(M), "unfiltered markers have been read"))

## sanitize column (sample) names so that they match with the IDs in the file of phenotypes
colnames(M) <- gsub("_.*$","",colnames(M))

## match between genotypes ans phenotypes
writeLines(" - matching with phenotypes ...")
phenotypes <- fread(phenotype_file)
vec <- colnames(M) %in% phenotypes$sample
M <- M[,vec]
print(paste(ncol(M), "samples left after matching genotypes and phenotypes"))

## get marker position on contig
writeLines(" - getting the position of SSR - start position")
M <- unfiltered_genotypes[,c(1:4,7:10)] %>% bind_cols(M)
M <- M %>% 
  separate(`....3`, into = c("start_pos","end_pos"), sep = "-") %>% 
  rename(contig = `LOCUS...2`, sequence = `....5`, length = lenght) %>% 
  select(-c(`LOCUS...1`,`LOC-MOT`))

## data cleaning
writeLines(" - cleaning data: removing NODE contig and SSR with missing position")
M <- M %>% filter(contig != "NODE", !is.na(start_pos))
print(paste(nrow(M), "markers left after cleaning"))

## filtering for frequency
vec <- rowSums(M[,-c(1:7)]) > 3
M <- M[vec,]
print(paste(nrow(M),"markers left after filtering for frequency"))

## make marker name
M <- mutate(M, marker = paste(contig,start_pos,end_pos,sequence, sep = "_")) %>% relocate(marker)

## write out data
writeLines(" - writing out the data ... ")
fwrite(x = M, file = "filtered_genotypes.csv", sep = ",", col.names = TRUE)

print("DONE!")

