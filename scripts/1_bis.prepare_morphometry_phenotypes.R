#!/usr/bin/Rscript
## R script to prepare the file of phenotypes for GWAS  

library("tidyverse")
library("data.table")

## parameters
exclude_target_locality = NULL


## read and transform
print("reading in the phenotypic data ...")
fname = "paper/gen-trait1.xlsx"
phenos = readxl::read_xlsx(fname)
names(phenos) <- gsub("_","-",names(phenos))

print(paste("N. phenotypic records", nrow(phenos)))

print("reading in the genotype data ...")
fname = "filtered_genotypes.csv"
genos = fread(fname)

print("filter phenotype records that have a corresponding genotype record")
phenos <- phenos |>
  filter(sample %in% as.numeric(names(genos)[-c(1:8)]))

print(paste("N. phenotypic records", nrow(phenos)))

print("remove records where all trait values are missing")
phenos <- phenos |>
  mutate(across(everything(), as.numeric))

ntraits = ncol(phenos)-1
vec <- rowSums(is.na(phenos)) < ntraits
phenos <- phenos[vec,]

print(paste("N. phenotypic records", nrow(phenos)))

## write out file
print("writing out the file ...")
fwrite(x = phenos, file = "phenotypes_morphometry.csv", sep = ",", col.names = TRUE)

print("DONE!")
