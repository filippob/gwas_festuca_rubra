#!/usr/bin/Rscript
## R script to prepare the file of genotypes for GWAS  

library("tidyverse")
library("data.table")

################################
## parameters
################################
prjfolder = "/home/filippo/Documents/zuzana_festuca_rubra/2nd_project"
fname = "Festuci-sequence-Alpy.xlsx"
phenotype_file = "phenotypes.csv"
sheet_num = 1
exclude_target_locality = NULL
subset_snps = "results/filtered_genotypes.csv_Age_GWAS_multinomial.results" ## file with SNPs to consider for the analysis (subset of the total n. of snps/alleles/markers)
# exclude_target_locality = c("Lavisdalen","Ulvehaugen") ## vector or NULL
keep_shifted_only = TRUE
excl_sample = FALSE
################################
################################

## read and transform
writeLines(" - reading in the genotypic data ...")
unfiltered_genotypes <- readxl::read_xlsx(file.path(prjfolder, fname), sheet = sheet_num, col_names = TRUE, skip = 1)

## get matrix of genotypes (0/1 coded)
M <- unfiltered_genotypes[,-c(1:7)]
print(paste(nrow(M), "unfiltered markers have been read"))

## match between genotypes ans phenotypes
writeLines(" - matching with phenotypes ...")
phenotypes <- fread(file.path(prjfolder, phenotype_file))
vec <- colnames(M) %in% phenotypes$sample
M <- M[,vec]
print(paste(ncol(M), "samples left after matching genotypes and phenotypes"))

## recode from n. of reads mapped to presence/absence
## probably using the n. of reads would introduce a bias: amplification, depth of sequencing, chance of mapping reads to one location etc.
M_recoded = ifelse(M>0,1,0)

## get marker position on contig
writeLines(" - getting the contig of SSR")
unfiltered_genotypes$contig = gsub(":.*$","",unfiltered_genotypes$Locus)

writeLines(" - getting the position of SSR")
unfiltered_genotypes$position = unfiltered_genotypes |>
  separate(Locus, into = c("contig", "start-end", "variant"), sep = ":") |>
  separate(`start-end`, into = c("start", "end"), sep = "-") |>
  mutate(Length = gsub(" bases","",Length), pos = as.integer(start) + as.integer(Length)) |>
  pull(pos)

# unfiltered_genotypes |>
#   group_by(contig) |>
#   summarise(N = n()) |>
#   ggplot(x = N) + geom_histogram(aes(x = N))
# 
# unfiltered_genotypes |>
#   group_by(contig, Length) |>
#   summarise(N = n()) |>
#   pull(N) |>
#   table()

M <- unfiltered_genotypes |>
  select(Allelle, contig,position) |>
  rename(allele = Allelle) |>
  bind_cols(M_recoded)

## take subset of markers
if (subset_snps != "" & !is.null(subset_snps)) {
  
  fname = file.path(prjfolder, subset_snps)
  subset = fread(fname)
  subset <- subset |>
    filter(pvalue <= 0.05) |>
    select(SNP, Chr, Pos) |>
    rename(allele = SNP, contig = Chr, position = Pos)
  
  M <- inner_join(M, subset, by = c("allele", "contig", "position"))
  
}

## data cleaning
writeLines(" - cleaning data: removing alleles with zero count in all samples")
temp <- M[,-c(1:3)]
nsample = ncol(temp)
freqs = rowSums(temp)/nsample
vec <- (freqs < 0.95 & freqs > 0.05)
M <- M[vec,]
print(paste(nrow(M), "markers left after cleaning"))

## write out data
writeLines(" - writing out the data ... ")
fwrite(x = M, file = file.path(prjfolder, "filtered_genotypes.csv"), sep = ",", col.names = TRUE)

print("DONE!")

