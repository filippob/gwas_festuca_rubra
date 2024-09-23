#!/usr/bin/Rscript
## R script to create the kinship matrix for GWAS  

library("madbito")
library("cluster")
library("data.table")

################################
## parameters
################################
prjfolder = "/home/filippo/Documents/zuzana_festuca_rubra/2nd_project"
genotype_file = "filtered_genotypes.csv"
phenotype_file = "phenotypes.csv"
sheet_num = 1
exclude_target_locality = NULL
# exclude_target_locality = c("Lavisdalen","Ulvehaugen") ## vector or NULL
keep_shifted_only = TRUE
excl_sample = FALSE
################################
################################

### support functions ##########################
modified_van_raden <- function (data)  {
  
  ## matrix of marker alleles in the format (n,m)
  ## coded: -1/0/1
  p = colSums(data)/(nrow(data)) ## calculated as haploid genome (presence/absence data)
  p.scaled = (p - 0.5) ## rescaling frequencies based on haploid genome
  Z = as.matrix((data) - matrix(rep(p.scaled, nrow(data)), 
                                nrow = nrow(data), byrow = TRUE))
  kin = (Z %*% t(Z))/(sum(p * (1 - p))) ## kinship based on haploid genome (removing the coefficient 2)
  return(kin)
}

## function to convert lower triangular matrix to symmetric matrix
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
################################################

## read and transform
geno <- fread(file.path(prjfolder, genotype_file))
M <- geno[,-c(1:3)]


## 1) modified Van Raden
print("- calculating modified Van Raden kinship")
K <- modified_van_raden(t(M))
fwrite(x = K, file = file.path(prjfolder, "kinship_modified_vr.csv"), sep = ",", col.names = TRUE)

## 2) à la Van Raden (2008)
print("- calculating Van Raden kinship")
Kvr = kinship.VR(t(M))
fwrite(x = Kvr, file = file.path(prjfolder, "kinship_vr.csv"), sep = ",", col.names = TRUE)

## 3) Gower similarities
print("- calculating Gower kinship")
dd <- daisy(t(M), metric = "gower")
dd <- makeSymm(as.matrix(dd))
k = (1-dd)
fwrite(x = k, file = file.path(prjfolder, "kinship_gower.csv"), sep = ",", col.names = TRUE)

print("DONE!")

