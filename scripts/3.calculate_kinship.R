#!/usr/bin/Rscript
## R script to create the kinship matrix for GWAS  

library("madbito")
library("cluster")
library("data.table")

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
geno <- fread("filtered_genotypes.csv")
M <- geno[,-c(1:7)]


## 1) modified Van Raden
print("- calculating modified Van Raden kinship")
K <- modified_van_raden(t(M))
fwrite(x = K, file = "kinship_modified_vr.csv", sep = ",", col.names = TRUE)

## 2) à la Van Raden (2008)
print("- calculating Van Raden kinship")
Kvr = kinship.VR(t(M))
fwrite(x = Kvr, file = "kinship_vr.csv", sep = ",", col.names = TRUE)

## 3) Gower similarities
print("- calculating Gower kinship")
dd <- daisy(t(M), metric = "gower")
dd <- makeSymm(as.matrix(dd))
k = (1-dd)
fwrite(x = k, file = "kinship_gower.csv", sep = ",", col.names = TRUE)

print("DONE!")

