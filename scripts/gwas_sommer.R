## R script to carry out a GWAS analysis with the package rrBLUP
## kinship matrix used to account for population structure in the data
## input: Plink .raw and .map files + phenotype file
# run as Rscript --vanilla gwas_sommer.R genotype_file=path_to_genotypes snp_map=path_to_map phenotype_file=path_to_phenotypes trait=trait_name_in_phenotype_file trait_label=label_to_use_for_trait

# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
  #loading the parameters
  source(args[1])
}else{
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/zuzana_festuca_rubra',
    genotype_file = 'filtered_genotypes.csv',
    phenotype_file = 'phenotypes.csv',
    trait = 'warmer',
    npc = 4, ## n. of PCs to include
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("qqman")
library("sommer")
library("tidyverse")
library("data.table")

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

print("GWAS using the sommer package")

print(paste("genotype file name:",config$genotype_file))
print(paste("phenotype file name:",config$phenotype_file))
print(paste("trait:",config$trait))
print(paste("number of principal components to include:",config$npc))

dataset = basename(config$genotype_file)

## READING DATA
print("now reading in the data ...")
### genotypes
genotypes <- fread(config$genotype_file, header = TRUE)
print(paste(nrow(genotypes),"records read from the genotype file",sep=" "))
genotypes <- mutate(genotypes, marker_name = paste(contig,start_pos,end_pos,sequence, sep = "_"))
SNP_INFO <- genotypes %>%
  dplyr::select(c(marker_name,contig,start_pos,end_pos)) %>%
  mutate(Position = (start_pos+end_pos)/2) %>%
  dplyr::select(-c(start_pos,end_pos)) %>%
  rename(marker = marker_name, Chrom = contig)

temp <- genotypes[,-c(1:7,ncol(genotypes)), with = FALSE]
matg <- t(as.matrix(temp))
colnames(matg) <- SNP_INFO$marker
rownames(matg) <- colnames(temp)
rm(temp)

################
## subsampling
################
vec <- sample(1:ncol(matg),2000)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]
summary(colSums(matg)/nrow(matg))
####

################
## filtering
################
vec <- which(colSums(matg)/nrow(matg) < 0.95 & colSums(matg)/nrow(matg) > 0.025)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]

### phenotypes
phenotypes <- fread(config$phenotype_file)
# phenotypes <- phenotypes[,c(1,3)]
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$id %in% snp_matrix$IID,]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## PC's
genotypes <- fread("filtered_genotypes.csv", header = TRUE)
snp_matrix = genotypes[,-c(1:7)]
matg <- t(as.matrix(snp_matrix))
pc <- prcomp(matg)
names(pc)
pc$x[1:20,1:10]


# SNP_INFO <- bind_cols(SNP_INFO,as.data.frame(t(X)))

print(paste(nrow(SNP_INFO),"SNPs read from the map file",sep=" "))

if ((ncol(snp_matrix)-6) != nrow(SNP_INFO)) {
  
  stop("!! N. of SNPs in the map file not equal to the number of genotyped SNPs in the genotype file")
  
} else print("N. of SNPs in the map and genotype files is the same: this is correct!!")



## kinship matrix
print("Calculating the kinship matrix")
K <-A.mat(X)

vec <- colnames(K) %in% phenotypes$id
K <- K[vec,vec]

# SNP_INFO <- as.data.frame(SNP_INFO)
# SNP_INFO <- SNP_INFO[,c(TRUE,TRUE,TRUE,vec)]

print("producing the heatmap kinship matrix ...")
pdf(paste(dataset,"_kinship_heatmap",".pdf",sep=""))
heatmap(K)
dev.off()

###################
## Running the GWAS
###################
phenotypes <- phenotypes %>% dplyr::rename(phenotype = !!as.name(trait))

fmod <- as.formula(
  paste("phenotype",
        gsub(",","+",covariates),
        sep = " ~ "))

mix_mod <- GWAS(fmod,
                random = ~vs(id, Gu=K),
                rcov = ~units,
                data = phenotypes,
                M = X,
                gTerm = "u:id", 
                verbose = TRUE)

###########
### RESULTS
###########
print("writing out results and figures ...")
ms <- as.data.frame(mix_mod$scores)
ms$snp <-gsub("_[A-Z]{1}$","",rownames(ms))

# convert -log(p) back to p-values
p <- 10^((-ms$phenotype))

mapf <-merge(SNP_INFO, ms, by="snp", all.x = TRUE);
mapf$pvalue <- p
mapf <- filter(mapf, phenotype < Inf)
png(paste(dataset,trait_label,"manhattan_sommer.png",sep="_"))
sommer::manhattan(mapf, pch=20,cex=.75, PVCN = "phenotype")
dev.off()

## rename P to log_p (as it is) and add the column with p-values
names(mapf)[4] <- "log_p"
fname <- paste(dataset,trait_label,"GWAS_sommer.results", sep="_")
fwrite(x = mapf, file = fname)

## qq-plot
png(paste(dataset,trait_label,"qqplot_sommer.png",sep="_"), width = 600, height = 600)
qqman::qq(mapf$pvalue)
dev.off()

print("#########")
print("## END ##")
print("#########")



