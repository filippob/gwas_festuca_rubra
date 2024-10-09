## R script to carry out a GWAS analysis with base GLM
# run as Rscript --vanilla gwas_glm.R

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
    base_folder = '~/Documents/zuzana_festuca_rubra/2nd_project',
    genotype_file = 'filtered_genotypes.csv',
    phenotype_file = 'phenotypes.csv',
    trait = 'Age',
    kinship_file = 'kinship_gower.csv',
    npc = 4, ## n. of PCs to include
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("afex")
library("qqman")
library("lmtest")
library("sommer")
library("ggfortify")
library("tidyverse")
library("data.table")

print("GWAS using LM for continuous traits")

print(paste("genotype file name:",config$genotype_file))
print(paste("phenotype file name:",config$phenotype_file))
print(paste("trait:",config$trait))
print(paste("kinship matrix:",config$kinship_file))
print(paste("number of principal components to include:",config$npc))

dataset = basename(config$genotype_file)

## READING DATA
print("now reading in the data ...")
### genotypes
genotypes <- fread(file.path(config$base_folder, config$genotype_file), header = TRUE)
snp_matrix = genotypes[,-c(1:3)]
print(paste(nrow(snp_matrix),"records read from the phenotype file",sep=" "))
SNP_INFO <- genotypes[,1:3]
names(SNP_INFO) <- c("SNP","Chr","Pos")

matg <- t(as.matrix(snp_matrix))
colnames(matg) <- SNP_INFO$SNP
rownames(matg) <- colnames(snp_matrix)

somme = colSums(matg)
min_count = round(nrow(matg)/20, digits = 0)
vec = (somme > min_count)
matg = matg[,vec]
SNP_INFO <- SNP_INFO[vec,]

print(paste(nrow(SNP_INFO),"SNPs read from the map file",sep=" "))

if ((ncol(matg)) != nrow(SNP_INFO)) {
  
  stop("!! N. of SNPs in the map file not equal to the number of genotyped SNPs in the genotype file")
  
} else print("N. of SNPs in the map and genotype files is the same: this is correct!!")

################
## subsampling
################
# vec <- sample(1:ncol(matg),2000)
# matg <- matg[,vec]
# SNP_INFO <- SNP_INFO[vec,]
# summary(colSums(matg)/nrow(matg))
####

################
## filtering
################
# vec <- which(colSums(matg)/nrow(matg) < 0.95 & colSums(matg)/nrow(matg) > 0.025)
# matg <- matg[,vec]
# SNP_INFO <- SNP_INFO[vec,]

### phenotypes
print("now reading in the binary phenotype ...")
phenotypes <- fread(file.path(config$base_folder, config$phenotype_file))
phenotypes <- dplyr::select(phenotypes, c(sample, Block, !!as.name(config$trait)))
head(phenotypes)
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

## remove missing phenotype records
phenotypes <- na.omit(phenotypes)
phenotypes <- phenotypes[phenotypes$sample %in% row.names(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## kinship matrix
if(config$kinship_file == "") {
  
  print("Calculating the kinship matrix")
  K <-A.mat(matg)
  
  vec <- colnames(K) %in% phenotypes$sample
  K <- K[vec,vec]
  
} else {
  
  writeLines(' - reading and preparing the kinship matrix')
  K <- fread(file.path(config$base_folder, config$kinship_file), sep=",", header = TRUE)
  
  vec <- colnames(K) %in% phenotypes$sample
  K <- K[vec,vec,with=FALSE]
  K <- as.matrix(K)
}

writeLines(" - heatmap of the kinship matrix")
fname = file.path(config$base_folder, "kinship.png")
png(filename = fname, width = 800, height = 800)
heatmap(K)
dev.off()


#########################
## Principal Components
#########################
writeLines(' - calculating principal components')
pc <- prcomp(matg)
n <- config$npc ## n. of principal components to use for GWAS
phenotypes <- cbind.data.frame(phenotypes,pc$x[,1:n])

###################
## Running the GWAS
###################
writeLines(' - running the GWAS')
df <- phenotypes[,-1] ## remove ID column
columns = names(df)
vec <- columns != "Age"
covariates = paste(columns[vec], collapse = ",")

fmod <- as.formula(
  paste("Age",
        gsub(",","+",covariates),
        sep = " ~ "))

phenotypes <- phenotypes |> rename(id = sample)
X <- as.matrix(matg)
rownames(K) <- phenotypes$id
summary(phenotypes)

mix_mod <- sommer::GWAS(fmod,
                random = ~vsr(id, Gu=K),
                rcov = ~units,
                data = phenotypes,
                M = X,
                gTerm = "u:id", 
                verbose = TRUE)

# save(mix_mod, file="res_sommer.RData")
# load("res_sommer.RData")

###########
### RESULTS
###########
print("writing out results and figures ...")
ms <- as.data.frame(mix_mod$scores)
ms$SNP <-gsub("_[A-Z]{1}$","",rownames(ms))
names(ms)[1] <- "phenotype"

# convert -log(p) back to p-values
p <- 10^((-ms$phenotype))

SNP_INFO <- mutate(SNP_INFO, SNP = as.character(SNP))

mapf <-merge(SNP_INFO, ms, by="SNP", all.x = TRUE);
mapf$pvalue <- p
mapf <- filter(mapf, phenotype < Inf)
mapf <- mapf |> rename(Chrom = Chr, Position = Pos, p.val = pvalue)

fname = paste(dataset,config$trait,"manhattan_sommer.png",sep="_")
png(fname)
sommer::manhattan(mapf, pch=20,cex=.75, PVCN = "phenotype")
dev.off()

## rename P to log_p (as it is) and add the column with p-values
names(mapf)[4] <- "log_p"
fname = paste(dataset,config$trait,"GWAS_sommer.results",sep="_")
fwrite(x = mapf, file = fname)

## qq-plot
fname = paste(dataset,config$trait,"qqplot_sommer.png",sep="_")
png(fname, width = 600, height = 600)
qqman::qq(mapf$p.val)
dev.off()


print("#########")
print("## END ##")
print("#########")



