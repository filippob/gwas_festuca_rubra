## R script to carry out a GWAS analysis with the package rrBLUP
## kinship matrix used to account for population structure in the data
## input: Plink .raw and .map files + phenotype file
# run as Rscript --vanilla gwas_multinomial.R <config_file>

###################################
## read arguments from config file
###################################
# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
  #loading the parameters
  source(args[1])
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    base_folder = '~/Documents/zuzana_festuca_rubra',
    genotype_file = 'filtered_genotypes.csv',
    phenotype_file = 'phenotypes.csv',
    trait = 'warm_wet',
    npc = 4, ## n. of PCs to include
    plots = TRUE, ## should plots be plotted out
    force_overwrite = FALSE
  ))
}

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("MASS")
library("nnet")
library("qqman")
library("lmtest")
# library("ggfortify")
library("tidyverse")
library("data.table")

print("GWAS using GLM for multinomial traits")

print(paste("genotype file name:",config$genotype_file))
print(paste("phenotype file name:",config$phenotype_file))
print(paste("trait:",config$trait))
print(paste("number of principal components to include:",config$npc))

dataset = basename(config$genotype_file)

## READING DATA
writeLines(" - now reading in the genotypic data ...")
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
writeLines(" - now subsampling the marker loci ...")
vec <- sample(1:ncol(matg),2000)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]
summary(colSums(matg)/nrow(matg))
####

################
## filtering
################
writeLines(" - now filtering the marker data for frequency ...")
vec <- which(colSums(matg)/nrow(matg) < 0.95 & colSums(matg)/nrow(matg) > 0.025)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]
print(paste("n. of markers left:", nrow(SNP_INFO)))

#################
### phenotypes
print("now reading in the binary phenotype ...")
phenotypes <- fread(config$phenotype_file)
phenotypes <- dplyr::select(phenotypes, c(sample, !!as.name(config$trait)))
head(phenotypes)
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% row.names(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## convert score to factor
phenotypes$warm_wet = factor(phenotypes$warm_wet)

#################
## kinship matrix
writeLines(' - reading and preparing the kinship matrix')
K <- fread("kinship_gower.csv", sep=",", header = TRUE)
vec <- colnames(K) %in% phenotypes$sample
K <- K[vec,vec, with=FALSE]

#########################
## Principal Components
#########################
writeLines(' - calculating principal components')
pc <- prcomp(matg)
phenotypes <- phenotypes %>% dplyr::rename(phenotype = !!as.name(config$trait)) %>% mutate(phenotype = as.factor(phenotype))
n <- config$npc ## n. of principal components to use for GWAS
phenotypes <- cbind.data.frame(phenotypes,pc$x[,1:n])


###################
## Running the GWAS
###################
target <- row.names(matg)
phenotypes <- phenotypes[match(target, phenotypes$sample),]
phenotypes$phenotype <- relevel(phenotypes$phenotype, ref = "0")

results = data.frame("SNP"=NULL,"odds_0.1"=NULL,"odds_1.0"=NULL,"odds_1.1"=NULL,
                     "p_0.1"=NULL,"p_1.0"=NULL,"p_1.1"=NULL,"p_overall"=NULL)

for(i in 1:ncol(matg)) {
  
  phenotypes$snp <- matg[,i]
  snp_name <- colnames(matg)[i]
  
  print(paste("SNP n.", i, snp_name))
  
  multinom_model <- multinom(phenotype ~ ., data = dplyr::select(phenotypes,-sample), Hess = TRUE, model = TRUE)
  
  ## overall p-value
  obj = lrtest(multinom_model, "snp")
  p_value = obj$`Pr(>Chisq)`[2]
  
  ## per-class p-values (relative to reference)
  coefs = summary(multinom_model)$coefficients
  std_errs = summary(multinom_model)$standard.errors
  chisq_stats = (coefs^2)/(std_errs^2)
  p <- 1-pchisq(q = chisq_stats, df = 1)
  p_0.1 = p["0.1","snp"]
  p_1.0 = p["1","snp"]
  p_1.1 = p["1.1","snp"]
  
  ## per-class odds (snp effects)
  odds_0.1 = coef(multinom_model)["0.1","snp"]
  odds_1.0 = coef(multinom_model)["1","snp"]
  odds_1.1 = coef(multinom_model)["1.1","snp"]
  
  ## results
  results = rbind.data.frame(results,data.frame(
    "SNP"=snp_name,"odds_0.1"=odds_0.1,"odds_1.0"=odds_1.0,"odds_1.1"=odds_1.1,
    "p_0.1"=p_0.1,"p_1.0"=p_1.0,"p_1.1"=p_1.1,"p_overall"=p_value))
}


###########
### RESULTS
###########
print("writing out results and figures ...")
results <- results %>% inner_join(SNP_INFO, by = c("SNP" = "marker")) %>% relocate(any_of(c("SNP","Chrom","Position")))

fname <- paste(dataset,config$trait,"GWAS_multinomial.results", sep="_")
fwrite(x = results, file = paste(config$base_folder,"/results/",fname,sep=""))


res <- dplyr::select(results, c(SNP,Chrom,Position,p_overall)) %>% rename(BP = Position, CHR = Chrom, P = p_overall)

if(config$plots == TRUE) {
  
  # manhattan plot
  temp <- filter(res, P < 0.05)
  pos <- temp %>%
    group_by(CHR) %>%
    summarise(N =n()) %>%
    arrange(desc(N))
  temp <- temp %>% inner_join(pos, by = "CHR")
  temp <- arrange(temp,desc(N))
  temp$CHR <- factor(temp$CHR, levels = unique(temp$CHR))
  
  p <- ggplot(temp, aes(x = CHR, y = -log(P)))  + geom_jitter(aes(color = P), width = 0.5)
  p <- p + theme(
    axis.text.x = element_text(angle = 90), text = element_text(size = 6)
  )
  
  fname = paste(dataset,config$trait,"manhattan_glm.pdf",sep="_")
  ggsave(filename = paste(config$base_folder,"/results/",fname,sep=""), plot = p, device = "pdf", width = 9, height = 9)

  ## qq-plot
  fname = paste(dataset,config$trait,"qqplot_glm.png",sep="_")
  png(paste(config$base_folder,"/results/",fname,sep=""), width = 600, height = 600)
  qq(res$P)
  dev.off()
}

print("Done!!")

print("#########")
print("## END ##")
print("#########")



