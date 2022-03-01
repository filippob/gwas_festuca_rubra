## R script to carry out a GWAS analysis with base GLM
# run as Rscript --vanilla 4b.gwas_glm.R

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
    trait = 'targetTemp',
    npc = 4, ## n. of PCs to include
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("afex")
library("qqman")
library("lmtest")
library("ggfortify")
library("tidyverse")
library("data.table")

print("GWAS using LM for continuous traits")

print(paste("genotype file name:",config$genotype_file))
print(paste("phenotype file name:",config$phenotype_file))
print(paste("trait:",config$trait))
print(paste("number of principal components to include:",config$npc))

dataset = basename(config$genotype_file)

## READING DATA
print("now reading in the data ...")
### genotypes
genotypes <- fread(config$genotype_file, header = TRUE)
snp_matrix = genotypes[,-c(1:8)]
print(paste(nrow(snp_matrix),"records read from the phenotype file",sep=" "))
SNP_INFO <- genotypes[,1:6]
SNP_INFO <- SNP_INFO %>% dplyr::select(c(contig,marker,start_pos))
names(SNP_INFO) <- c("Chr","SNP","Pos")

matg <- t(as.matrix(snp_matrix))
colnames(matg) <- SNP_INFO$SNP
rownames(matg) <- colnames(snp_matrix)

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
vec <- which(colSums(matg)/nrow(matg) < 0.95 & colSums(matg)/nrow(matg) > 0.025)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]

### phenotypes
print("now reading in the binary phenotype ...")
phenotypes <- fread(config$phenotype_file)
phenotypes <- dplyr::select(phenotypes, c(sample, !!as.name(config$trait)))
head(phenotypes)
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% row.names(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

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
phenotypes <- phenotypes %>% dplyr::rename(phenotype = !!as.name(config$trait))
n <- config$npc ## n. of principal components to use for GWAS
phenotypes <- cbind.data.frame(phenotypes,pc$x[,1:n])


###################
## Running the GWAS
###################
writeLines(' - running the GWAS')
df <- phenotypes[,-1] ## remove ID column
res = data.frame("SNP"=NULL, "effect"=NULL,"pvalue"=NULL)

for(i in 1:ncol(matg)) {
  
  if (i%%100 == 0) print(paste("analysing SNP n.",i))
  df$snp <- matg[,i]
  snp_name <- colnames(matg)[i]
  
  fit <- glm(phenotype ~ ., data = df, family = gaussian(link = "identity"))
  
  vv <- car::Anova(fit, type = "III")
  pvalue <- vv["snp","Pr(>Chisq)"]
  snp_effect <- as.numeric(fit$coefficients["snp"])
  
  res = rbind.data.frame(res, data.frame("SNP"=snp_name, "effect"=snp_effect,"pvalue"=pvalue))
}

###########
### RESULTS
###########
print(" - writing out results and figures ...")
res <- res %>% inner_join(SNP_INFO, by = "SNP")
res <- dplyr::select(res,c(SNP,Chr,Pos,effect,pvalue))

fname <- paste(dataset,config$trait,"GWAS_lm.results", sep="_")
fwrite(x = res, file = paste(config$base_folder,"/results/",fname,sep=""))

res$effect <- NULL
names(res) <- c("SNP","CHR","BP","P")

temp <- filter(res, P < 0.05)
pos <- temp %>%
  group_by(CHR) %>%
  summarise(N =n()) %>%
  arrange(desc(N))
temp <- temp %>% inner_join(pos, by = "CHR")
temp <- arrange(temp,desc(N))
temp$CHR <- factor(temp$CHR, levels = unique(temp$CHR))

p <- ggplot(temp, aes(x = CHR, y = -log(P)))  + geom_jitter(aes(color = P), width = 0.5)
# p <- p + facet_wrap(~CHR)
p <- p + theme(
  axis.text.x = element_text(angle = 90), text = element_text(size = 6)
        )
# p

fname = paste(dataset,config$trait,"manhattan_lm.pdf",sep="_")
ggsave(filename = paste(config$base_folder,"/results/",fname,sep=""), plot = p, device = "pdf", width = 9, height = 9)

## qq-plot
fname = paste(dataset,config$trait,"qqplot_lm.png",sep="_")
png(paste(config$base_folder,"/results/",fname,sep=""), width = 600, height = 600)
qq(res$P)
dev.off()

print("#########")
print("## END ##")
print("#########")



