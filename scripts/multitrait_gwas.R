
### R script for bivariate GWAS models
### using the sommer package

# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1){
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
    traits = 'targetMois,targetTemp',
    npc = 4, ## n. of PCs to include,
    use_kinship = FALSE,
    force_overwrite = FALSE
  ))
  
}

# SETUP -------------------------------------------------------------------
library("sommer")
library("polycor")
library("ggplot2")
library("tidyverse")
library("data.table")

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

print("GWAS using the sommer package")

print(paste("genotype file name:",config$genotype_file))
print(paste("phenotype file name:",config$phenotype_file))
print(paste("trait:",config$trait))
print(paste("number of principal components to include:",config$npc))
print(paste("Use kinship:",config$use_kinship))

dataset = basename(config$genotype_file)

## READING DATA
writeLines(" - now reading in the genotypic data ...")
### genotypes
genotypes <- fread(config$genotype_file, header = TRUE)
print(paste(nrow(genotypes),"records read from the genotype file",sep=" "))
SNP_INFO <- genotypes %>%
  dplyr::select(c(marker,contig,start_pos,end_pos)) %>%
  mutate(Position = (start_pos+end_pos)/2) %>%
  dplyr::select(-c(start_pos,end_pos)) %>%
  rename(Chrom = contig)

temp <- genotypes[,-c(1:8,ncol(genotypes)), with = FALSE]
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

### phenotypes
phenotypes <- fread(config$phenotype_file)
# phenotypes <- phenotypes[,c(1,3)]
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% rownames(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## create sample IDs (pure numbers are misinterpreted by sommer: numeric instead of strings, wrong matrix indexing)
phenotypes <- phenotypes %>% mutate(id = paste("s",sample,sep="_"))

## correlation between target temperature and humidity
p <- ggplot(phenotypes, aes(x = targetTemp, y = targetMois)) + geom_jitter()
print(p)

print(paste("Pearson correlation is:", cor(phenotypes$targetTemp, phenotypes$targetMois, method = "pearson")))
print(paste("Spearman correlation is:", cor(phenotypes$targetTemp, phenotypes$targetMois, method = "spearman")))
print(paste("Kendall correlation is:", cor(phenotypes$targetTemp, phenotypes$targetMois, method = "kendall")))
## polychoric correlation
print(paste("Polychoric correlation (ordinal categorical variables) is:", polychor(phenotypes$targetTemp, phenotypes$targetMois)))

## PC's
writeLines(" - calculating principal components of the genotype matrix")
pc <- prcomp(matg)

## add first n principal components to the file of phenotypes
phenotypes <- cbind.data.frame(phenotypes,pc$x[,1:config$npc])

## marker matrix
X <- as.matrix(matg)
rownames(X) <- phenotypes$id

## subset phenotypes
# P <- phenotypes %>% dplyr::rename(phenotype = !!as.name(config$trait))
traits = strsplit(config$traits, split = ",")[[1]]
covs = paste("PC", seq(1,config$npc),sep="")
P <- dplyr::select(phenotypes, c(id, all_of(traits), all_of(covs)))

###################
## Running the GWAS
###################
writeLines(" - running the GWAS")
## build model (cbind(trait1, trait2))
fmod <- as.formula(
  paste(paste("cbind(",config$traits,")",sep=""),
        paste(covs, collapse = "+"),
        sep = " ~ "))

if (config$use_kinship) {
  print("using kinship matrix in the GWAS model")
  
  ## kinship matrix
  writeLines(" - calculating the kinship matrix")
  K <- fread("kinship_vr.csv", header = TRUE)
  vec <- colnames(K) %in% phenotypes$sample
  K <- K[vec,vec, with=FALSE]
  names(K) <- phenotypes$id
  K <- as.matrix(K)
  rownames(K) <- colnames(K)
  
  mix_mod <- GWAS(fmod,
                  # random = ~vs(id, Gu=K),
                  random = ~vs(id, Gtc=unsm(2)),
                  rcov = ~units,
                  data = P,
                  M = X,
                  gTerm = "u:id", 
                  verbose = TRUE)
  
} else {
  print("NOT using kinship matrix in the GWAS model")
  mix_mod <- GWAS(fmod,
                  random = ~vs(id, Gtc=unsm(2)),
                  rcov = ~vs(units, Gtc=diag(2)),
                  data = P,
                  M = X,
                  gTerm = "u:id", 
                  verbose = TRUE) 
}

###########
### RESULTS
###########
writeLines(" - writing out results and figures ...")
ms <- as.data.frame(mix_mod$scores)
ms$marker <-gsub("_[A-Z]{1}$","",rownames(ms))

# convert -log(p) back to p-values
p_temp <- 10^((-ms$targetTemp))
p_mois <- 10^((-ms$targetMois))

mapf <-merge(SNP_INFO, ms, by="marker", all.x = TRUE);
mapf$p_temp <- p_temp
mapf$p_mois <- p_mois
mapf <- filter(mapf, targetTemp < Inf, targetMois < Inf)

fname <- paste(dataset,config$trait,"GWAS_sommer_multitrait.results", sep="_")
fwrite(x = rename(mapf, log_p_temp = targetTemp, log_p_mois = targetMois), file = paste(config$base_folder,"/results/",fname,sep=""))

## plot
temp <- filter(mapf, p_mois < 0.05)
pos <- temp %>%
  group_by(Chrom) %>%
  summarise(N =n()) %>%
  arrange(desc(N))
temp <- temp %>% inner_join(pos, by = "Chrom")
temp <- arrange(temp,desc(N))
temp$Chrom <- factor(temp$Chrom, levels = unique(temp$Chrom))

p <- ggplot(temp, aes(x = Chrom, y = -log(p_temp)))  + geom_jitter(aes(color = p_mois), width = 0.5)
# p <- p + facet_wrap(~CHR)
p <- p + theme(
  axis.text.x = element_text(angle = 90), text = element_text(size = 6)
)
# p

fname = paste(dataset,traits[1],"manhattan_sommer.pdf",sep="_")
ggsave(filename = paste(config$base_folder,"/results/",fname,sep=""), plot = p, device = "pdf", width = 9, height = 9)

## qq-plot
fname = paste(dataset,traits[1],"qqplot_sommer.png",sep="_")
png(paste(config$base_folder,"/results/",fname,sep=""), width = 600, height = 600)
qqman::qq(mapf$p_mois)
dev.off()

print("#########")
print("## END ##")
print("#########")



