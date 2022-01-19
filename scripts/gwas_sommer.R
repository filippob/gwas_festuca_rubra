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
    trait = 'wetter',
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

### phenotypes
writeLines(" - now reading in the phenotypes ...")
phenotypes <- fread(config$phenotype_file)
# phenotypes <- phenotypes[,c(1,3)]
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% rownames(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## create sample IDs (pure numbers are misinterpreted by sommer: numeric instead of strings, wrong matrix indexing)
phenotypes <- phenotypes %>% mutate(id = paste("s",sample,sep="_"))

## PC's
writeLines(" - calculating principal components of the genotype matrix")
pc <- prcomp(matg)

## add first n principal components to the file of phenotypes
phenotypes <- cbind.data.frame(phenotypes,pc$x[,1:config$npc])

## kinship matrix
writeLines(" - calculating the kinship matrix")
K <- fread("kinship_vr.csv", header = TRUE)
names(K) <- phenotypes$id
K <- as.matrix(K)
rownames(K) <- colnames(K)

## marker matrix
X <- as.matrix(matg)
rownames(X) <- phenotypes$id

## subset phenotypes
P <- phenotypes %>% dplyr::rename(phenotype = !!as.name(config$trait))
covs = paste("PC", seq(1,config$npc),sep="")
P <- dplyr::select(P, c(id, phenotype, all_of(covs)))

###################
## Running the GWAS
###################
writeLines(" - running the GWAS")
## build model
fmod <- as.formula(
  paste("phenotype",
        paste(covs, collapse = "+"),
        sep = " ~ "))

mix_mod <- GWAS(fmod,
                random = ~vs(id, Gu=K),
                rcov = ~units,
                data = P,
                M = X,
                gTerm = "u:id", 
                verbose = TRUE)

###########
### RESULTS
###########
writeLines(" - writing out results and figures ...")
ms <- as.data.frame(mix_mod$scores)
ms$marker <-gsub("_[A-Z]{1}$","",rownames(ms))

# convert -log(p) back to p-values
p <- 10^((-ms$phenotype))

mapf <-merge(SNP_INFO, ms, by="marker", all.x = TRUE);
mapf$pvalue <- p
mapf <- filter(mapf, phenotype < Inf)

fname <- paste(dataset,config$trait,"GWAS_sommer.results", sep="_")
fwrite(x = rename(mapf, log_pval = phenotype), file = paste(config$base_folder,"/results/",fname,sep=""))

## plot
temp <- filter(mapf, pvalue < 0.05)
pos <- temp %>%
  group_by(Chrom) %>%
  summarise(N =n()) %>%
  arrange(desc(N))
temp <- temp %>% inner_join(pos, by = "Chrom")
temp <- arrange(temp,desc(N))
temp$Chrom <- factor(temp$Chrom, levels = unique(temp$Chrom))

p <- ggplot(temp, aes(x = Chrom, y = -log(pvalue)))  + geom_jitter(aes(color = pvalue), width = 0.5)
# p <- p + facet_wrap(~CHR)
p <- p + theme(
  axis.text.x = element_text(angle = 90), text = element_text(size = 6)
)
# p

fname = paste(dataset,config$trait,"manhattan_sommer.pdf",sep="_")
ggsave(filename = paste(config$base_folder,"/results/",fname,sep=""), plot = p, device = "pdf", width = 9, height = 9)

## qq-plot
fname = paste(dataset,config$trait,"qqplot_sommer.png",sep="_")
png(paste(config$base_folder,"/results/",fname,sep=""), width = 600, height = 600)
qqman::qq(mapf$pvalue)
dev.off()

print("#########")
print("## END ##")
print("#########")



