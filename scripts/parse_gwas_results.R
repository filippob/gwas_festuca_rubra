## parse GWAS results

library("tidyverse")
library("data.table")

significance_threshold = 1e-3

## POWER
source("~/Dropbox/cursos/gwas_berlin2020/introduction_to_gwas/5.power_and_significance/power_calc_functions.R")
# power_n_hsq(n = seq(1.7,2,0.1)*1000, qsq = seq(1,15,1)/100, pval=5e-5)
pow = power_n_hsq(n = seq(150,250,10), qsq = seq(1,15,1)/100, pval=significance_threshold)
print(pow)
pow["200",]
# power_plot(pow, "n", "q-squared")

library("ldDesign")
## Power calculation Tool
marker.h2=0.08
marker.p=0.5
qtl.p=0.5
sample_size=200
LD.D=0.25
Bayes.Factor=20
ld.power(B=Bayes.Factor, D=LD.D, p=marker.p, q=qtl.p,n=sample_size, h= marker.h2, phi=0)

## Bayesian odds
n_qtl = 25
n_snp = 16039
prior = n_qtl/n_snp # if 10 true associations are expected before the GWAS
power = 0.7 # (based on 500 samples, 0.03 heritability)
p_value = significance_threshold ## threshold

posterior_odds = prior*power/p_value ## P(true_assoc)/(1-P(true_assoc))

# top_pvals <- results %>% arrange(pvalue) %>% head() %>% select(pvalue) %>% pull()
# (prior*power)/top_pvals

### get significant results ###
### read pvalues
# warmer & wetter
results <- fread("results/filtered_genotypes.csv_warm_wet_GWAS_multinomial.results")
results <- select(results, c(SNP,Chrom,Position,p_overall)) %>% rename(pvalue = p_overall)
filtered_results <- filter(results, pvalue < significance_threshold)
filtered_results$posterior_odds = (prior*power)/filtered_results$pvalue
filtered_results$model = "warmer_wetter"

# warmer
results <- fread("results/filtered_genotypes.csv_warmer_GWAS_glm.results")
temp <- filter(results, pvalue < significance_threshold)
temp <- select(temp, -effect) %>% rename(Chrom = Chr, Position = Pos)
temp$posterior_odds = (prior*power)/temp$pvalue
temp$model = "warmer"

filtered_results <- bind_rows(filtered_results,temp)

# wetter
results <- fread("results/filtered_genotypes.csv_wetter_GWAS_glm.results")
temp <- filter(results, pvalue < significance_threshold)
temp <- select(temp, -effect) %>% rename(Chrom = Chr, Position = Pos)
temp$posterior_odds = (prior*power)/temp$pvalue
temp$model = "wetter"

filtered_results <- bind_rows(filtered_results,temp)

library("ggVennDiagram")

warmer <- filter(filtered_results, model == "warmer") %>% pull(SNP)
wetter <- filter(filtered_results, model == "wetter") %>% pull(SNP)
warmer_wetter <- filter(filtered_results, model == "warmer_wetter") %>% pull(SNP)

snps <- list("warmer" = warmer, "wetter" = wetter, "warmer_wetter" = warmer_wetter)

p <- ggplot(filtered_results, aes( x = Position, y = -log(pvalue))) + geom_jitter(aes(color = model), width = 0.5, size=2)
p <- p + facet_wrap(~Chrom)
p <- p + theme(axis.text.x = element_text(angle = 90))
# p

ggsave(filename = "results/significant_snps.png", plot = p, device = "png", width = 10, height = 7)

dd <-filtered_results %>%
  group_by(Chrom, model) %>%
  summarise(pvalue, `-log10` = -log(pvalue), posterior_odds) %>%
  fwrite("results/table_significant_ssr.csv", sep = ",")

fwrite(x = filtered_results, file = "results/significant_ssr_extended.csv", sep=",")


####################################
### target temperature and humidity
####################################

significance_threshold = 1e-3

## Bayesian odds
n_qtl = 25
n_snp = 16039
prior = n_qtl/n_snp # if 10 true associations are expected before the GWAS
power = 0.7 # (based on 500 samples, 0.03 heritability)
p_value = significance_threshold ## threshold

posterior_odds = prior*power/p_value ## P(true_assoc)/(1-P(true_assoc))

### get significant results ###
### read pvalues
# target_temp
results <- fread("results/filtered_genotypes.csv_targetTemp_GWAS_glm.results")
filtered_results <- filter(results, pvalue < significance_threshold) %>% rename(Chrom = Chr, Position = Pos)
filtered_results$posterior_odds = (prior*power)/filtered_results$pvalue
filtered_results$model = "temperature"

# target_mois
results <- fread("results/filtered_genotypes.csv_targetMois_GWAS_lm.results")
temp <- filter(results, pvalue < significance_threshold) %>% rename(Chrom = Chr, Position = Pos)
# temp <- select(temp, -effect) %>% rename(Chrom = Chr, Position = Pos)
temp$posterior_odds = (prior*power)/temp$pvalue
temp$model = "humidity"

filtered_results <- bind_rows(filtered_results,temp)

temperature <- filter(filtered_results, model == "temperature") %>% pull(SNP)
humidity <- filter(filtered_results, model == "humidity") %>% pull(SNP)

snps <- list("temperature" = temperature, "humidity" = humidity)

temp <- group_by(filtered_results, Chrom) %>% summarise(N = n())
filtered_results$N <- temp$N[match(filtered_results$Chrom,temp$Chrom)]

p <- ggplot(filter(filtered_results, N>1), aes( x = Position, y = -log(pvalue))) + geom_jitter(aes(color = model), width = 0.5, size=2)
p <- p + facet_wrap(~Chrom)
p <- p + theme(axis.text.x = element_text(angle = 90))
# p

ggsave(filename = "results/significant_snps_temp_humidity.png", plot = p, device = "png", width = 10, height = 7)

dd <-filtered_results %>%
  group_by(Chrom, model) %>%
  summarise(pvalue, `-log10` = -log(pvalue), posterior_odds) %>%
  fwrite("results/table_significant_ssr_temp_mois.csv", sep = ",")

fwrite(x = filtered_results, file = "results/significant_ssr_extended_temp_mois.csv", sep=",")

