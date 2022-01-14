library("qqman")
library("ggplot2")
library("tidyverse")
library("data.table")

base_folder = "~/Documents/zuzana_festuca_rubra"
trait = "warmer"
model = "glm"

res <- fread("results/filtered_genotypes.csv_warmer_GWAS_glm.results")

temp <- filter(res, pvalue < 0.99)
pos <- temp %>%
  group_by(Chr) %>%
  summarise(N =n()) %>%
  arrange(desc(N))
temp <- temp %>% inner_join(pos, by = "Chr")
temp <- arrange(temp,desc(N))
temp$Chr <- factor(temp$Chr, levels = unique(temp$Chr))

p <- ggplot(temp, aes(x = Chr, y = -log(pvalue)))  + geom_jitter(aes(color = Chr), width = 0.5)
p <- p + theme(
  axis.text.x = element_text(angle = 90), text = element_text(size = 3)
) + guides(color="none")
p

fname = paste(trait,model,"manhattan.pdf",sep="_")
ggsave(filename = paste(base_folder,"/results/",fname,sep=""), plot = p, device = "pdf", width = 11, height = 8)

## qq-plot
fname = paste(trait,model,"qqplot.png",sep="_")
png(paste(base_folder,"/results/",fname,sep=""), width = 600, height = 600)
qq(res$pvalue)
dev.off()
