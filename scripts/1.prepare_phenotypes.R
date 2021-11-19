#!/usr/bin/Rscript
## R script to prepare the file of phenotypes for GWAS  

library("tidyverse")
library("data.table")

## read and transform
print("reading in the phenotypic data ...")
metadata <- readxl::read_xlsx("288_samples_STR_ALL-final2.xlsx", sheet = 2)

## create the columns of traits: wetter, warmer (binary), warmer and wetter (multinomial)
metadata <- metadata %>% 
  rename(local_foreign = `local0/foreign1`,
         original_temperature = `original Temperature`,
         original_moisture = `original moisture`,
         locality_target = `locality target`) %>% 
  filter(`shifted populations` == 4) %>%
  mutate(local_foreign = as.factor(local_foreign),
         wetter = ifelse(wetter == 2,1,wetter),
         warm_wet = interaction(warmer,wetter),
         temp_moist = interaction(original_temperature, original_moisture))

metadata <- metadata %>%
  select(`sample no.`, warmer, wetter, warm_wet,origin) %>%
  rename(sample = `sample no.`)

## write out file
print("writing out the file ...")
fwrite(x = metadata, file = "phenotypes.csv", sep = ",", col.names = TRUE)

print("DONE!")