#!/usr/bin/Rscript
## R script to prepare the file of phenotypes for GWAS  

library("tidyverse")
library("data.table")

## parameters
exclude_target_locality = c("Lavisdalen","Ulvehaugen") ## vector or NULL
keep_shifted_only = FALSE

## read and transform
print("reading in the phenotypic data ...")
metadata <- readxl::read_xlsx("288_samples_STR_ALL-final2.xlsx", sheet = 2)

## create the columns of traits: wetter, warmer (binary), warmer and wetter (multinomial)
metadata <- metadata %>% 
  rename(local_foreign = `local0/foreign1`,
         original_temperature = `original Temperature`,
         original_moisture = `original moisture`,
         locality_target = `locality target`) %>% 
  mutate(local_foreign = as.factor(local_foreign),
         wetter = ifelse(wetter == 2,1,wetter),
         warm_wet = interaction(warmer,wetter),
         temp_moist = interaction(original_temperature, original_moisture))

if (keep_shifted_only) metadata <- filter(metadata, `shifted populations` == 4)
  
if (length(exclude_target_locality) > 0) {
  
  writeLines(" - removing unnecessary samples ")
  print(paste("samples to remove: ", exclude_target_locality, collapse = ","))
  metadata = filter(metadata, !(locality_target %in% exclude_target_locality))
  print(paste("n. of records left:", nrow(metadata)))
}

metadata <- metadata %>%
  select(`sample no.`, warmer, wetter, warm_wet,origin, targetTemp, targetMois) %>%
  rename(sample = `sample no.`)

## write out file
print("writing out the file ...")
fwrite(x = metadata, file = "phenotypes.csv", sep = ",", col.names = TRUE)

print("DONE!")