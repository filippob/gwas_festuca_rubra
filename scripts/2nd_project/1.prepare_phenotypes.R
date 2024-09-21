#!/usr/bin/Rscript
## R script to prepare the file of phenotypes for GWAS  
## phenotypes here are environmental conditions
## in particular, for this project, these are the duration of drought in years

library("tidyverse")
library("data.table")

################################
## parameters
################################
prjfolder = "/home/filippo/Documents/zuzana_festuca_rubra/2nd_project"
fname = "Festuci-sequence-Alpy.xlsx"
sheet_num = 3
exclude_target_locality = NULL
# exclude_target_locality = c("Lavisdalen","Ulvehaugen") ## vector or NULL
keep_shifted_only = TRUE
excl_sample = FALSE
################################
################################

## read and transform
print("reading in the phenotypic data ...")
metadata <- readxl::read_xlsx(file.path(prjfolder, fname), sheet = sheet_num)

if (excl_sample) {
  
  ## samples to be excluded
  excl_samples <- readxl::read_xlsx("samples to be excluded.xlsx")
  metadata <- metadata |> filter(!(`sample no.` %in% excl_samples$`samples to be excluded`))  
}


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
