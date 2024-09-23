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

metadata <- metadata %>%
  select(-Sample) %>%
  rename(sample = `Sample code`)

## write out file
print("writing out the file ...")
fwrite(x = metadata, file = file.path(prjfolder, "phenotypes.csv"), sep = ",", col.names = TRUE)

print("DONE!")
