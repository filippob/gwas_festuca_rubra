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
sheet_num = 4
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
  select(-set) %>%
  rename(sample = code) %>%
  filter(measurement == 2) %>%
  select(-measurement)

X <- metadata[,-1]
varnames = names(X)
std = X |> summarise(across(everything(), sd))
avg = colMeans(X)
X <- as.matrix(X)

X = sweep(x = X, MARGIN = 2, avg)
X = X %*% diag(1 / std)

X <- as_tibble(X)
names(X) = varnames
X |> summarise(across(where(is.numeric), list(avg = mean,std = sd)))

metadata = cbind.data.frame(metadata[,1], X)

## write out file
print("writing out the file ...")
fwrite(x = metadata, file = file.path(prjfolder, "phenotypes_morphology.csv"), sep = ",", col.names = TRUE)

print("DONE!")
