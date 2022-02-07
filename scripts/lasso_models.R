## R script to carry out a Lasso-penalised predictive model
## genotype data can be normalised or not (default -and better option- is not!)
# run as Rscript --vanilla lasso_models.R <config_file>

###################################
## read arguments from config file
###################################
# INPUT CONFIGURATION MANAGEMENT ------------------------------------------
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
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
    trait = 'wetter',
    npc = 0, ## n. of PCs to include
    normalise = FALSE, ## n. of PCs to include
    plots = TRUE, ## should plots be plotted out
    test_split = 0.8, ## proportion used as training/testing
    nfolds = 5,
    force_overwrite = FALSE,
    model = "caret"
  ))
}

print("-------------------------------")
print("## CONFIGURATION PARAMETERS ##")
print("-------------------------------")
writeLines(paste(" - trait is:", config$trait))
writeLines(paste(" - model implementation is:", config$model))
writeLines(paste(" - n. of PCs:", config$npc))
writeLines(paste(" - normalisation:", ifelse(config$normalise, "yes", "no")))
print("-------------------------------")

# SETUP -------------------------------------------------------------------
#if not installed, use package remotes, function install_gitlab()
library("caret")
library("glmnet")
library("tidyverse")
library("data.table")

## READING DATA ------------------------------------------------------------
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

print(paste("N. of markers read in: ", nrow(SNP_INFO)))
temp <- genotypes[,-c(1:7,ncol(genotypes)), with = FALSE]
matg <- t(as.matrix(temp))
colnames(matg) <- SNP_INFO$marker
rownames(matg) <- colnames(temp)
rm(temp)

################
## subsampling
################
# writeLines(" - now subsampling the marker loci ...")
# vec <- sample(1:ncol(matg),2000)
# matg <- matg[,vec]
# SNP_INFO <- SNP_INFO[vec,]
# summary(colSums(matg)/nrow(matg))
####

################
## filtering
################
writeLines(" - now filtering the marker data for frequency ...")
vec <- which(colSums(matg)/nrow(matg) < 0.95 & colSums(matg)/nrow(matg) > 0.025)
matg <- matg[,vec]
SNP_INFO <- SNP_INFO[vec,]
print(paste("n. of markers left:", nrow(SNP_INFO)))

nsamples = nrow(matg)
nmarkers = nrow(SNP_INFO)

#################
### phenotypes
writeLines(" - now reading in the phenotype ...")
print(paste("selected trait:", config$trait))
phenotypes <- fread(config$phenotype_file)
phenotypes <- dplyr::select(phenotypes, c(sample, !!as.name(config$trait)))
head(phenotypes)
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% row.names(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## convert score to factor
phenotypes <- mutate(phenotypes, !!as.name(config$trait) := factor(!!as.name(config$trait)))

#################
## kinship matrix
writeLines(' - reading and preparing the kinship matrix')
K <- fread("kinship_gower.csv", sep=",", header = TRUE)
vec <- colnames(K) %in% phenotypes$sample
K <- K[vec,vec, with=FALSE]

#########################
## Principal Components
#########################
if (config$npc > 0) {
  
  writeLines(' - calculating principal components')
  pc <- prcomp(matg)
  # n <- config$npc ## n. of principal components to use for Lasso
  # rownames(pc$x) == rownames(matg)
  ## add PCs to the feature matrix
  matg <- cbind(matg,pc$x[,1:config$npc])
}


## PREPARING THE DATA ----------------------------------------------------------
# Scale data
if (config$normalise == TRUE) {
  
  print("normalising the marker data")
  preprocessParams <- preProcess(matg, method = c("center", "scale"))
  matg <- predict(preprocessParams, matg) 
}

### CARET
### data partition
if (config$model == "caret") {
  
  writeLines(paste(" - running model with", config$model))
  inTrain <- createDataPartition(
    y = dplyr::select(phenotypes, !!as.name(config$trait)) %>% pull(),
    ## the outcome data are needed
    p = config$test_split,
    ## The percentage of data in the
    ## training set
    list = FALSE
  )
  
  X_train <- matg[inTrain, ]
  X_test <- matg[-inTrain, ]
  y_train <- select(phenotypes[inTrain,], all_of(config$trait)) %>% pull()
  y_test <- select(phenotypes[-inTrain], all_of(config$trait)) %>% pull()
  
  ## FIT THE MODEL ----------------------------------------------------------
  # Create and fit Lasso and Ridge objects
  parameters <- seq(0, 1.5, 0.01)
  
  lasso_fit <- train(y= y_train,
                     x = X_train,
                     method = 'glmnet', 
                     tuneGrid = expand.grid(alpha = 1, lambda = parameters),
                     metric = "Accuracy"
  )
  
  print(paste0('Lasso best parameters: ' , lasso_fit$finalModel$lambdaOpt))
  # names(lasso_fit$finalModel)
  
  ## MAKE PREDICTIONS ----------------------------------------------------------
  predictions_lasso <- lasso_fit %>% predict(X_test)
  confm = confusionMatrix(data = predictions_lasso, reference = y_test)
  
  best_lambda = lasso_fit$finalModel$lambdaOpt
  accuracy = confm$overall['Accuracy']
  kappa = confm$overall['Kappa']
  TPR = confm$byClass['Sensitivity']
  TNR = confm$byClass['Specificity']
  print(paste("overall accuracy", accuracy))
  print(confm$table)
  
  res = data.frame("trait"=config$trait,
                   "train_split"=config$test_split,
                   "nsamples"=nsamples,
                   "nmarkers"=nmarkers,
                   "normalise"=ifelse(config$normalise,"yes","no"),
                   "lambda"=best_lambda,
                   "accuracy"=accuracy,
                   "kappa"=kappa,
                   "TPR"=TPR,
                   "TNR"=TNR,
                   "n_pc"=config$npc,
                   "model"="caret")
  
  writeLines(" - saving results to file")
  fname = paste(config$base_folder, "results_lasso.csv", sep="/")
  if(file.exists(fname)) {
    
    fwrite(x = res, file = fname, append = TRUE)
  } else fwrite(x = res, file = fname)
  
  ## MODEL COEFFICIENTS ----------------------------------------------------------
  writeLines(" - extracting model coefficients")
  lasso_coefs = as.data.frame.matrix(coef(lasso_fit$finalModel, lasso_fit$finalModel$lambdaOpt))
  lasso_coefs = filter(lasso_coefs, s1 != 0, !(row.names(lasso_coefs) %in% c("(Intercept)")))
  
  writeLines(" - saving model coefficients to 'dictionary'")
  fname = paste(config$base_folder, "dict_coefs.RData", sep="/")
  
  for (name in rownames(lasso_coefs)) {
    
    if(file.exists(fname)) {
      
      load(fname)
      if (name %in% names(dict_coefs)) {
        
        dict_coefs[name] = dict_coefs[name] + 1
      } else dict_coefs[name] = 1
    } else {
      
      dict_coefs = c(NULL)
      if (name %in% names(dict_coefs)) {
        
        dict_coefs[name] = dict_coefs[name] + 1
      } else dict_coefs[name] = 1
    }
    save(dict_coefs, file = fname)
  }
}


####
## ALERNATIVE APPROACH
###
# lasso_model <- glmnet(x = matg, y = y, alpha = 1, family = "binomial")
# plot(lasso_model)

if (config$model == "glmnet") {
  
  writeLines(" - running alternative implementation of Lasso-penalised logistic regression")
  writeLines(paste(" - running model with", config$model))
  ### data partition
  inTrain <- createDataPartition(
    y = dplyr::select(phenotypes, !!as.name(config$trait)) %>% pull(),
    ## the outcome data are needed
    p = config$test_split,
    ## The percentage of data in the
    ## training set
    list = FALSE
  )
  
  X_train <- matg[inTrain, ]
  X_test <- matg[-inTrain, ]
  y_train <- select(phenotypes[inTrain,], all_of(config$trait)) %>% pull()
  y_test <- select(phenotypes[-inTrain], all_of(config$trait)) %>% pull()
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x = X_train, y = y_train, alpha = 1, nfolds = config$nfolds, type.measure = "class", family = "binomial")
  plot(cv_model)
  
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  print(paste("best lambda value is: ", best_lambda))
  
  #find coefficients of best model
  best_model <- glmnet(x = X_train, y = y_train, alpha = 1, lambda = best_lambda, family = "binomial")
  lasso_coefs = as.data.frame.matrix(coef(best_model))
  lasso_coefs = filter(lasso_coefs, s0 != 0, !(row.names(lasso_coefs) %in% c("(Intercept)")))
  0
  y_predicted <- predict(best_model, s = best_lambda, newx = X_test, type = "class")
  y_pred = as.factor(y_predicted[,1])
  confm <- confusionMatrix(data = y_pred, reference = y_test)
  
  best_lambda = best_lambda
  accuracy = confm$overall['Accuracy']
  kappa = confm$overall['Kappa']
  TPR = confm$byClass['Sensitivity']
  TNR = confm$byClass['Specificity']
  print(confm$table)
  
  res = data.frame("trait"=config$trait,
                   "train_split"=config$test_split,
                   "nsamples"=nsamples,
                   "nmarkers"=nmarkers,
                   "normalise"=ifelse(config$normalise,"yes","no"),
                   "lambda"=best_lambda,
                   "accuracy"=accuracy,
                   "kappa"=kappa,
                   "TPR"=TPR,
                   "TNR"=TNR,
                   "n_pc"=config$npc,
                   "model"="glmnet")
  
  fname = paste(config$base_folder, "results_lasso.csv", sep="/")
  if(file.exists(fname)) {
    
    fwrite(x = res, file = fname, append = TRUE)
  } else fwrite(x = res, file = fname)
  
}

print("DONE!!")
