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
    trait = 'warmer',
    normalise = FALSE, ## n. of PCs to include
    plots = TRUE, ## should plots be plotted out
    test_split = 0.8, ## proportion used as training/testing
    force_overwrite = FALSE
  ))
}

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

#################
### phenotypes
print("now reading in the binary phenotype ...")
phenotypes <- fread(config$phenotype_file)
phenotypes <- dplyr::select(phenotypes, c(sample, !!as.name(config$trait)))
head(phenotypes)
print(paste(nrow(phenotypes),"records read from the phenotype file",sep=" "))

phenotypes <- phenotypes[phenotypes$sample %in% row.names(matg),]
print(paste(nrow(phenotypes),"records read from the phenotype file after alignment with genotypes",sep=" "))

## convert score to factor
phenotypes <- mutate(phenotypes, !!as.name(config$trait) := factor(!!as.name(config$trait)))


## PREPARING THE DATA ----------------------------------------------------------
# Scale data
if (normalise == TRUE) {
 
  preprocessParams <- preProcess(matg, method = c("center", "scale"))
  matg <- predict(preprocessParams, matg) 
}

### data partition
inTrain <- createDataPartition(
  y = phenotypes$warmer,
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
parameters <- seq(0, 1, 0.01)

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
names(confm)

best_lambda = lasso_fit$finalModel$lambdaOpt
accuracy = confm$overall['Accuracy']
kappa = confm$overall['Kappa']
TPR = confm$byClass['Sensitivity']
TNR = confm$byClass['Specificity']
confm$table

## MODEL COEFFICIENTS ----------------------------------------------------------
lasso_coefs = as.data.frame.matrix(coef(lasso_fit$finalModel, lasso_fit$finalModel$lambdaOpt))
lasso_coefs = filter(lasso_coefs, s1 != 0)
nrow(lasso_coefs)

dict_coefs = c("q" = 1, "w" = 0, "d" = 2)
dict_coefs["q"] = dict_coefs["q"] + 1
dict_coefs["h"] = 1



for (name in rownames(lasso_coefs)) {
  
    
    dict_coefs[]
}

####
## ALERNATIVE APPROACH
###
lasso_model <- glmnet(x = matg, y = y, alpha = 1, family = "binomial")
plot(lasso_model)

### data partition
inTrain <- createDataPartition(
  y = phenotypes$warmer,
  ## the outcome data are needed
  p = .80,
  ## The percentage of data in the
  ## training set
  list = FALSE
)

X_train <- matg[inTrain, ]
X_test <- matg[-inTrain, ]
y_train <- select(phenotypes[inTrain,], all_of(config$trait)) %>% pull()
y_test <- select(phenotypes[-inTrain], all_of(config$trait)) %>% pull()

#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x = X_train, y = y_train, alpha = 1, nfolds = 5, type.measure = "class", family = "binomial")
plot(cv_model)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda

#find coefficients of best model
best_model <- glmnet(x = X_train, y = y_train, alpha = 1, lambda = best_lambda, family = "binomial")
lasso_coefs = as.data.frame.matrix(coef(best_model))
lasso_coefs = filter(lasso_coefs, s0 != 0)

y_predicted <- predict(best_model, s = best_lambda, newx = X_test, type = "class")
y_pred = as.factor(y_predicted[,1])
confusionMatrix(data = y_pred, reference = y_test)

