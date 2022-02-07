#!/bin/sh

project_folder='Documents/zuzana_festuca_rubra'
path_to_script=${HOME}/${project_folder}/gwas_festuca_rubra/scripts
rscript=lasso_models.R
config=config.txt

for i in `seq 1 2`;
do
	echo "########################"
	echo "running replicate n. $i";
	echo "########################"
	Rscript --vanilla ${path_to_script}/${rscript} $config;
done

