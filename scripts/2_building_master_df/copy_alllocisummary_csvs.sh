#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
outdir=/Users/rachel/divdiv/data/methodological/input_and_working/all_locisummary_csvs

mkdir -p $outdir

#grab allparams-locisummary_stats-bioprj_XXX.csv for each dataset 
#(.csv contains values for all Stacks param combos, will filter to just r80 in R script that builds methodological preds df)
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	
	echo "run_name: $run_name"
	echo "outdir: $outdir"
	
	rclone copy divdiv_datafiles:$run_name/popgen --include "allparams-locisummary_stats*" $outdir
	
	echo " "

done < ../master_keys/$list_of_datasets