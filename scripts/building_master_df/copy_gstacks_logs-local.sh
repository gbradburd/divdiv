#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
outdir=/Users/rachel/divdiv/data/all_r80_gstacksoutlogs

mkdir -p $outdir

#grab gstacks logfile for r80 params for each dataset
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	link_name=$(echo "${run_name//bioprj_/}")
	
	r80dir=$storagenode/$run_name/popgen
	r80paramstemp=$(rclone cat divdiv_datafiles:$run_name/popgen/r80params.txt)
	r80params=$(echo "${r80paramstemp//_n[0-9]/}")
	
	echo "run_name: $run_name"
	echo "link_name: $link_name"
	echo "r80params: $r80params"
	echo "outdir: $outdir"
	
	if [ -z "$r80params" ]
	then
		echo "r80params variable is empty"
	else 
		rclone copy divdiv_datafiles:$run_name/${link_name}-LOGS/logfiles_gstacks --include "*${r80params}*.out" $outdir
	fi
	
	echo " "

done < ../master_keys/$list_of_datasets




