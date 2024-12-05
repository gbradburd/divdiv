#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
outdir=/Users/rachel/divdiv/data/methodological/input_and_working/all_r80_bpstats
storagenode=divdiv_datafiles: 
#note - storagenode needs a / at end if it's a non Google Drive file path


mkdir -p $outdir

#grab gstacks logfile for r80 params for each dataset
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	link_name=$(echo "${run_name//bioprj_/}")
	
	r80dir=$storagenode/$run_name/popgen
	r80paramstemp=$(rclone cat $r80dir/r80params.txt)
	r80params=$(echo $r80paramstemp | sed 's/_n[0-9]\{1,\}//g') #remove the _n[#]_ chunk so it matches stacks params formatting in file names elsewhere
	
	echo "run_name: $run_name"
	echo "link_name: $link_name"
	echo "r80params: $r80params"
	echo "outdir: $outdir"
	
	if [ -z "$r80params" ]
	then
		echo "r80params variable is empty"
	else 
		
		littlem=$(echo $r80params | cut -d _ -f 2)
		bigm=$(echo $r80params | cut -d _ -f 4)
		nname=$(echo $r80params | cut -d _ -f 5)
	
		if [ $nname = "nis1lessthanM" ]; then n=$(echo $(($bigm-1))); fi
		if [ $nname = "nequalM" ]; then n=$(echo $bigm); fi
		if [ $nname = "nis1morethanM" ]; then n=$(echo $(($bigm+1))); fi
	
		fullr80params="littlem_"${littlem}"_bigm_"${bigm}"_n"${n}"_"${nname}
		echo "fullr80params: $fullr80params"
		
		rclone copy ${storagenode}${run_name}/popgen $outdir --include "bpstats.0.5.*${fullr80params}*.Robj" --progress
	fi
	
	echo " "

done < ../master_keys/$list_of_datasets