#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
#list_of_datasets=list-test.txt
storagenode=divdiv_datafiles: 
#note - storagenode needs a / at end if it's a non Google Drive file path
outdirA=/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_gendiv_data
outdirB=/Users/rachel/divdiv/data/popgen/input_and_working/ALL_r80_popgen_data

mkdir -p $outdirA
mkdir -p $outdirB

#grab gstacks logfile for r80 params for each dataset
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	link_name=$(echo "${run_name//bioprj_/}")
	
	r80dir=${storagenode}$run_name/popgen
	r80paramstemp=$(rclone cat $r80dir/r80params.txt)
	r80params=$(echo $r80paramstemp | sed 's/_n[0-9]\{1,\}//g') #remove the _n[#]_ chunk so it matches stacks params formatting in file names elsewhere
		
	echo "run_name: $run_name"
	echo "link_name: $link_name"
	echo "r80params: $r80params"
	
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

		rclone copy divdiv_datafiles:$run_name/gendiv_data --include "WM*${fullr80params}*pars*" $outdirA
		rclone copy divdiv_datafiles:$run_name/gendiv_data --include "WM*${fullr80params}*initPars*" $outdirA
		rclone copy divdiv_datafiles:$run_name/popgen --include "popgenstats.0.5.*${fullr80params}*" $outdirB
	fi
	
	echo " "

done < ../master_keys/$list_of_datasets




