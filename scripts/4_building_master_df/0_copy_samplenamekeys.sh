#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
outdir=/Users/rachel/divdiv/data/all_samplenamekeys
storagenode=divdiv_datafiles: 
#note - storagenode needs a / at end if it's a non Google Drive file path


mkdir -p $outdir

#grab samplenamekey txt file for each dataset (and append run_name, since files not uniquely named on Google Drive)
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
		
	echo "run_name: $run_name"
	echo "outdir: $outdir"
	
	rclone copyto divdiv_datafiles:$run_name/samplenamekey.txt $outdir/samplenamekey-${run_name}.txt
	
	echo " "

done < ../master_keys/$list_of_datasets




