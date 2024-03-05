#!/bin/bash

#define variables
list_of_datasets=list-allfinal115.txt
storagenode=divdiv_datafiles: 
#note - storagenode needs a / at end if it's a non Google Drive file path
#note - divdiv_datafiles: is name of rclone remote that references/points to this Shared Drive on UMich Google Drive
outdir1=/Users/rachel/divdiv/data/abiotic/input_and_working/marmapdists-output
outdir2=/Users/rachel/divdiv/data/abiotic/input_and_working/ecoregions-output
outdir3=/Users/rachel/divdiv/data/abiotic/input_and_working/postmarmap_coords

mkdir -p $outdir1
mkdir -p $outdir2
mkdir -p $outdir3

#grab marmap output and ecoregion output for each dataset in list specified above (these files are nested within bioprj folders)
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	link_name=$(echo "${run_name//bioprj_/}")
	
	echo "run_name: $run_name"
	echo "link_name: $link_name"
	
	#marmap
	rclone copy divdiv_datafiles:$run_name/location_data --include "*max_and_pw*" $outdir1 -P
	#ecoregion
	rclone copy divdiv_datafiles:$run_name/location_data --include "*number_of_ecoregions*" $outdir2 -P
	#final points after moving points on land to water
	rclone copy divdiv_datafiles:$run_name/location_data --include "*coordspostmarmap*" $outdir3 -P
	
	echo " "

done < ../../master_keys/$list_of_datasets

