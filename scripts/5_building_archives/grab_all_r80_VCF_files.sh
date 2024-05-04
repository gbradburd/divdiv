#!/bin/bash

#run on umich cluster
#define variables
list_of_datasets=list-allfinal115.txt
storagenode=/nfs/turbo/lsa-bradburd/divdiv_datafiles 
outdirA=$storagenode/dataforpub/all_r80_VCFs
outdirB=$storagenode/dataforpub/all_r80_full_loci_fasta
gdrivenode=divdiv_datafiles:


mkdir -p $outdirA
mkdir -p $outdirB

#grab gstacks logfile for r80 params for each dataset
while read dataset
do
	
	run_name=$(echo $dataset | cut -d " " -f 1)
	link_name=$(echo "${run_name//bioprj_/}")
	
	r80dir=$gdrivenode/$run_name/popgen
	r80paramsfull=$(rclone cat $r80dir/r80params.txt)
	#r80params=$(echo $r80paramsfull | sed 's/_n[0-9]\{1,\}//g') #remove the _n[#]_ chunk so it matches stacks params formatting in file names elsewhere
	
	echo "run_name: $run_name"
	echo "link_name: $link_name"
	echo "r80paramsfull: $r80paramsfull"
	echo "outdirs are: $outdirA $outdirB"
	
	if [ -z "$r80paramsfull" ]
	then
		echo "r80paramsfull variable is empty"
	else 
		
		rsync -av $storagenode/$run_name/genetic_data/*${r80paramsfull}*.snps.vcf $outdirA
		rsync -av $storagenode/$run_name/genetic_data/*${r80paramsfull}*.samples.fa $outdirB
		
	fi
	
	echo " "

done < ./$list_of_datasets
#end run on umich cluster

#run locally
#rsync -avP rhtoczyd@greatlakes.arc-ts.umich.edu:/nfs/turbo/lsa-bradburd/divdiv_datafiles/dataforpub /Users/rachel/divdiv/data/popgen/input_and_working/dataforpub

