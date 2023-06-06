#!/bin/bash

# This script grabs the following output files from the final tarball that Stacks pipeline produces
# populations.log.distribs
# .vcf
# populations.samples.fa
# populations.log


##############################################################################################################################
############ DEFINE VARIABLES ################################################################################################



# OUTPUT_POPDISTRIBS is the 
# FILE_POPDISTRIBS is the name that Stacks gives to the log file that we want that populations writes out - this has the total number of loci in it (vcf only has polymorphic loci)
# OUTPUT_VCF is the folder where all of the extracted .vcf files from each assembly parameter combo will be written to and access from later

list_of_datasets=list-extract.txt #name of dataset that we want to process

tempdir=tempplacetoextractto

inpath=/mnt/scratch/rhtoczyd/temp #main path to place where all input directories live
indir=stacks_final_output_from_populations	#directory that holds all of the .tar.gz output from Stacks' module populations
outdir=genetic_data	#folder where all of the extracted populations.XXX files from each assembly parameter combo will be written to
outpath=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen	#main path to place where all output directories should live
file_popdistrib=populations.log.distribs	#name that Stacks gives to file we want that populations writes out
file_poplog=populations.log	#name that Stacks gives to file we want that populations writes out
file_vcf=populations.snps.vcf	#name that Stacks gives to the vcf file that populations writes out
file_fasta=populations.samples.fa	#name that Stacks gives to file we want that populations writes out

# USERS SHOULD NORMALLY NOT NEED TO EDIT BELOW THIS POINT
##############################################################################################################################
##############################################################################################################################


while read dataset
do 
	
	# Get name of dataset folder that tarballs live inside
	run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)
	
	echo Extracting files from dataset $run_name
	
	#navigate to temp dir so files will extract here
	if [ ! -d $inpath/$tempdir ]; then mkdir $inpath/$tempdir; fi
	if [ ! -d $outpath/$run_name/$outdir ]; then mkdir $outpath/$run_name/$outdir; fi	
	cd $inpath/$tempdir
	echo pwd is $(pwd)
	
	# Extract files we want from each tarball, rename each with run info (directory name), move them to an output directory
	for f in $inpath/$run_name/$indir/*.tar.gz; do
	
		#get prefix name for output files (aka name of tarball file)
		prefix=$(echo $f | rev | cut -d "/" -f1 | rev | cut -d "." -f 1)
		echo Processing $prefix
	
		#extract files
		dir=${prefix//$run_name/}
		dir=${dir:1}
		
		tar -xf $f $dir"/"$file_popdistrib
		mv $dir"/"$file_popdistrib $prefix"_"$file_popdistrib
		rm -rf $dir
		
		tar -xf $f $dir"/"$file_poplog
		mv $dir"/"$file_poplog $prefix"_"$file_poplog
		rm -rf $dir
		
		tar -xf $f $dir"/"$file_vcf
		mv $dir"/"$file_vcf $prefix"_"$file_vcf
		rm -rf $dir
		
		tar -xf $f $dir"/"$file_fasta
		mv $dir"/"$file_fasta $prefix"_"$file_fasta
		rm -rf $dir
		
		#move extracted files to outputdir
		mv bioprj_* $outpath/$run_name/$outdir
		wait

	done

	echo ALL DONE EXTRACTING REQUIRED FILES
	#say go copy files

done < ./master_keys/$list_of_datasets
