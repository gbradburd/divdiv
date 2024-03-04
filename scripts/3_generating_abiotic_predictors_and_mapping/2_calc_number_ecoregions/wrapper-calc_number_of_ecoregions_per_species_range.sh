#!/bin/bash

#define variables:
storagenode=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen #path to main node where input files live

gbifdir=$storagenode/speciesLocDat_Sep08_2021_withFilters-cleaned #name of directory where input files with GBIF lats and longs live

generalinputdir=/mnt/home/rhtoczyd/divdiv/data/abiotic/input_and_working #name of directory where misc. needed inputs for all jobs live

jobname=run-ecoregions #label for SLURM book-keeping
executable=calc_number_of_ecoregions_per_species_range.sbatch #script to run

logfilesdir=logfiles_ecoregions #name of directory to create and then write log files to

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=60G #amount of RAM to request/use per CPU
time=168:00:00

list_of_datasets=list-allfinal115.txt #name of dataset that we want to process

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


#for each line of the list_of_datasets file: 
#	execute job that removes low quality reads from each seq. file inside the directories specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1)
	species=$(echo $run_name | cut -d "_" -f3 | sed -r 's/-/ /g')

	#define path to where downloaded raw/unprocessed sequence data live for this dataset
	#and where clean/processed reads should be written to
	indir=$storagenode/$run_name
	outdir=$storagenode/$run_name/location_data

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,GENERALINPUTDIR=$generalinputdir,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,GBIFDIR=$gbifdir,SPECIES="$species" \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$time \
			$executable
			
	echo submitted dataset $run_name from $list_of_datasets
		
done < ../../master_keys/$list_of_datasets

