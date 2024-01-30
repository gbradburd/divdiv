#!/bin/bash

#define variables:
storagenode=/mnt/scratch/$USER

jobname=run-searchfor_common_adapters_firsttime #label for SLURM book-keeping
executable=3-run-searchfor_common_adapters.sbatch #script to run

logfilesdir=logfiles_searchforcommonadapter #name of directory to create and then write log files to

resultsfilesdir=/home/rhtoczyd/divdiv/data/bioinformatics/summary_files_of_results #directory to store summary results file for each dataset

cpus=24 #number of CPUs to request/use per dataset
ram_per_cpu=250M #amount of RAM to request/use per CPU
time=168:00:00 #max amount of time (wall clock) to let job run for

list_of_datasets=list_of_datasets-TEST.txt #name of dataset that we want to process

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi

#check if output directory for summary results .csvs has been created in submit dir yet; if not, make one
if [ ! -d ./$resultsfilesdir ]; then mkdir ./$resultsfilesdir; fi


#for each line of the list_of_datasets file: 
#	execute job that "removes" adapters from each seq. file inside the directories specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)
	
	#get whether data are paired end or single end reads
	read_type=$(echo $dataset | cut -d " " -f2)
	
	#redefine dataset to be just bioprj_<prj>_<species>.txt so executable will work correctly
	dataset=$(echo $dataset | cut -d " " -f1)

	#define path to where sequence data with low quality reads removed live for this dataset
	#and where clean/processed reads should be written to
	cleanreadsdir=$storagenode/$run_name/2_lowqualrmvd_reads
	finalreadsdir=$storagenode/$run_name/3_adapterrmvd_reads
	
	##if directory to downloaded files doesn't contain at least 1 .gz file; print warning, otherwise process the files that are there
	n_gzfiles=($(ls $cleanreadsdir/*.gz | wc -l))
	if [ $n_gzfiles = 0 ]
	then echo WARNING - there are no .gz files in $cleanreadsdir, go investigate
	else


	#submit adapter removal job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,CLEANREADSDIR=$cleanreadsdir,FINALREADSDIR=$finalreadsdir,LOGFILESDIR=$logfilesdir,JOBNAME=$jobname,RESULTSFILESDIR=$resultsfilesdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$time \
			$executable
	
	
	echo submitted ID list with name $dataset from $list_of_datasets
	fi		
	
done < ../../master_keys/$list_of_datasets

