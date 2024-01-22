#!/bin/bash

#define variables:
storagenode=/mnt/scratch/$USER #path to main node where fastq files will get written to / stored

jobname=run-downloadfastqs #label for SLURM book-keeping
executable=1-run-download_fastqs_from_ncbi_sra.sbatch #script to run

logfilesdir=logfiles_sampledownload #name of directory to create and then write log files to

cpus=12 #number of CPUs to request/use per batch download
ram_per_cpu=2G #amount of RAM to request/use per CPU
runtime=168:00:00 #max amount of time to let job run for

masterkeydir=/home/rhtoczyd/divdiv/scripts/master_keys #file path to master_keys folder (that holds text list of datasets to process)
listdir=/home/rhtoczyd/divdiv/data/bioinformatics/lists_to_download #file path to folder that contains text lists of reads to download for each dataset

list_of_datasets=list_of_datasets-TEST.txt #name of file that contains the name of all of the datasets that we want to download, aka each of files in this file contains sequence run accession #'s to download fo a bioprj/species combo
expectedreadcountlist=allexpectedreadcounts-6-3-2021.txt #name of files that contains the expected number of reads that each downloaded file should have according to NCBI metadata

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


#for each line of the list_of_datasets file: 
#	create a directory on storage node, 
#	execute job that downloads the SRAs numbers listed in one of the files (id_list) specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)
	
	#get whether data are paired end or single end reads
	read_type=$(echo $dataset | cut -d " " -f2)
	
	#redefine dataset to be just bioprj_<prj>_<species>.txt so executable will work correctly
	dataset=$(echo $dataset | cut -d " " -f1)

	#define path to write downloaded sequence data to
	rawreadsdir=$storagenode/$run_name/1_unprocessed_reads

	#check if output directory to download files to exists; if not, make one
	if [ ! -d $storagenode/$run_name ]; then mkdir $storagenode/$run_name; fi
	if [ ! -d $rawreadsdir ]; then mkdir $rawreadsdir; fi
	

	#submit job to cluster, one job per dataset
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,READ_TYPE=$read_type,STORAGENODE=$storagenode,SRATOOLS=$sratools,RAWREADSDIR=$rawreadsdir,DATASET=$dataset,EXPECTEDREADCOUNTLIST=$expectedreadcountlist,MASTERKEYDIR=$masterkeydir,LISTDIR=$listdir \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$runtime \
			$executable
			
	echo submitted ID list with name $dataset from $list_of_datasets and it has ${read_type}-end data
			
done < $masterkeydir/$list_of_datasets
