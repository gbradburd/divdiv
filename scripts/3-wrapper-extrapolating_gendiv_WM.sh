#!/bin/bash

#define variables:
storagenode=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen #path to main node where input files live

executable=3-run-extrapolating_gendiv_WM.sbatch #script to run

logfilesdir=logfiles_gendivextrap_WM #name of directory to create and then write log files to

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=8G #amount of RAM to request/use per CPU
time=168:00:00

list_of_datasets=list-allfinal136.txt #name of dataset that we want to process

copy_files_to_execute_node=yes #yes to copy input folder to tmp dir on execute node and load files into R from there, no to load files into R directly from where they live on cluster aka $indir below in this file
minpropindivsscoredin=0.5 #percent of indivs that a locus must be present in to save

model_flavor=wishart #value of wishart or cmplnl for which version/flavor of WM model we want to run

jobname=run-gendivextrapWM-${model_flavor} #label for SLURM book-keeping

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


#for each line of the list_of_datasets file: 
#	execute job that removes low quality reads from each seq. file inside the directories specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1)

	#define path to where downloaded raw/unprocessed sequence data live for this dataset
	#and where clean/processed reads should be written to
	indir=$storagenode/$run_name
	outdir=$storagenode/$run_name/gendiv_data
	figdir=$storagenode/gendiv-figures
	
	#if directory to popgen r80 input files is empty; print warning, otherwise process the files that are there
	n_files=($(ls $indir/r80_outputs/popgenstats* | wc -l))
	if [ $n_files = 0 ]
	then echo WARNING - there are no popgen input files for $run_name, go investigate
	else

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,FIGDIR=$figdir,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,MINPROPINDIVSSCOREDIN=$minpropindivsscoredin,MODEL_FLAVOR=$model_flavor \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$time \
			$executable
			
	echo submitted dataset $run_name from $list_of_datasets
	fi
	
done < ./master_keys/$list_of_datasets
