#!/bin/bash --login

#define variables:
storagenode=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen #path to main node where outputs will end up
#storagenode=/mnt/scratch/$USER

keysdir=/mnt/home/rhtoczyd/divdiv/data/all_samplenamekeys #place where samplenamekey live
alllatlong=/mnt/home/rhtoczyd/divdiv/data/final_genetic_latlongs.csv #file with all the final lat/longs in it for all samples
outdirall=/mnt/scratch/rhtoczyd/ALL_fst #place to write all outputs to across all datasets (also writing to each indiv r80 output folder per dataset)

jobname=run-calcfst #label for SLURM book-keeping
executable=1-run-calc_fst_and_related_metrics.sbatch #script to run

logfilesdir=logfiles_calcfst #name of directory to create and then write log files to

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=66G #amount of RAM to request/use per CPU
time=168:00:00

list_of_datasets=list-allfinal93.txt #name of dataset that we want to process

minpropindivsscoredin=0.5 #percent of indivs that a locus must be present in to save

copy_files_to_execute_node=yes #yes to copy input folder to tmp dir on execute node and load files into R from there, no to load files into R directly from where they live on cluster aka $indir below in this file

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi
if [ ! -d ./$outdirall ]; then mkdir ./$outdirall; fi


#for each line of the list_of_datasets file: 
#	execute job that removes low quality reads from each seq. file inside the directories specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)

	#define file paths
	r80dir=/mnt/scratch/rhtoczyd/$run_name/r80_outputs
	outdir=/mnt/scratch/rhtoczyd/$run_name/r80_outputs

	##if directory to downloaded files doesn't contain at least 1 .gz file; print warning, otherwise process the files that are there
	n_popgenfiles=($(ls $indir/popgen* | wc -l))
	if [ $n_popgenfiles = 0 ]
	then echo WARNING - there are no VCF files in $indir, go investigate
	else

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,KEYSDIR=$keysdir,ALLLATLONG=$alllatlong,OUTDIRALL=$outdirall,R80DIR=$r80dir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,MINPROPINDIVSSCOREDIN=$minpropindivsscoredin,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$time \
			$executable
			
	echo submitted ID list with name $run_name from $list_of_datasets
	fi	
		
done < ../master_keys/$list_of_datasets
