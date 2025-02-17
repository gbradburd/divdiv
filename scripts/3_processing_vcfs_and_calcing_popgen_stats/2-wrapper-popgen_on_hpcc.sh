#!/bin/bash

#define variables:
storagenode=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen #path to main node where outputs will end up
#storagenode=/mnt/scratch/$USER

jobname=run-popgenstats #label for SLURM book-keeping
executable=2-run-popgen_on_hpcc.sbatch #script to run

logfilesdir=logfiles_popgenstats #name of directory to create and then write log files to

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=66G #amount of RAM to request/use per CPU
#ram_per_cpu=250G #amount of RAM to request/use per CPU
#ram_per_cpu=20G #amount of RAM to request/use per CPU
#ram_per_cpu=491G #amount of RAM to request/use per CPU
time=168:00:00

#list_of_datasets=list-allfinal115.txt #name of dataset that we want to process
list_of_datasets=list-final88.txt #name of dataset that we want to process

nPCs=4 #number of principal component axes to save
minpropindivsscoredin=0.5 #percent of indivs that a locus must be present in to save

copy_files_to_execute_node=yes #yes to copy input folder to tmp dir on execute node and load files into R from there, no to load files into R directly from where they live on cluster aka $indir below in this file

manualsampstodrop=/mnt/home/rhtoczyd/divdiv/scripts/master_keys/additional_samples_tossed-alsolistedin_master_tossed.csv #list of additional samples to drop, with sampleX names (rather than SRR IDs)

#---------------------------------------------------------

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


#for each line of the list_of_datasets file: 
#	execute job that removes low quality reads from each seq. file inside the directories specified in the list_of_datasets file  
while read dataset
do 
	
	#define label to give dataset
	run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)

	#define file paths
	indir=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen/$run_name/genetic_data
	outdir=/mnt/scratch/rhtoczyd/$run_name/popgen
	figdir=/mnt/scratch/rhtoczyd/popgen-figures
	keysdir=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen/$run_name #place where samplenamekey and lat long tables of genetic pts live

	##if directory to downloaded files doesn't contain at least 1 .gz file; print warning, otherwise process the files that are there
	n_gzfiles=($(ls $indir/*.snps.vcf | wc -l))
	if [ $n_gzfiles = 0 ]
	then echo WARNING - there are no VCF files in $indir, go investigate
	else

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,INDIR=$indir,OUTDIR=$outdir,FIGDIR=$figdir,LOGFILESDIR=$logfilesdir,NPCS=$nPCs,MINPROPINDIVSSCOREDIN=$minpropindivsscoredin,KEYSDIR=$keysdir,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,MANUALSAMPSTODROP=$manualsampstodrop \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--time=$time \
			--constraint="NOAUTO:lac|amr|skl" \
			$executable
			
	echo submitted ID list with name $run_name from $list_of_datasets
	fi	
		
done < ../master_keys/$list_of_datasets
