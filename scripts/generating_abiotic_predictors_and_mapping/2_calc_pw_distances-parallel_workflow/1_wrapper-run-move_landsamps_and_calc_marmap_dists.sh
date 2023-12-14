#!/bin/bash

#define variables:
storagenode=/mnt/research/Fitz_Lab/bradburd/rht/divdiv_working_popgen #path to main node where input files live

gbifdir=$storagenode/speciesLocDat_Sep08_2021_withFilters-cleaned #name of directory where input files with GBIF lats and longs live

bathydir=/mnt/home/rhtoczyd/divdiv/data/abiotic/input_and_working/bathy_Robjs #name of directory where world bathy Robjs live - resistance/depth layers

logfilesdir=logfiles_marmapdists #name of directory to create and then write log files to

cpus=2 #number of CPUs to request/use per dataset
ram_per_cpu=28G #amount of RAM to request/use per CPU
time=168:00:00

n_landsamps_per_job=10 #number of landsamp points to move to water per hpcc array job

list_of_datasets=list-allfinal115.txt #name of dataset that we want to process


copy_files_to_execute_node=yes #yes to copy input folder to tmp dir on execute node and load files into R from there, no to load files into R directly from where they live on cluster aka $indir below in this file

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
	echo "starting to submit jobs for dataset: $run_name"

	#define path to where downloaded raw/unprocessed sequence data live for this dataset
	#and where clean/processed reads should be written to
	indir=$storagenode/$run_name
	outdir=$storagenode/$run_name/location_data
	
	#count number of unique biosamples to look up
	n_landsamps=($(wc -l $outdir/landsamps-gbif_and_genetic-${run_name}.csv))
	n_landsamps=$((n_landsamps - 1)) #don't count header
	echo "I have $n_landsamps samples on land that I need to move to water"

	#calculate number of jobs needed for moving landsamps parallelization
	n_movelandsamp_jobs=$(((n_landsamps + n_landsamps_per_job - 1) / n_landsamps_per_job))
	echo -e "will submit $n_movelandsamp_jobs jobs to move landsamps to water with $n_landsamps landsamp points per job for this dataset"


	#--------------------------------------
	jobname=run-marmap-movelandsamps #label for SLURM book-keeping
	executable=move_landsamps_to_water.sbatch #script to run

	#submit job to cluster
	jid_move=$(sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,BATHYDIR=$bathydir,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,GBIFDIR=$gbifdir,SPECIES="$species",N_LANDSAMPS=$n_landsamps,N_LANDSAMPS_PER_JOB=$n_landsamps_per_job \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--array=[0-$((n_movelandsamp_jobs-1))]%$n_movelandsamp_jobs \
			--output=./$logfilesdir/${jobname}_${run_name}_%A_batch_%a.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A_batch_%a.err \
			--time=$time \
			--constraint="NOAUTO:lac|amr|skl" \
			$executable \
			|awk -v "nID=$n_movelandsamp_jobs" '{for (i=0; i<nID; i++) printf(":%s",$4"_"i)}')
			#| cut -d " " -f 4)
	#jid_move=:${jid_move}
	echo "jid_move is $jid_move"		
	echo " "	

	# keep tacking all step6 jobIDs onto one variable to pass to error checker as dependency
	#declare "jid_allforerrorcheck=${jid_step7}${jid_allforerrorcheck}"

	
	#--------------------------------------
	jobname=run-marmap-calcdists #label for SLURM book-keeping
	executable=calc_pw_gcd_and_sea_dists_and_plot.sbatch #script to run

	#submit job to cluster
	sbatch --job-name=$jobname \
			--export=CPUS=$cpus,RUN_NAME=$run_name,STORAGENODE=$storagenode,BATHYDIR=$bathydir,INDIR=$indir,OUTDIR=$outdir,LOGFILESDIR=$logfilesdir,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,GBIFDIR=$gbifdir,SPECIES="$species" \
			--cpus-per-task=$cpus \
			--mem-per-cpu=$ram_per_cpu \
			--output=./$logfilesdir/${jobname}_${run_name}_%A.out \
			--error=./$logfilesdir/${jobname}_${run_name}_%A.err \
			--dependency=afterok$(eval echo \$jid_move) \
			--kill-on-invalid-dep=yes \
			--time=$time \
			--constraint="NOAUTO:lac|amr|skl" \
			$executable
			
	echo submitted dataset $run_name from $list_of_datasets
	echo " "
	echo " "
		
done < ../../master_keys/$list_of_datasets
