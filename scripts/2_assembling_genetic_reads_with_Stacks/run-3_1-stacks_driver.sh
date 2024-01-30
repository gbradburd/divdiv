#!/bin/bash

##############################################################################################################################
############ DEFINE VARIABLES ################################################################################################

# Note - edit variables in this section but should not need to 
#		 edit any other sections of this script (unless maybe when switching to a new cluster)

# Note - if moving to a new cluster, check that file paths for the below variables are correct in .sbatch scripts:
#			WORKDIR; INDIR; OUTDIR are all currently /tmp/local/... 	#file path on execute node where work is done

run_name=testdataset	#name of project, all files will live in top directory with this name

vector_of_littlem_values=( 3 )	#the little m (m) ustacks assembly parameter that you want to test here (likely the optimal value identified in Section 1)
vector_of_bigm_values=( 3 6 9 )	#list all values of the big m (M) ustacks assembly parameter that you want to test here
max_loc=3	#ustacks parameter
vector_of_subname_values=( nequalM nis1lessthanM nis1morethanM )	#list which of three relational values of little n that you want to test
																	#pick one or more of this list: nequalM; nis1lessthanM; nis1morethanM

storagenode=/mnt/scratch/$USER	#file path to storage node on cluster - aka where input and output tarballs are stored in between jobs on execute nodes

samp_storage_dir=$storagenode/$run_name/samples	#file path to directory that contains gzipped fastq files for ustacks

n_samps_per_ustacks=30	#number of samples to run through ustacks per job
n_samps_per_cstacks=30	#number of samples to add to the cstacks catalog per iteration of running cstacks
n_samps_per_sstacks=30	#number of samples to run through sstacks per job


copy_files_to_execute_node=yes # can take value yes or no, determines whether input files are copied directly to the execute node, and then back to storage node (yes) 
							   # or always left on storage node while Stacks is running (no)

paired_end=no	#define if samples are paired end (yes) or single end (no)

modules_to_load=(GCC/9.3.0 OpenMPI/4.0.3 Stacks/2.54) #list of all modules to load, listed in order to load them, separated by spaces


##############################################################################################################################
############ SET COMPUTING RESOURCE REQUESTS #################################################################################

#ustacks
cpus_ustacks=6
ram_ustacks=16G
time_ustacks=168:00:00

#gather sample tarballs
cpus_gather=2
ram_gather=8G
time_gather=168:00:00

#cstacks0
cpus_cstacks0=2
ram_cstacks0=82G
time_cstacks0=168:00:00

#cstacks1
cpus_cstacks1=2
ram_cstacks1=250G
time_cstacks1=168:00:00

#split sample tarballs
cpus_split=2
ram_split=8G
time_split=168:00:00

#sstacks + tsv2bam
cpus_sstacks=10
ram_sstacks=84G
time_sstacks=168:00:00

#gstacks
cpus_gstacks=12
ram_gstacks=132G
time_gstacks=168:00:00

#populations
cpus_populations=6
ram_populations=24G
time_populations=168:00:00

#error check
cpus_errorcheck=1
ram_errorcheck=8G
time_errorcheck=3:00:00


# USERS SHOULD NOT NEED TO EDIT BELOW THIS POINT UNDER MOST CIRCUMSTANCES
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


#count number of input fastq files
n_samps=($(printf "%s\n" $samp_storage_dir/*.gz | grep -v "\.2\." | wc -l))
echo I have $n_samps input .fq.gz sample files \(for ustacks\)

#calculate number of jobs needed for ustacks and sstacks, round up to whole number
n_ustacks_jobs=$(((n_samps + n_samps_per_ustacks - 1) / n_samps_per_ustacks))
n_sstacks_jobs=$(((n_samps + n_samps_per_sstacks - 1) / n_samps_per_sstacks))

echo -e "will submit $n_ustacks_jobs ustacks jobs with $n_samps_per_ustacks samples per job for each value of big_m and little_m"
echo -e "will submit $n_sstacks_jobs sstacksandtsv2bam jobs with $n_samps_per_sstacks samples per job for each value of big_m and little_m"


### split popmap up into chunks for cstacks step ### -----------------------------------------------------------------------------------------------------------------------

#This section finds the popmap file in $storagenode/$run_name/popmap (file should be called "popmap") ...
#and then splits it up into sub popmaps named popmap0, popmap1, etc. ...
#using the variable n_samps_per_cstacks to determine how many samples should be in each sub popmap.
#the number of files in each sub popmap are the number of files that will be added to the cstacks catalog each time (e.g. during cstacks0.sbatch job)

num_cstacks1jobs=$(perl -s   -e '
      
#### perl program starts here
use Carp;
use File::Basename;

my $popMapFile = "$storagenode/$run_name/popmap";
my $outdir = "$storagenode/$run_name";

## call subroutine to split the popmap file
my $numSplitPopmapFiles = splitPopmap($popMapFile,$n_samps_per_cstacks,$outdir);
my $numCstacks1Jobs = $numSplitPopmapFiles-1;
print  $numCstacks1Jobs;

### splitPopmap subroutine
sub splitPopmap {
    my $popmapFile = shift;
    my $samplesPerJob = shift;
    my $outdir = shift;

    open POPMAP, $popmapFile
	or croak "Cannot open popmap file $popmapFile";
    my @samples;
    while(<POPMAP>) {
        ## parse popmap file, allowing comments and blank lines;
	chomp;
	next if /^#/;
	next if /^\s*$/;
	push @samples, $_;
    }
    close POPMAP;
    my $basename = basename($popmapFile);
    my $index = 0;
    my $numFiles;
    while(1) {
        ## write each chunk to a separate file;
	my $popmapOutfile = "$outdir/$basename.$index";
	open POPMAP, ">$popmapOutfile"
	    or croak "Cannot write popmap file to $popmapOutfile";
	$numFiles ++;
	my @out = splice(@samples,0,$samplesPerJob);
	print POPMAP join("\n", @out), "\n";
	close POPMAP;
	last if (scalar @samples == 0);
	$index ++;
    } ## while
    return $numFiles;
};'  -- -storagenode=$storagenode -run_name=$run_name -n_samps_per_cstacks=$n_samps_per_cstacks)

### end of perl program

#slurm variable key:
# %A = SLURM_ARRAY_JOB_ID
# %a = SLURM_ARRAY_TASK_ID

### run ustacks ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-ustacks
logfilesdir=logfiles_ustacks
executable=run-ustacks.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
	jid_ustacks=$(sbatch --job-name=$jobname \
			--export=N_SAMPS_PER_USTACKS=$n_samps_per_ustacks,CPUS=$cpus_ustacks,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_ustacks \
			--mem=$ram_ustacks \
			--array=[0-$((n_ustacks_jobs-1))]%$n_ustacks_jobs \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_batch%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_batch%a.err \
			--time=$time_ustacks \
			$executable \
			|awk -v "nID=$n_ustacks_jobs" '{for (i=0; i<nID; i++) printf(":%s",$4"_"i)}')
	declare "jid_ustacks_${little_m}_${big_m}=$jid_ustacks" #this stores the job id assigned to the job with little m = 1 and big m = 2 to jid_ustacks_1_2
	eval echo submitted ustacks with jobID \$jid_ustacks_${little_m}_${big_m}
	done
done


### run gathertarballs ### -----------------------------------------------------------------------------------------------------------------------

#define variables: 
jobname=run-gather_indiv_samps_pre-cstacks
logfilesdir=logfiles_gathertarballs
executable=run-gather_indiv_samps_pre-cstacks.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
	jid_gathertarball=$(sbatch --job-name=$jobname \
	--export=CPUS=$cpus_gather,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
	--cpus-per-task=$cpus_gather \
	--mem=$ram_gather \
	--array=[0]%1 \
	--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_%a.out \
	--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_%a.err \
	--dependency=afterok$(eval echo \$jid_ustacks_${little_m}_${big_m}) \
	--kill-on-invalid-dep=yes \
	--time=$time_gather \
	$executable \
	|cut -d " " -f 4)
	declare "jid_gathertarball_${little_m}_${big_m}=${jid_gathertarball}_0"
	eval echo submitted gathertarballs with jobID \$jid_gathertarball_${little_m}_${big_m}
	done
done


### run cstacks0 ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-cstacks0
logfilesdir=logfiles_cstacks
executable=run-cstacks0.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
		jid_cstacks0=$(sbatch --job-name=$jobname \
			--export=CPUS=$cpus_cstacks0,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_cstacks0 \
			--mem=$ram_cstacks0 \
			--array=[0]%1 \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.err \
			--dependency=afterok:$(eval echo \$jid_gathertarball_${little_m}_${big_m}) \
			--kill-on-invalid-dep=yes \
			--time=$time_cstacks0 \
			$executable \
			|cut -d " " -f 4)
		declare "jid_cstacks0_${little_m}_${big_m}_${sub_name}=${jid_cstacks0}_0"
		eval echo submitted cstacks0 with jobID \$jid_cstacks0_${little_m}_${big_m}_${sub_name}
		done
	done
done


### run cstacks1 ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-cstacks1
logfilesdir=logfiles_cstacks
executable=run-cstacks1.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
			for ((jobIndex = 1; jobIndex <= num_cstacks1jobs; jobIndex++)); #num_cstacks1jobs assigned/calculated previously in part of script that writes out many popmaps
			do
		
			if [ $jobIndex = 1 ]; then job_dependency=$(eval echo \$jid_cstacks0_${little_m}_${big_m}_${sub_name}); else job_dependency=$previous_cstacksID; fi

        	jid_cstacks1=$(sbatch --job-name=$jobname \
				--export=CPUS=$cpus_cstacks1,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,JOBINDEX=$jobIndex,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
             	--cpus-per-task=$cpus_cstacks1 \
             	--mem=$ram_cstacks1 \
             	--array=[0]%1 \
  				--output=./$logfilesdir/${jobname}-${jobIndex}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.out \
  				--error=./$logfilesdir/${jobname}-${jobIndex}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.err \
             	--dependency=afterok:$job_dependency \
             	--kill-on-invalid-dep=yes \
             	--time=$time_cstacks1 \
             	$executable \
             	|cut -d " " -f 4)
             declare "previous_cstacksID=${jid_cstacks1}_0"
             echo submitted cstacks1, iteration $jobIndex
             declare "jid_cstacks1_${little_m}_${big_m}_${sub_name}=$jid_cstacks1"
         	done
     	done
     done
done

echo I am totally done submitting cstacks1 jobs


### run split_cstacks_output_by_sample ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-split_cstacks_output_by_sample
logfilesdir=logfiles_splitintosampsaftercstacks
executable=run-split_cstacks_output_by_sample.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
		jid_splitsamps=$(sbatch --job-name=$jobname \
			--export=CPUS=$cpus_split,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_split \
			--mem=$ram_split \
			--array=[0]%1 \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.err \
			--dependency=afterok:$(eval echo \$jid_cstacks1_${little_m}_${big_m}_${sub_name}) \
			--kill-on-invalid-dep=yes \
			--time=$time_split \
			$executable \
			|cut -d " " -f 4)
		declare "jid_splitsamps_${little_m}_${big_m}_${sub_name}=${jid_splitsamps}_0"
		eval echo submitted splitsamps with jobID \$jid_splitsamps_${little_m}_${big_m}_${sub_name}
		done
	done
done


### run sstacks and tsv2bam ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-sstacksandtsv2bam
logfilesdir=logfiles_sstacksandtsv2bam
executable=run-sstacksandtsv2bam.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}" 
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
		jid_sstacks=$(sbatch --job-name=$jobname \
			--export=N_SAMPS_PER_SSTACKS=$n_samps_per_sstacks,CPUS=$cpus_sstacks,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_sstacks \
			--mem=$ram_sstacks \
			--array=[0-$((n_sstacks_jobs-1))]%$n_sstacks_jobs \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_batch%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_batch%a.err \
			--dependency=afterok:$(eval echo \$jid_splitsamps_${little_m}_${big_m}_${sub_name}) \
			--kill-on-invalid-dep=yes \
			--time=$time_sstacks \
			$executable \
			|awk -v "nID=$n_sstacks_jobs" '{for (i=0; i<nID; i++) printf(":%s",$4"_"i)}')
		declare "jid_sstacks_${little_m}_${big_m}_${sub_name}=$jid_sstacks"
		eval echo submitted sstacks with jobID \$jid_sstacks_${little_m}_${big_m}_${sub_name}
		done
	done
done


### run gstacks ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-gstacks
logfilesdir=logfiles_gstacks
executable=run-gstacks.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}" 
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
		jid_gstacks=$(sbatch --job-name=$jobname \
			--export=CPUS=$cpus_gstacks,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_gstacks \
			--mem=$ram_gstacks \
			--array=[0]%1 \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.err \
			--dependency=afterok$(eval echo \$jid_sstacks_${little_m}_${big_m}_${sub_name}) \
			--kill-on-invalid-dep=yes \
			--time=$time_gstacks \
			$executable \
			|cut -d " " -f 4)
		declare "jid_gstacks_${little_m}_${big_m}_${sub_name}=${jid_gstacks}_0"
		eval echo submitted gstacks with jobID \$jid_gstacks_${little_m}_${big_m}_${sub_name} 
		done
	done
done


### run populations ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-populations
logfilesdir=logfiles_populations
executable=run-populations.sbatch

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


for little_m in "${vector_of_littlem_values[@]}"
do
	for big_m in "${vector_of_bigm_values[@]}"
	do
		for sub_name in "${vector_of_subname_values[@]}"
		do
		jid_populations=$(sbatch --job-name=$jobname \
			--export=CPUS=$cpus_populations,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,BIG_M=$big_m,SUB_NAME=$sub_name,COPY_FILES_TO_EXECUTE_NODE=$copy_files_to_execute_node,PAIRED_END=$paired_end,MODULES_TO_LOAD="${modules_to_load[*]}" \
			--cpus-per-task=$cpus_populations \
			--mem=$ram_populations \
			--array=[0]%1 \
			--output=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.out \
			--error=./$logfilesdir/${jobname}_%A_${run_name}_littlem_${little_m}_bigm_${big_m}_${sub_name}_%a.err \
			--dependency=afterok:$(eval echo \$jid_gstacks_${little_m}_${big_m}_${sub_name}) \
			--kill-on-invalid-dep=yes \
			--time=$time_populations \
			$executable \
			|cut -d " " -f 4)
		declare "all_populations_jids=:${jid_populations}_0${all_populations_jids}"
		done
	done
done
eval echo submitted populations with jobID \$all_populations_jids


### run error check ### -----------------------------------------------------------------------------------------------------------------------

#define variables:
jobname=run-populations #note - this needs to be the same name as in the previous section for this to run after ALL previous jobs finish, for all param values
logfilesdir=logfiles_error_check
executable=run-error_check.sbatch

# set variable to pass to error check
total_little_m_num=$(echo "${#vector_of_littlem_values[@]}")
total_big_M_num=$(echo "${#vector_of_bigm_values[@]}")
total_n_num=$(echo "${#vector_of_subname_values[@]}")

echo using $total_little_m_num values of little m
echo using $total_big_M_num values of big m
echo using $total_n_num relational values of n

#check if logfiles directory has been created in submit dir yet; if not, make one
if [ ! -d ./$logfilesdir ]; then mkdir ./$logfilesdir; fi


	jid_error=$(sbatch --job-name=$jobname \
		--export=N_USTACKS_JOBS=$n_ustacks_jobs,N_SSTACKS_JOBS=$n_sstacks_jobs,CPUS=$cpus_errorcheck,RUN_NAME=$run_name,LITTLE_M=$little_m,MAX_LOC=$max_loc,STORAGENODE=$storagenode,SAMP_STORAGE_DIR=$samp_storage_dir,N_SAMPS=$n_samps,PAIRED_END=$paired_end,BIG_M=$big_m,SUB_NAME=$sub_name,TOTAL_LITTLE_M_NUM=$total_little_m_num,TOTAL_BIG_M_NUM=$total_big_M_num,TOTAL_N_NUM=$total_n_num,NUM_CSTACKS1JOBS=$num_cstacks1jobs \
		--cpus-per-task=$cpus_errorcheck \
		--mem=$ram_errorcheck \
		--array=[0]%1 \
		--output=./$logfilesdir/${jobname}_%A_${run_name}%a.out \
		--error=./$logfilesdir/${jobname}_%A_${run_name}%a.err \
		--dependency=afterany$(eval echo \$all_populations_jids) \
		--kill-on-invalid-dep=yes \
		--time=$time_errorcheck \
		--mail-type=END \
		$executable \
		|cut -d " " -f 4)
	eval echo submitted error check with jobID \$jid_error


