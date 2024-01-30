#!/bin/bash

#this is an array job - runs one job per chunk of sample/s


#--------------- EXECUTABLE ---------------

# load all of the programs that we want to use
MODULES_TO_LOAD=(${MODULES_TO_LOAD}) #redefine this as an array variable
module purge
module load ${MODULES_TO_LOAD[*]}

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#define file paths
DIR=stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #general name for lowest level dir, where most of work happens
STORAGEDIR=$STORAGENODE/$RUN_NAME/$DIR #output storage location (where files live when not in active use)

#check if the sample directory has been created on /scratch yet
#if not, make one
if [ ! -d $STORAGEDIR ]; then mkdir $STORAGEDIR; fi



#do math to figure out what set of samples we should process in this job based on array ID
START_SAMP_ID=$((SLURM_ARRAY_TASK_ID * N_SAMPS_PER_USTACKS))
END_SAMP_ID=$((START_SAMP_ID + N_SAMPS_PER_USTACKS - 1))
#make sure we run right number of jobs for final array ID
if [ $END_SAMP_ID -ge $N_SAMPS ]; then END_SAMP_ID=$((N_SAMPS - 1)); fi

echo starting with sample${START_SAMP_ID} and ending with sample${END_SAMP_ID}




#for each sample - run ustacks and do all the file tranferring and error checking
for SAMPLE in $(seq $START_SAMP_ID $END_SAMP_ID)
do

	#define output tarball
	OUTTARBALL=${RUN_NAME}_${DIR}_sample${SAMPLE}.tar.gz #name for output tarball/s to be written to

	#define more file paths and stage files depending on if you want to copy input files directly to execute node and work on them there
	#or leave input files on $STORAGENODE and work on them there
	if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
	then
		INDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME/samples #input file path on execute node
		OUTDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME/${DIR}_sample$SAMPLE #output file path on execute node

		#make working directories on execute node
		mkdir /tmp/local/$SLURM_JOB_ID/$RUN_NAME
		mkdir $INDIR
		mkdir $OUTDIR

		#get and stage one input sample tarball from storage node to process in this job via ustacks
		cp ${SAMP_STORAGE_DIR}/sample${SAMPLE}*.fq.gz $INDIR
		wait
	fi
	
	
	if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
	then
		INDIR=$SAMP_STORAGE_DIR #input file path on execute node
		OUTDIR=$STORAGEDIR/working/$UNIQUEJOBID #output file path on execute node
		if [ ! -d $STORAGEDIR/working ]; then mkdir $STORAGEDIR/working; fi
		if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi	
	fi
	
	echo I am processing sample with ID $SAMPLE of $N_SAMPS samples from dataset $RUN_NAME 




	echo Starting ustacks
	#run ustacks
	if [ $PAIRED_END = no ]
	then
	ustacks -t gzfastq \
	-o $OUTDIR \
	-i $SAMPLE \
	-m $LITTLE_M \
	-M $BIG_M \
	--max_locus_stacks $MAX_LOC \
	-p $CPUS \
	-d \
	-f $INDIR/sample${SAMPLE}.fq.gz

	#If exit code ($?) of previous command is not zero, there was an error
	#print error message and environment vars for debugging, requeue the job in a held state, and exit this script
	if [ $? -ne 0 ]
	then
		echo ERROR! There was an error thrown when running Stacks - did job run out of memory?
		echo ----------------------------------------------------------------------------------------
		echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
        	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
		scontrol requeuehold $UNIQUEJOBID
		exit 1
	fi
	fi


	if [ $PAIRED_END = yes ]
	then
	ustacks -t gzfastq \
	-o $OUTDIR \
	-i $SAMPLE \
	-m $LITTLE_M \
	-M $BIG_M \
	--max_locus_stacks $MAX_LOC \
	-p $CPUS \
	-d \
	-f $INDIR/sample${SAMPLE}.1.fq.gz

	#If exit code ($?) of previous command is not zero, there was an error
	#print error message and environment vars for debugging, requeue the job in a held state, and exit this script
	if [ $? -ne 0 ]
	then
		echo ERROR! There was an error thrown when running Stacks - did job run out of memory?
		echo ----------------------------------------------------------------------------------------
		echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
        	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
		scontrol requeuehold $UNIQUEJOBID
		exit 1
	fi
	fi


	#tar up job output
	cd $OUTDIR
	tar -czvf $OUTTARBALL sample${SAMPLE}*
	wait
	
	#count how many files are in tarball we just made
	FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)
	wait
	
	#print how many files are present
	echo I have $FILECOUNT files in tarball

	#check if all files are present (should be 3 [.alleles,.snps,.tags per sample])
	#if all files are present, copy to storage, otherwise,
	#print error message and environment vars for debugging, requeue the job in a held state
	#then exit this script and do not write over storage node contents
	if [ $FILECOUNT = 3 ]
	then 
		echo I have all the files tarred
		cp $OUTTARBALL $STORAGEDIR
		wait
	else 
		echo ----------------------------------------------------------------------------------------
		echo ERROR! I do not have all the files
		echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
		(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
		scontrol requeuehold $UNIQUEJOBID
		exit 1
	fi


	if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
	then
		#clean up execute node
		cd /tmp/local/$SLURM_JOB_ID
		rm -rf $RUN_NAME
		echo DONE WITH USTACKS JOB FOR SAMPLE $SAMPLE
	fi

	
	if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
	then
		cd ../
		rm -rf $UNIQUEJOBID
		echo DONE WITH USTACKS JOB FOR SAMPLE $SAMPLE
	fi


cd ${SLURM_SUBMIT_DIR}


done


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
