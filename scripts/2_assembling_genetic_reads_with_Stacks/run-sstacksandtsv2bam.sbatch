#!/bin/bash

#this is an array job - runs one job per chunk of sample/s (within unique assembly param combos)


#--------------- EXECUTABLE ---------------

# load all of the programs that we want to use
MODULES_TO_LOAD=(${MODULES_TO_LOAD}) #redefine this as an array variable
module purge
module load ${MODULES_TO_LOAD[*]}

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#define what value of n should be for cstacks based on/relative to value of big m used in ustacks
if [ $SUB_NAME = "nequalM" ]; then N=$BIG_M; else echo; fi
if [ $SUB_NAME = "nis1lessthanM" ]; then N=$((BIG_M-1)); else echo; fi
if [ $SUB_NAME = "nis1morethanM" ]; then N=$((BIG_M+1)); else echo; fi


#define file paths
DIR=stacks_littlem_${LITTLE_M}_bigm_${BIG_M}_n${N}_${SUB_NAME} #general name for lowest level dir, where most of work happens

STORAGEDIR=$STORAGENODE/$RUN_NAME/stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #output storage location (where files live when not in active use)



#do math to figure out what set of samples we should process in this job based on array ID
START_SAMP_ID=$((SLURM_ARRAY_TASK_ID * N_SAMPS_PER_SSTACKS))
END_SAMP_ID=$((START_SAMP_ID + N_SAMPS_PER_SSTACKS - 1))
#make sure we run right number of jobs for final array ID
if [ $END_SAMP_ID -ge $N_SAMPS ]; then END_SAMP_ID=$((N_SAMPS - 1)); fi



#for each sample - run ustacks and do all the file tranferring and error checking
for SAMPLE in $(seq $START_SAMP_ID $END_SAMP_ID)
do

	#define input and output tarball names
	INTARBALL=${RUN_NAME}_${DIR}_sample$SAMPLE.tar.gz #name of input tarball/s
	OUTTARBALL=${RUN_NAME}_${DIR}_sample$SAMPLE.tar.gz #name for output tarball/s to be written to

	#define more file paths and stage files depending on if you want to copy input files directly to execute node and work on them there
	#or leave input files on $STORAGENODE and work on them there
	if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
	then
		WORKDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME #input and output file path on execute node
	fi
	
	if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
	then
		WORKDIR=$STORAGEDIR/${DIR}_sample$SAMPLE #input and output file path on storagenode	
	fi
	
	#make working directory
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
	
	
	#get cstacks catalog tarball 
	cp $STORAGEDIR/${RUN_NAME}_${DIR}_catalog.tar.gz $WORKDIR
	wait

	#get one sample to process through sstacks and tsv2bam
	cp $STORAGEDIR/$INTARBALL $WORKDIR
	wait
	
	#untar
	cd $WORKDIR
	for f in *.tar.gz; do
  	tar xzvf $f &
	done
	wait
	rm *.tar.gz

	#make output dir to write sstacks output into
	if [ ! -d ${DIR}_sample$SAMPLE ]; then mkdir ${DIR}_sample$SAMPLE; fi
	## note:  mv * ${DIR}_sample$SAMPLE generates an error "cannot move .. to a subdirectory of itself" 
	## instead use the find cmd
	find .  -type f -mindepth 1 -maxdepth 1  -exec mv {} ${DIR}_sample$SAMPLE \;



	#run sstacks
	sstacks -o $WORKDIR/${DIR}_sample$SAMPLE \
	-p $CPUS \
	-c $WORKDIR/${DIR}_sample$SAMPLE \
	-s $WORKDIR/${DIR}_sample$SAMPLE/sample$SAMPLE

	#note: must include a -s for each sample that you want to process aka match to the cstacks catalog

	#If exit code ($?) of previous command is not zero, there was an error
	#print error message and environment vars for debugging, requeue the job in a held state, and exit this script
	if [ $? -ne 0 ]
	then
		echo "ERROR! There was an error thrown when running (sstacks) Stacks - did job run out of memory?"
		echo ----------------------------------------------------------------------------------------
		echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
        	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
		scontrol requeuehold $UNIQUEJOBID
		exit 1
	fi



	#run tsv2bam
	tsv2bam -P $WORKDIR/${DIR}_sample$SAMPLE \
	-t $CPUS \
	-s sample$SAMPLE

	#If exit code ($?) of previous command is not zero, there was an error
	#print error message and environment vars for debugging, requeue the job in a held state, and exit this script
	if [ $? -ne 0 ]
	then
		echo "ERROR! There was an error thrown when running (tsv2bam) Stacks - did job run out of memory?"
		echo ----------------------------------------------------------------------------------------
		echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
       	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
		scontrol requeuehold $UNIQUEJOBID
		exit 1
	fi



	#give tsv2bam.log a unique name (so there aren't clashes when we merge all samples later)
	mv ${DIR}_sample$SAMPLE/tsv2bam.log ${DIR}_sample$SAMPLE/sample${SAMPLE}_tsv2bam.log

	#remove cstacks catalog files so we don't copy them back in the tarball we are about to build
	rm ${DIR}_sample$SAMPLE/catalog*

	#tar up job output (.alleles,.matches,.snps,.tags for sampleX)
	tar -czvf $OUTTARBALL ${DIR}_sample$SAMPLE
	wait

	#count how many files are in tarball we just made
	FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)

	#print how many files are present
	echo I have $FILECOUNT files in $OUTTARBALL tarball

	#check if all files are present (should be 7, .alleles,.matches.tsv.gz,.matches.bam,.snps,.tags,tsv2bam.log for sampleX and directory itself)
	#if all files are present, copy to storage, otherwise, 
	#print error message and environment vars for debugging, requeue the job in a held state
	#then exit this script and do not write over storage node contents
	if [ $FILECOUNT = "7" ]
	then 
		echo I have all the files tarred and am copying them to storage node
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
	else
		rm -rf $WORKDIR
		wait
	fi

	echo DONE WITH SSTACKS AND TSV2BAM FOR SAMPLE $SAMPLE


cd ${SLURM_SUBMIT_DIR}


done


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
