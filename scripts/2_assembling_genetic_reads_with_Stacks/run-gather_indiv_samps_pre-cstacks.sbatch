#!/bin/bash

#this is an array job - runs one job per assembly param combo (for all samples)


#--------------- EXECUTABLE ---------------

# load all of the programs that we want to use
MODULES_TO_LOAD=(${MODULES_TO_LOAD}) #redefine this as an array variable
module purge
module load ${MODULES_TO_LOAD[*]}

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#define file pathss and names
DIR=stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #general name for lowest level dir, where most of work happens

INTARBALL=${RUN_NAME}_${DIR}_sample
OUTTARBALL=${RUN_NAME}_${DIR}.tar.gz #name for output tarball/s to be written to

STORAGEDIR=$STORAGENODE/$RUN_NAME/$DIR #output storage location (where files live when not in active use)



#define more file paths and stage files depending on if you want to copy input files directly to execute node and work on them there
#or leave input files on $STORAGENODE and work on them there
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	WORKDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME/$DIR #input and output file path on execute node
	#make working directories on execute node
	mkdir /tmp/local/$SLURM_JOB_ID/$RUN_NAME
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
fi

if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
then
	WORKDIR=$STORAGEDIR/$UNIQUEJOBID/$DIR #input and output file path on storage node
	rm -rf $STORAGEDIR/working #and remove (now empty) working directory from previous ustack jobs	
	mkdir $STORAGEDIR/$UNIQUEJOBID
	mkdir $WORKDIR
fi



#get all files from previous step for X value of big m (individual samples for each value of big m)
cp -r $STORAGEDIR/$INTARBALL* $WORKDIR
wait
cd $WORKDIR


#untar all individual tarballs using a loop
for f in *.tar.gz; do
  tar xzvf $f &
done
wait

rm ${RUN_NAME}*
wait


#tar up job output
cd ../
tar -czvf $OUTTARBALL $DIR
wait

#count how many files are in tarball we just made
FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)

#print how many files are present
echo I have $FILECOUNT files in tarball
echo I count $N_SAMPS fastq samples in /samples directory on storage node

#check if all files are present (should be number of samples * 3 [.alleles,.snps,.tags per sample] + 1 [the directory])
#if all files are present, copy to storage, otherwise, 
#print error message and environment vars for debugging, requeue the job in a held state
#then exit this script and do not write over storage node contents
if [ $FILECOUNT = $((N_SAMPS*3+1)) ]
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
	cd /tmp/local/
	rm -rf $SLURM_JOB_ID
else
	rm -rf $STORAGEDIR/$UNIQUEJOBID
	wait
fi

echo DONE WITH GATHERING INDIVIDUAL SAMPLES INTO 1 TARBALL PRE-CSTACKS

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
