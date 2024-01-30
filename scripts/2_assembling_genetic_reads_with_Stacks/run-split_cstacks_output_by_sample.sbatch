#!/bin/bash

#this is an array job - runs one job per assembly param combo


#--------------- EXECUTABLE ---------------

# load all of the programs that we want to use
MODULES_TO_LOAD=(${MODULES_TO_LOAD}) #redefine this as an array variable
module purge
module load ${MODULES_TO_LOAD[*]}


#define what value of n should be for cstacks based on/relative to value of big m used in ustacks
if [ $SUB_NAME = "nequalM" ]; then N=$BIG_M; else echo; fi
if [ $SUB_NAME = "nis1lessthanM" ]; then N=$((BIG_M-1)); else echo; fi
if [ $SUB_NAME = "nis1morethanM" ]; then N=$((BIG_M+1)); else echo; fi

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#define file paths
DIR=stacks_littlem_${LITTLE_M}_bigm_${BIG_M}_n${N}_${SUB_NAME} #general name for lowest level dir, where most of work happens

INTARBALL=${RUN_NAME}_${DIR}.tar.gz #name of input tarball/s

STORAGEDIR=$STORAGENODE/$RUN_NAME/stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #output storage location (where files live when not in active use)



#define more file paths and stage files depending on if you want to copy input files directly to execute node and work on them there
#or leave input files on $STORAGENODE and work on them there
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	WORKDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME #input and output file path on execute node
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
fi

if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
then
	WORKDIR=$STORAGEDIR/$UNIQUEJOBID/$DIR #input and output file path on storagenode
	mkdir $STORAGEDIR/$UNIQUEJOBID
	mkdir $WORKDIR
fi


#get tar file from completed cstacks/catalog build step
cp ${STORAGEDIR}/$INTARBALL $WORKDIR
wait
cd $WORKDIR

#untar 
tar -xzvf $INTARBALL
wait
rm $INTARBALL


#split tarball created by cstacks into individual tarballs by sample (and one for catalog)

#build a tarball for the cstacks catalog without including $DIR/ in the path
cd $DIR
echo about to build tarball for catalog
tar -czvf ../${RUN_NAME}_${DIR}_catalog.tar.gz catalog*
wait

#build a tarball for each unique sample 
START=0
END=$((N_SAMPS-1))

for (( sampleIndex=$START; sampleIndex<=$END; sampleIndex++ ))

do
    echo about to build tarball for sample $sampleIndex
    tar -czvf ../${RUN_NAME}_${DIR}_sample$sampleIndex.tar.gz sample$sampleIndex.*
	wait
done


#move to the previous dir
cd -

#count how many tarball files we just made
FILECOUNT=$(ls -lh | grep .tar.gz | wc -l)

#print how many files are present
echo I have $FILECOUNT .tar.gz files in $(pwd)

#check if all files are present (should be 1 file per sample + catalog file)
#if all files are present, copy to storage, otherwise, 
#print error message and environment vars for debugging, requeue the job in a held state
#then exit this script and do not write over storage node contents
if [ $FILECOUNT = $((N_SAMPS+1)) ]
then 
	echo I have all the files and am copying them to storage node
	cp -r *.tar.gz $STORAGEDIR
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

echo DONE BUILDING ONE TARBALL PER SAMPLE POST CSTACKS

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
