#!/bin/bash

#this is an array job - runs one job per assembly param combo (for selected samples to add to catalog)


#--------------- EXECUTABLE ---------------

# load all of the programs that we want to use
MODULES_TO_LOAD=(${MODULES_TO_LOAD}) #redefine this as an array variable
module purge
module load ${MODULES_TO_LOAD[*]}


#define what value of n should be for cstacks based on/relative to value of big m used in ustacks
if [ $SUB_NAME = "nequalM" ]; then N=$BIG_M; else echo; fi
if [ $SUB_NAME = "nis1lessthanM" ]; then N=$((BIG_M-1)); else echo; fi
if [ $SUB_NAME = "nis1morethanM" ]; then N=$((BIG_M+1)); else echo; fi
echo Running cstacks with n equal to $N using dataset with big M from ustacks equal to $BIG_M

#define unique slurm jobid
UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

#define file paths
DIR=stacks_littlem_${LITTLE_M}_bigm_${BIG_M}_n${N}_${SUB_NAME} #general name for lowest level dir, where most of work happens

INTARBALL=${RUN_NAME}_stacks_littlem_${LITTLE_M}_bigm_${BIG_M}.tar.gz #name of input tarball/s
OUTTARBALL=${RUN_NAME}_${DIR}.tar.gz #name for output tarball/s to be written to

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



#get files from previous step (now grouped together at level of unique littlem/bigm combos)
cp -r ${STORAGEDIR}/$INTARBALL $WORKDIR
wait
cd $WORKDIR

#untar
tar -xzvf $INTARBALL
wait
rm $INTARBALL

#rename working dir (this will be input and output dir)
mv stacks_littlem_${LITTLE_M}_bigm_${BIG_M} $DIR

#get popmap file that specifies which samples to use to start building the catalog of loci
cp $STORAGENODE/$RUN_NAME/popmap.0 $DIR 

echo Running cstacks with n equal to $N using dataset $RUN_NAME with big M from ustacks equal to $BIG_M
echo "These are the samples that I have in the popmap (hard coded to be popmap.0)"
cat $DIR/popmap.0
printf "\nend of popmap\n"

#run cstacks
cstacks -n $N \
-p $CPUS \
--in_path $WORKDIR/$DIR/ \
--popmap $WORKDIR/$DIR/popmap.0

#note: include an -s flag for every sample that you want to add to the catalog OR --popmap with list of samples to add
# -s $WORKDIR/$DIR/sample0 (and remove --in_path flag and --popmap flag)

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

#remove popmap - we don't want it to be included in the output tarball
rm $WORKDIR/$DIR/popmap.0

#tar up job output
tar -czvf $OUTTARBALL $DIR
wait

#count how many files are in tarball we just made
FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)

#print how many files are present
echo I have $FILECOUNT files in tarball

#check if all files are present (should be number of samples * 3 [.alleles,.snps,.tags per sample] + 1 [the directory] + 3 [catalog.alleles,catalog.snps,catalog.tags from cstacks catalog])
#if all files are present, copy to storage, otherwise, 
#print error message and environment vars for debugging, requeue the job in a held state
#then exit this script and do not write over storage node contents
if [ $FILECOUNT = $((N_SAMPS*3+1+3)) ] 
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

echo DONE WITH CSTACKS0

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
