#!/bin/bash

#this is an array job - runs one job per sample (within unique assembly param combos)


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

OUTTARBALL=${RUN_NAME}_${DIR}.tar.gz #name for output tarball/s to be written to

STORAGEDIR=$STORAGENODE/$RUN_NAME/stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #output storage location (where files live when not in active use)



#define more file paths and stage files depending on if you want to copy input files directly to execute node and work on them there
#or leave input files on $STORAGENODE and work on them there
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	WORKDIR=/tmp/local/$SLURM_JOB_ID/$RUN_NAME #input and output file path on execute node

else
	WORKDIR=$STORAGEDIR/$DIR #input and output file path on storagenode
fi

#make working directories on execute node
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi

#get cstacks catalog tarball 
cp $STORAGEDIR/${RUN_NAME}_${DIR}_catalog.tar.gz $WORKDIR
wait

#get all samples from a unique combination of littlem/bigm/n to process through gstacks
cp $STORAGEDIR/${RUN_NAME}_${DIR}_sample*.tar.gz $WORKDIR
wait

#untar
cd $WORKDIR
for f in *.tar.gz; do
  tar xzvf $f &
done
wait
rm *.tar.gz

#move all sample files out of their individual sampleX folders to this dir level
find . -mindepth 2 -type f -exec mv -i -- {} . \;
#and delete all the empty sampleX directories that are left over
find . -depth -mindepth 1 -type d -empty -exec rmdir {} \;

#make working dir for all untarred samples and move all sampleX files, catalog, and popmap into it
mkdir $DIR
## note:  mv * $DIR generates an error "cannot move .. to a subdirectory of itself" 
## instead use the find cmd
find .  -type f -mindepth 1 -maxdepth 1  -exec mv {} $DIR \;

#copy popmap over from storage node
cp $STORAGENODE/$RUN_NAME/popmap $DIR


#run gstacks
gstacks -P $WORKDIR/$DIR \
-t $CPUS \
-M $WORKDIR/$DIR/popmap

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


#tar up job output (.alleles,.matches,.snps,.tags for sampleX)
tar -czvf $OUTTARBALL $DIR
wait

#count how many files are in tarball we just made
FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)

#print how many files are present
echo I have $FILECOUNT files in $OUTTARBALL tarball

#check if all files are present (should be 6 files per sample [.alleles,.matches.tsv.gz,.matches.bam,.snps,.tags,tsv2bam.log] 
#											+ popmap + gstacks.log + gstacks.log.distribs + 5 catalog files[.alleles,.calls,.fa,.snps,.tags] 
#											+ directory)
#if all files are present, copy to storage, otherwise, 
#print error message and environment vars for debugging, requeue the job in a held state
#then exit this script and do not write over storage node contents
if [ $FILECOUNT = $((N_SAMPS*6+3+5+1)) ]
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
	cd /tmp/local/
	rm -rf $SLURM_JOB_ID
else
	rm -rf $WORKDIR
	wait
fi

echo DONE WITH GSTACKS

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
