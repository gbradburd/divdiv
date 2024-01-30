#!/bin/bash

# Shell script to run R script for calcing pairwise distances (great circle and via sea) between sampled genetic points and GBIF points

#load programs
module purge
module load foss/2020a
module load R/4.0.3
module load rgdal/1.5-16-R-4.0.3

#set up work environment for job depending on if you want to load files from where they live on cluster or from tmp dir on execute node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	#make unique working directory on execute node
	UNIQUEJOBID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
	WORKDIR=/tmp/local/$UNIQUEJOBID #working dir on (remote) execute node
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
	
	#copy input files to tmpdir to speed up run time
	cp ${OUTDIR}/landsamps-gbif_and_genetic-"${RUN_NAME}".csv $WORKDIR
			
	cp $BATHYDIR/world_NOAA_bathy_res4.Robj $WORKDIR
	cp $BATHYDIR/world_NOAA_bathy_res4-resistancematrix-mindepth0.Robj $WORKDIR
	cp $BATHYDIR/world_NOAA_bathy_res4-antimeridian.Robj $WORKDIR
	cp $BATHYDIR/world_NOAA_bathy_res4-resistancematrix-mindepth0-antimeridian.Robj $WORKDIR
	
	#get directory where general inputs live aka one level above where we submitted scripts from
	dropdir=$(echo $SLURM_SUBMIT_DIR | rev | cut -d / -f 1 | rev)
	GENERAL_SUBMIT_DIR=$( echo "$SLURM_SUBMIT_DIR" | sed -e "s/$dropdir$//" )

	cp $GENERAL_SUBMIT_DIR/map_functions.R $WORKDIR
	cp $SLURM_SUBMIT_DIR/move_landsamps_to_water.R $WORKDIR
	
	wait
	
	#move to execute node
	cd $WORKDIR
	
	#define input directory for R for this job
	INDIR=$WORKDIR
	
	#define output directory for R for this job
	OUTDIR_FINAL=$OUTDIR
	OUTDIR_SUFFIX=$(echo $OUTDIR | rev | cut -d "/" -f1 | rev)
	OUTDIR=$WORKDIR/$OUTDIR_SUFFIX
	if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi
		
fi

if [ $COPY_FILES_TO_EXECUTE_NODE = no ]
then
	echo "This job needs to copy files to execute node, change COPY_FILES_TO_EXECUTE_NODE variable to yes and try again"
fi

#print files currently on execute node
echo "files currently where I am working from are"
du -a 

#do math to figure out what set of samples we should process in this job based on array ID
STARTROW=$((SLURM_ARRAY_TASK_ID * N_LANDSAMPS_PER_JOB))
STARTROW=$((STARTROW + 1))
STOPROW=$((STARTROW + N_LANDSAMPS_PER_JOB - 1))
#make sure we run right number of jobs for final array ID
if [ $STOPROW -ge $N_LANDSAMPS ]; then STOPROW=$((N_LANDSAMPS)); fi

echo "I am array job $SLURM_ARRAY_TASK_ID"
echo "I am starting with landsamp row $STARTROW and ending with landsamp row $STOPROW of $N_LANDSAMPS total"
echo "I am running R script directly from execute node"

#run R, with the name of my R script
Rscript $WORKDIR/move_landsamps_to_water.R $RUN_NAME $STORAGENODE $INDIR $OUTDIR $STARTROW $STOPROW

#If exit code ($?) of previous command is not zero, there was an error
#print error message and environment vars for debugging, requeue the job in a held state, and exit this script
if [ $? -ne 0 ]
then
    echo ERROR! There was an error thrown when running the R script - did job run out of memory?
    echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
    (set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
    scontrol requeuehold $SLURM_JOBID
    exit 1
fi


#print files currently on execute node
echo "files currently where I am working from are"
du -a

#copy output files back to storage node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	cp -r $OUTDIR/* $OUTDIR_FINAL
	wait
		
	rm -rf $WORKDIR
	
fi


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

#end