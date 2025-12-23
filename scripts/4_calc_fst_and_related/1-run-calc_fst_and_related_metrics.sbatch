#!/bin/bash

# Shell script to run R script for calcing popgen stats and doing some initial data visualization

#load programs
module purge
module load R/4.2.2-foss-2022b

#create output dir if it doesn't already exist
if [ ! -d $OUTDIR ]; then mkdir -p $OUTDIR; fi


#set up work environment for job depending on if you want to load files from where they live on cluster or from tmp dir on execute node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	#make unique working directory on execute node
	WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
	if [ ! -d $WORKDIR ]; then mkdir -p $WORKDIR; fi
	
	#copy input files to tmpdir to speed up run time
	cp -r $R80DIR/popgenstats* $WORKDIR
	cp -r $R80DIR/gt* $WORKDIR
	
	cp $ALLLATLONG $WORKDIR
	
	cp -r $KEYSDIR/*${RUN_NAME}* $WORKDIR
	wait
	
	#move to execute node
	cd $WORKDIR
	
fi

#copy R scripts to execute node (so it runs faster)
cp $SLURM_SUBMIT_DIR/calc_fst_and_related_metrics.R $WORKDIR
wait

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#run R, with the name of my R script
Rscript $WORKDIR/calc_fst_and_related_metrics.R $RUN_NAME $STORAGENODE $WORKDIR $OUTDIR $MINPROPINDIVSSCOREDIN
wait

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#copy all output files back to the output dir on storage node for this dataset,
#and also put a copy of just the figure .pdfs into one master output dir on storage node
#and also put a copy of all genetic related files for r80 params in a new folder
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	cp -r $WORKDIR/fststats* $R80DIR
	cp -r $WORKDIR/fststats* $OUTDIRALL
	wait
	
	rm -rf $WORKDIR
fi


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

