#!/bin/bash

# Shell script to run R script for calcing popgen stats and doing some initial data visualization

#load programs
module purge
module load foss/2020a
module load R/4.0.3

#create output dir if it doesn't already exist
if [ ! -d $OUTDIR ]; then mkdir $OUTDIR; fi
#create output dir to copy all figures into
if [ ! -d $FIGDIR ]; then mkdir $FIGDIR; fi


#set up work environment for job depending on if you want to load files from where they live on cluster or from tmp dir on execute node
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	
	#make unique working directory on execute node
	WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
	if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
	
	#copy input files to tmpdir to speed up run time
	cp -r $INDIR $WORKDIR
	wait
	
	#move to execute node
	cd $WORKDIR
	
	#define input directory for R for this job
	INDIR=$WORKDIR/$(ls)
	
	#define output directory for R for this job
	OUTDIR_FINAL=$OUTDIR
	OUTDIR_SUFFIX=$(echo $OUTDIR | rev | cut -d "/" -f1 | rev)
	OUTDIR=$WORKDIR/$OUTDIR_SUFFIX
	if [ ! -d $OUTDIR ]; then mkdir -p $OUTDIR; fi
	
fi

#copy R scripts to execute node (so it runs faster)
cp $SLURM_SUBMIT_DIR/parsing.R $WORKDIR
cp $SLURM_SUBMIT_DIR/stats.R $WORKDIR
cp $SLURM_SUBMIT_DIR/full_popgen.R $WORKDIR
cp $SLURM_SUBMIT_DIR/comparing_assemblies_and_IDing_r80_params.R $WORKDIR
wait

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#run R, with the name of my R script
Rscript $WORKDIR/full_popgen.R $RUN_NAME $STORAGENODE $INDIR $OUTDIR $NPCS $MINPROPINDIVSSCOREDIN $KEYSDIR $WORKDIR $OUTDIR_FINAL
wait

Rscript $WORKDIR/comparing_assemblies_and_IDing_r80_params.R $RUN_NAME $STORAGENODE $INDIR $OUTDIR $MINPROPINDIVSSCOREDIN $WORKDIR
wait

#print files currently on execute node
echo "files currently where I am working from are"
du -a

#copy all output files back to the output dir on storage node for this dataset,
#and also put a copy of just the figure .pdfs into one master output dir on storage node
#and also put a copy of all genetic related files for r80 params in a new folder
if [ $COPY_FILES_TO_EXECUTE_NODE = yes ]
then
	cp -r $OUTDIR/* $OUTDIR_FINAL
	cp -r $OUTDIR/plots* $FIGDIR
	
	#create dir to store all genetic files associated with r80 optimal params in if it doesn't already exist
	R80DIR=$STORAGENODE/$RUN_NAME/r80_outputs
	if [ ! -d $R80DIR ]; then mkdir $R80DIR; fi
	R80PARAMS=$(less $OUTDIR/r80params.txt)
	
	#cp $INDIR/*${R80PARAMS}* $R80DIR
	cp $OUTDIR/*${R80PARAMS}* $R80DIR
	wait
	
	rm -rf $WORKDIR
fi


#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

