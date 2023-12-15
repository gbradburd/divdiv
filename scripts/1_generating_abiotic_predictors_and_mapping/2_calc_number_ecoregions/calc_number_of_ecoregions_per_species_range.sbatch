#!/bin/bash

#calculate the number of ecoregions that GBIF points for each species fall into

#load modules
module purge
module load foss/2020a
module load R/4.0.3
module load rgdal/1.5-16-R-4.0.3

#make unique working directory on execute node
WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi

#copy input files to execute node
cp $SLURM_SUBMIT_DIR/calc_number_of_ecoregions_per_species_range.R $WORKDIR
cp $SLURM_SUBMIT_DIR/calc_number_of_ecoregions_per_genetic_dataset.R $WORKDIR

cp $GENERALINPUTDIR/countries.outlines.Robj $WORKDIR

cp $GBIFDIR/GBIFLocations_"${SPECIES}".csv $WORKDIR

cp $STORAGENODE/$RUN_NAME"/lat_long_table-"${RUN_NAME}".txt" $WORKDIR

cp -r $GENERALINPUTDIR/Marine_Ecoregions_Of_the_World_\(MEOW\)-shp/ $WORKDIR
cp -r $GENERALINPUTDIR/MarineRealmsShapeFile/ $WORKDIR

wait


#move to execute node
cd $WORKDIR

#run R script - for GBIF points
Rscript calc_number_of_ecoregions_per_species_range.R $WORKDIR "$SPECIES" $RUN_NAME

#run R script - for genetic points
Rscript calc_number_of_ecoregions_per_genetic_dataset.R $WORKDIR "$SPECIES" $RUN_NAME

#copy output files back to storage node
cp -r number_of_ecoregions_* $OUTDIR
wait
	
#copy all output files to one place too (so we don't have to get them from each dataset subfolder)
ALLFILESOUTDIR=$STORAGENODE/ecoregions-output
if [ ! -d $ALLFILESOUTDIR ]; then mkdir $ALLFILESOUTDIR; fi
cp number_of_ecoregions_* $ALLFILESOUTDIR
wait
	
rm -rf $WORKDIR

#print some environment variables to stdout for records
datestamp=$(echo $(date '+%Y-%m-%d'))
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

#end