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

INTARBALL=${RUN_NAME}_${DIR}.tar.gz #name of input tarball/s
OUTTARBALL=${RUN_NAME}_${DIR}.tar.gz #name for output tarball/s to be written to

STORAGEDIR=$STORAGENODE/$RUN_NAME/stacks_littlem_${LITTLE_M}_bigm_${BIG_M} #output storage location (where files live when not in active use)

FINALOUTPUTDIR=$STORAGENODE/$RUN_NAME/stacks_final_output_from_populations #one common folder to store all final tarballs/outputs from populations in, so we can easily access all final output files in one place
																	#note that final tarballs are also copied to $STORAGEDIR, like all of the other/previous output tarballs in this pipeline



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


#get output from previous step (gstacks)
cp $STORAGEDIR/$INTARBALL $WORKDIR
wait

#untar
cd $WORKDIR
tar -xzvf $INTARBALL
wait
rm $INTARBALL

#copy popmap over from storage node
cp $STORAGENODE/$RUN_NAME/popmap $DIR

#run populations
populations -P $WORKDIR/$DIR \
-t $CPUS \
-M $WORKDIR/$DIR/popmap \
--verbose \
--log-fst-comp \
--hwe \
--fstats \
--fasta-loci \
--fasta-samples \
--vcf \
--plink \
--phylip \
--phylip-var \
--phylip-var-all \
--fasta-samples-raw

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


#other data filtering options/flags you could use:
#--lnl_lim [float] — filter loci with log likelihood values below this threshold.

#filtering that we do later in R and/or do not want to do because we will lose data/loci we care about for determining best parameters:
#that is, can filter data/loci out in R but cannot get filtered loci back without rerunning populations and turning off or modifying filters
#-p [int] — minimum number of populations a locus must be present in to process a locus.
#-r [float] — minimum percentage of individuals in a population required to process a locus for that population.
#--min_maf [float] — specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5).
#--max_obs_het [float] — specify a maximum observed heterozygosity required to process a nucleotide site at a locus.
#--write_single_snp — restrict data analysis to only the first SNP per locus.
#--write_random_snp — restrict data analysis to one random SNP per locus.

#best practice is usually to not apply above data filters to begin with, filter data manually in R and ID best assembly parameters, 
	#then rerun populations with all filters applied you would like to use in final dataset, then check that filters are working correctly
	#by running final vcf through R scripts again - to confirm everything that should be filtered out is
	#or just rely on R post scripts for filtering (lazier but still good way)

#tar up final output from full stacks pipeline
tar -czvf $OUTTARBALL $DIR
wait

#count how many files are in tarball we just made
FILECOUNT=$(tar -tzf $OUTTARBALL | wc -l)

#print how many files are present
echo I have $FILECOUNT files in $OUTTARBALL tarball

#check if all files are present (should be 	  6 files per sample [.alleles,.matches.tsv.gz,.matches.bam,.snps,.tags,tsv2bam.log] 
#											+ popmap + gstacks.log + gstacks.log.distribs 
#											+ 5 catalog files[.alleles,.calls,.fa,.snps,.tags] 
#											+ directory
#											+ 23 populations.X files [e.g. .haps.vcf, .sumstats.tsv, .snps.vcf, .log])
#if not, print error message and environment vars for debugging, requeue the job in a held state, and exit this script
#if all files are present, copy to storage, otherwise, 
#print error message and environment vars for debugging, requeue the job in a held state
#then exit this script and do not write over storage node contents
if [ $FILECOUNT = $((N_SAMPS*6+3+5+1+23)) ]
then 
	echo I have all the files tarred and am copying them to storage node
	cp $OUTTARBALL $STORAGEDIR 		#copy files to storage dir
	wait
	if [ ! -d $FINALOUTPUTDIR ]; then mkdir $FINALOUTPUTDIR; fi #make one common folder to store output tarballs here too
	cp $OUTTARBALL $FINALOUTPUTDIR 	#also copy files to final_files folder so we can easily access just the final files all in one place from all runs
	wait

else 
	echo ----------------------------------------------------------------------------------------
	echo ERROR! I do not have all the files
	echo Expected $((N_SAMPS*6+3+5+1+23)) files
	echo EXITED WITHOUT TRANSFERRING OUTPUT TARBALL
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

echo DONE WITH POPULATIONS

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
