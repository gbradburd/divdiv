#!/bin/bash

#--------------- EXECUTABLE --------------------------------------------------------------

#STEP 2

#this script filters out any low quality reads
#
# input: XXX.XX.gz file of nextgen genetic sequence data - downloaded from NCBI and gzipped, nothing else done to it yet
#			input files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/1_unprocessed_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/1_unprocessed_reads/SRR8625100.fastq.gz
#
# output: XXX.XX.gz file of nextgen genetic sequence data - any/all low quality reads removed
#			output files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/2_lowqualrmvd_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/2_lowqualrmvd_reads/SRR8625100.fq.gz



#-----------------------------------------------------------------------------------------
#--------------- SET UP JOB ENVIRONMENT --------------------------------------------------

#define variables
DATE=`date +%m-%d-%Y_%H.%M.%S` #date/time stamp
WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
UNIQUEJOBID=${SLURM_JOBID} #define unique slurm jobid
logfile=processradtags_readqualityfiltering_${RUN_NAME}_${DATE}.log #name of file to write log file to the process_radtags generates

#load programs
module purge
module load GCC/6.4.0-2.28 OpenMPI/2.1.2
module load powertools
module load Stacks/2.4

#make new storage directory for this step on storage node
if [ ! -d $CLEANREADSDIR ]; then mkdir $CLEANREADSDIR; fi
#remove any files already inside $CLEANREADSDIR (e.g. from previous failed jobs)
#because they won't get written over/removed otherwise and will mess up file counts at end of job and throw error)
if [ $(ls $CLEANREADSDIR | wc -l) -gt 0 ]; then rm $CLEANREADSDIR/*; fi

#make working directories on execute node
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
mkdir $WORKDIR/indir
mkdir $WORKDIR/outdir

cd $WORKDIR

#so that we can check later if we processed all of the files:
#count number of input files on storage node
n_input=($(ls $RAWREADSDIR/*.gz | wc -l))
echo there are $n_input .gz files on storage node to process


#-----------------------------------------------------------------------------------------
#--------------- REMOVE LOW QUALITY READS ------------------------------------------------

#run process_radtags - to search all sequences files for low quality reads and remove any found

### description of cutadapt args used above
#	process_radtags \	
#	-p indir/ \			# input directory
#	-o outdir/ \		# output directory
#	--clean \			# remove any reads with uncalled bases
#	--quality \			# remove low quality reads
#	-w 0.15 \			# use a sliding window of 15% of read length (default) to assess quality, if average quality falls below -s in any window, discard read
#	-s 10 \				# use a phred score of 10, if mean phred score falls below this value in sliding window of length -w, discard read
#	--disable_rad_check	# turn off checking that the restriction enzyme cut site is present and intact


for seqfile in $RAWREADSDIR/*.gz
do	
	
	#copy one sequence file over to execute node (working on one dataset per node aka cluster job)
	echo I am copying seqfile $seqfile to $WORKDIR/indir
	cp $seqfile indir/
	
	seqfile=$(echo $seqfile | rev | cut -d "/" -f 1 | rev)
	
	echo I am running process_radtags on $seqfile
	
	process_radtags \
	-f indir/$seqfile \
	-o outdir/ \
	--clean \
	--quality \
	-w 0.15 \
	-s 15 \
	--disable_rad_check
	
	#rename log file that is automatically written out by process_radtags and append results into it
	cat outdir/process_radtags.indir.log >> $logfile
	
	wait
	
	#delete input seq file to save space
	rm indir/$seqfile
	
	#move processed read file back to storage node
	echo I am copying processed sequence file $seqfile back to $CLEANREADSDIR
	cp outdir/*.gz $CLEANREADSDIR
	wait
	rm outdir/*.gz
	
done

wait

echo I am all done with using process_radtags to remove low quality reads



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#remove the (not renamed version of the) log file from the last iteration of the loop
rm outdir/process_radtags.indir.log

#check if we processed all of the files
#count number of output files on storage node
n_output=($(ls $CLEANREADSDIR/*.gz | wc -l))
echo there are $n_output files in outdir $CLEANREADSDIR
#compare storage node input to storage node output
#requeue job in held state here if file counts don't check out
if [ $n_input = $n_output ]
then 
	echo -e "GOOD - There are $n_input .gz files in starting dir and $n_output .gz files in ending dir"
else 
	echo ERROR! - The number of .gz files in input and output directories on storage node is not equal
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi 

#move cutadapt logfile back to submit node
echo I am copying the log file back to $SLURM_SUBMIT_DIR/$LOGFILESDIR
cp $logfile $SLURM_SUBMIT_DIR/$LOGFILESDIR



#-----------------------------------------------------------------------------------------
#--------------- BUILD SUMMARY FILE OF READ COUNTS PER SAMPLE ----------------------------

#define output file to store results
output_results_file=summary_of_lowqual_read_removal-${RUN_NAME}.csv

echo "File_ID,Retained_Reads,Low_quality,Total" > $output_results_file

#pull out the sample name, total number of reads in the file, number of reads dropped due to low quality, and new total number of reads into a .csv
grep -A 1 '^File' $logfile |grep -v '^File' | grep -v -- "^--$" | cut -f 1-3,6 | tr "\\t" "," >> $output_results_file

wait

echo I AM ALL DONE GENERATING SUMMARY CSV OF LOWQUAL REMOVAL RESULTS



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#check if we processed all of the files

#count number of lines in output results summary file (minus header)
n_outputlines=$(($(less $output_results_file | wc -l)-1))
echo there are $n_outputlines lines in $output_results_file on execute node

#compare
#requeue job in held state here if file counts don't check out
if [ $n_outputlines = $n_output ]
then 
	
	echo -e "I built summary results file for all $n_output .gz files in $CLEANREADSDIR"
	
	#copy summary of lowqual filtering .csv file back to submit node
	cp $output_results_file ${SLURM_SUBMIT_DIR}/$RESULTSFILESDIR

else 
	echo -e "ERROR! - There are $n_output .gz files in $CLEANREADSDIR but $n_outputlines lines in $output_results_file"
	cp $output_results_file ${SLURM_SUBMIT_DIR}/$RESULTSFILESDIR
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi 

echo THIS IS THE END OF THE JOB



#-----------------------------------------------------------------------------------------
#--------------- CLEAN UP ----------------------------------------------------------------

#clean up execute node
cd ../
rm -rf $WORKDIR 

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

