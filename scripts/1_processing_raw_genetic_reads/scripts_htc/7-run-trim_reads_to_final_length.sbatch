#!/bin/bash

#--------------- EXECUTABLE --------------------------------------------------------------

#STEP 7

#this script trims all reads to a specified length (and tosses reads that are shorter)
#
# input: XXX.fq.gz file of nextgen genetic sequence data - already filtered to remove low quality reads and "true" adapters
#			input files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/3_adapterrmvd_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/3_adapterrmvd_reads/readyforreadlengthcheck_SRR8625100.fq.gz
#
# output: XXX.fq.gz file of nextgen genetic sequence data - filtered to remove low quality reads, "true" adapters, and all one consistent read length (bases trimmed from R-hand aka 3' side of read)
#			output files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/4_cleanandtrimmedtoconsistentlength_reads
#			example: /mnt/home/rhtoczyd/processing_NCBI_data/summary_files_for_results_of_removingcommonadapter/summary_of_readlengths-bioprj_PRJNA524160_Impatiens-capensis.csv



#-----------------------------------------------------------------------------------------
#--------------- SET UP JOB ENVIRONMENT --------------------------------------------------

#define variables
DATE=`date +%m-%d-%Y_%H.%M.%S` #date/time stamp
WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
UNIQUEJOBID=${SLURM_JOBID} #define unique slurm jobid
output_results_file=summary_of_readcounts_posttrim-${RUN_NAME}.csv #define output file to store results

#load programs
module purge
module load GCC/7.3.0-2.30 
module load OpenMPI/3.1.1
module load powertools
module load cutadapt/2.1-Python-3.6.6 

#make new storage directory for this step on storage node
#move processed read files back to storage node
if [ ! -d $FINALREADSDIR ]; then mkdir $FINALREADSDIR; fi
#remove any files already inside $FINALREADSDIR (e.g. from previous failed jobs)
#because they might not get written over/removed otherwise and will mess up file counts at end of job and throw error
if [ $(ls $FINALREADSDIR | wc -l) -gt 0 ]; then rm $FINALREADSDIR/*; fi


#make working directories on execute node
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
mkdir $WORKDIR/indir
mkdir $WORKDIR/outdir

#copy master file over that lists read length to trim all reads to 
cp master_bookkeeping_sheet-preStacks.csv $WORKDIR

cd $WORKDIR

#so that we can check later if we processed all of the files:
#count number of input files on storage node
n_input=($(ls $INPUTREADSDIR/*.gz | wc -l))
echo there are $n_input .gz files on storage node to process
echo RUN_NAME is $RUN_NAME

echo files in workdir aka $WORKDIR are 
ls

#get length to trim to from the master sheet 
#the second cut is to remove the ^M return character that is present after length value otherwise
#trimlength=$(grep "^$RUN_NAME," master_bookkeeping_sheet-preStacks.csv | cut -d "," -f3 | cut -d$'\r' -f1)
echo Trimming all reads to length $TRIMLENGTH and tossing the rest



#-----------------------------------------------------------------------------------------
#--------------- TRIM READS TO ONE LENGTH ------------------------------------------------

#run cutadapt - to trim all reads to a consistent length and remove reads that are too short

for seqfile in $INPUTREADSDIR/*.gz
do	
	
	#copy one sequence file over to execute node (working on one dataset per node aka cluster job)
	echo I am copying seqfile $seqfile to $WORKDIR/indir
	cp $seqfile indir/
	wait
	
	seqfile=$(echo $seqfile | rev | cut -d "/" -f 1 | rev)
	

	#run cutadapt - to trim all reads to a consistent length and remove reads that are too short
	echo STARTING_SAMPLE
	echo I am running cutadapt on $seqfile
	
	echo SAMPLE: $seqfile
	
	cutadapt \
	--length $TRIMLENGTH \
	--minimum-length $TRIMLENGTH \
	--cores=$CPUS \
	-o outdir/trimmedtolength_${seqfile} \
	indir/$seqfile #\
#	>> cutadapt_commonadapterremoveal_${RUN_NAME}_${DATE}.log
	
	wait
	
	echo " "
	echo " "
	echo I am done processing sample $seqfile
	echo ENDING_SAMPLE

	wait
	
	#copy final .fq file with adapters trimmed back to execute node
	cp outdir/trimmedtolength_${seqfile} $FINALREADSDIR
	
	wait
	
	#remove samples we processed to save disk space
	rm indir/$seqfile
	rm outdir/trimmedtolength_${seqfile}
	
	wait
	
done

wait

echo I am all done with using cutadapt to remove common adapter sequences

echo " "
echo " "
echo " ------ "
echo " "
echo " "



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#check if we processed all of the files
#count the number of lines that have "ENDING_SAMPLE" in the .out log file for this job
#(since we processed seq files one at a time and deleted each when we were done to save space, we can't check by counting those like usual)

# --- FOR THIS JOB --- .out file
file="${SLURM_SUBMIT_DIR}/${LOGFILESDIR}/${JOBNAME}_${RUN_NAME}_${SLURM_JOB_ID}.out"

n_processed=($(grep -c "ENDING_SAMPLE" $file))
echo there are $n_processed lines with ENDING_SAMPLE in them in .out file $file

#compare storage node input to storage node output
#requeue job in held state here if file counts don't check out
if [ $n_input = $n_processed ]
then 
	echo -e "GOOD - There are $n_input .gz files in starting dir and $n_processed ENDING_SAMPLE messages in .out file"
else 
	echo -e "ERROR! - The number of .gz files in input and .out log file is not equal"
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi 


#count number of output files copied back to storage node
n_output=($(ls $FINALREADSDIR/*.gz | wc -l))
wait
echo there are $n_output .gz files in $FINALREADSDIR on storage node

#compare processed to output on storage node
#requeue job in held state here if file counts don't check out
if [ $n_output = $n_processed ]
then 
	echo -e "GOOD - There are $n_output .gz files in output dir on storage node and $n_processed ENDING_SAMPLE messages in .out file"
else 
	echo -e "ERROR! - The number of output .gz files processed and copied back to storage node are not equal"	
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi

#move cutadapt logfile back to submit node
#echo I am copying log file back to $SLURM_SUBMIT_DIR/$LOGFILESDIR
#cp cutadapt_commonadapterremoveal_${RUN_NAME}_${DATE}.log $SLURM_SUBMIT_DIR/$LOGFILESDIR



#-----------------------------------------------------------------------------------------
#--------------- GENERATE SUMMARY FILE OF ADAPTER REMOVAL RESULTS ------------------------

# pull adapter removal stats out of .out log file we just created into a .csv file we can use to visualize them later in R


#write header to a file to store results in
echo "run_name,sample_name,total_reads_processed,reads_that_were_too_short,total_reads_written_to_final_fastq" > $output_results_file


# --- FOR THIS JOB --- .out file
file="${SLURM_SUBMIT_DIR}/${LOGFILESDIR}/${JOBNAME}_${RUN_NAME}_${SLURM_JOB_ID}.out"

echo variable file is $file

cp $file .
wait
#remove variable info at end of file
sed -n '/^------------------/q;p' $file > run.txt	
#get run_name (one per .out file)	
run_name=$(grep ^"RUN_NAME" run.txt | cut -d " " -f 3) 	
echo run_name is $run_name

#get number of samples present in .out we're working on
n_samples=$(grep "STARTING_SAMPLE" run.txt | wc -l | sed 's/ //g')
n_samples=$(expr $n_samples - 1)
echo n_samples is $n_samples

#split .out into one samplefile_XXX per .fq.gz sample
csplit -n 3 -s -f samplefile_ run.txt /"STARTING_SAMPLE"/ {$n_samples}
rm samplefile_000

echo AN ITERATION OF logfile.out

	# --- FOR EACH --- sample in .out
	for sample in samplefile_*
	do
	#get sample_name and overall sample stats
	sample_name=$(grep "SAMPLE:" $sample | sed 's/SAMPLE: //g')
	wait
	total_reads_processed=$(grep "Total reads processed:" $sample | rev | cut -d " " -f 1 | rev | sed 's/,//g' | sed -n '1p')
	wait
	reads_that_were_too_short=$(grep "Reads that were too short:" $sample | rev | cut -d " " -f 2 | rev | sed 's/,//g' | sed -n '1p')
	wait
	total_reads_written_to_final_fastq=$(grep "Reads written (passing filters):" $sample | rev | cut -d " " -f 2 | rev | sed 's/,//g' | sed -n '1p')
	wait
	
	echo $run_name,$sample_name,$total_reads_processed,$reads_that_were_too_short,$total_reads_written_to_final_fastq >> $output_results_file
		
	echo AN ITERATION OF $sample_name
	
	wait 
		
	done
		
	rm samplefile_*

rm run.txt

wait

echo I AM ALL DONE GENERATING SUMMARY CSV OF ADAPTER REMOVAL RESULTS

#copy summary of adapter removal .csv file back to submit node
cp $output_results_file ${SLURM_SUBMIT_DIR}/$RESULTSFILESDIR
wait

echo THIS IS THE END OF THE JOB



#-----------------------------------------------------------------------------------------
#--------------- CLEAN UP AND FINISH JOB -------------------------------------------------

#clean up execute node
cd ../
rm -rf $WORKDIR 

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)



### description of cutadapt args used above
#	cutadapt \
#	---length $TRIMLENGTH \ #length to trim all reads to, "This shortens all reads from input.fastq.gz down to X bases. The removed bases are those on the 3’ end." - manual
#	--minimum-length $TRIMLENGTH \ #length of read to retain read, shorter reads are discarded
#	--cores=$cpus \ #number of cores cutadapt knows it has to use
#	-o outdir/readyforrenaming_${seqfile} \ #output
#	$seqfile \ #input file
#	>> cutadapt_commonadapterremoveal_${RUN_NAME}_${DATE}.log #name of input file > file to direct stats/log info to