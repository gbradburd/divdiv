#!/bin/bash

#--------------- EXECUTABLE --------------------------------------------------------------

#STEP 5

#this script searches for a list of common adapter sequences and removes sequences chunks on the 5' end that match an adapter sequence in the list
#
# NOTE! - this script doesn't bring any files back - this is just to generate a log file of the results so we can confirm that adapters are no longer present
#
# input: XXX.fq.gz file of nextgen genetic sequence data - already filtered to remove low quality reads and "true" adapters
#			input files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/3_adapterrmvd_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/3_adapterrmvd_reads/readyforreadlengthcheck_SRR8625100.fq.gz
#
# output: summary_of_adapter_hits-postadapterremoval-${RUN_NAME}.csv file of counts and lengths of adapters removed (from R-hand end of read aka 3')
#			output files are stored in directory at <submitnode>/<user>/<submitdir>/<resultsfilesdir>
#			example: /mnt/home/rhtoczyd/processing_NCBI_data/summary_files_for_results_of_removingcommonadapter/summary_of_adapter_hits-postadapterremoval-bioprj_PRJNA524160_Impatiens-capensis.csv

#list of all potential adapters downloaded from https://github.com/BioInfoTools/BBMap/tree/master/resources/adapters.fa
#removed all adapters having to do with RNA procedures from this list (based on their names), as we are focusing on DNA seqs.



#-----------------------------------------------------------------------------------------
#--------------- SET UP JOB ENVIRONMENT --------------------------------------------------

#define variables
DATE=`date +%m-%d-%Y_%H.%M.%S` #date/time stamp
WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
UNIQUEJOBID=${SLURM_JOBID} #define unique slurm jobid
output_results_file=summary_of_adapter_hits-postadapterremoval-${RUN_NAME}.csv #define output file to store summary results

#load programs
module purge
module load GCC/7.3.0-2.30 
module load OpenMPI/3.1.1
module load powertools
module load cutadapt/2.1-Python-3.6.6 

#make working directories on execute node
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi
mkdir $WORKDIR/indir
mkdir $WORKDIR/outdir

#copy list of potential adapters to search for over
cp adapters-noRNA-mergedseqs.fa $WORKDIR

cd $WORKDIR

#so that we can check later if we processed all of the files:
#count number of input files on storage node
n_input=($(ls $FINALREADSDIR/*.gz | wc -l))
echo there are $n_input .gz files on storage node to process
echo RUN_NAME is $RUN_NAME



#-----------------------------------------------------------------------------------------
#--------------- SEARCH FOR AND REMOVE ADAPTERS ------------------------------------------

#run cutadapt - to search for just the Illumina common adapter sequence (>PCR_Primer2_rc / PhiX_read1_adapter)

for seqfile in $FINALREADSDIR/*.gz
do	
	
	#copy one sequence file over to execute node (working on one dataset per node aka cluster job)
	echo I am copying seqfile $seqfile to $WORKDIR/indir
	cp $seqfile indir/
	wait
	
	seqfile=$(echo $seqfile | rev | cut -d "/" -f 1 | rev)
	
	
	#run cutadapt - to search for just the Illumina common adapter sequence (>PCR_Primer2_rc / PhiX_read1_adapter)
	echo STARTING_SAMPLE
	echo I am running cutadapt on $seqfile
	
	echo SAMPLE: $seqfile
	
	cutadapt \
	-a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG \
	-a file:adapters-noRNA-mergedseqs.fa \
	-e 0.1 \
	--overlap 4 \
	--cores=$CPUS \
	-o outdir/adaptercheck_${seqfile} \
	indir/$seqfile #\
#	>> cutadapt_commonadapterremoveal_${RUN_NAME}_${DATE}.log
	
	wait
	
	echo " "
	echo " "
	echo I am done processing sample $seqfile
	echo ENDING_SAMPLE

	wait
	
	#remove samples we processed to save disk space
	rm indir/$seqfile
	rm outdir/adaptercheck_${seqfile}
	
	wait
	
	#NOTE - not moving these output adapter trimmed files back to execute node because we 
	#are just checking if any adapters are still present and only need the log files
	
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
	echo ERROR! - The number of .gz files in input and .out log file is not equal
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
echo "run_name,sample_name,total_reads_processed,reads_with_any_adapters,total_basepairs_processed,basepairs_in_any_adapters,adapter_name,reads_with_this_adapter,length_trimmed,times_trimmed,times_expected" > $output_results_file


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
#echo n_samples is $n_samples

#split .out into one samplefile_XXX per .fq.gz sample
csplit -n 3 -s -f samplefile_ run.txt /"STARTING_SAMPLE"/ {$n_samples}
rm samplefile_000

echo AN ITERATION OF .out

	# --- FOR EACH --- sample in .out
	for sample in samplefile_*
	do
	#get sample_name and overall sample stats
	sample_name=$(grep "SAMPLE:" $sample | sed 's/SAMPLE: //g')
	total_reads_processed=$(grep "Total reads processed:" $sample | rev | cut -d " " -f 1 | rev | sed 's/,//g' | sed -n '1p')
	reads_with_any_adapters_ca=$(grep "Reads with adapters:" $sample | rev | cut -d " " -f 2 | rev | sed 's/,//g' | sed -n '1p')
	reads_with_any_adapters_other=$(grep "Reads with adapters:" $sample | rev | cut -d " " -f 2 | rev | sed 's/,//g' | sed -n '2p')
	reads_with_any_adapters=$((reads_with_any_adapters_ca + reads_with_any_adapters_other))
	total_basepairs_processed=$(grep "Total basepairs processed:" $sample | rev | cut -d " " -f 2 | rev | sed 's/,//g' | sed -n '1p') 
	basepairs_in_any_adapters_ca=$(grep "Total written (filtered):" $sample | rev | cut -d " " -f 3 | rev | sed 's/,//g' | sed -n '1p')
	basepairs_in_any_adapters_other=$(grep "Total written (filtered):" $sample | rev | cut -d " " -f 3 | rev | sed 's/,//g' | sed -n '2p')
	basepairs_in_any_adapters=$((basepairs_in_any_adapters_ca + basepairs_in_any_adapters_other))

	#get number of adapters present in samplefile_XXX we're working on
	n_adapters=$(grep "=== Adapter " $sample | wc -l | sed 's/ //g')
	n_adapters=$(expr $n_adapters - 1)
	#echo n_adapters $n_adapters
	
	#split samplefile_XXX file into one file per Adapter
	csplit -n 3 -s -f adapterfile_ $sample /"===\ Adapter"/ {$n_adapters}		
	rm adapterfile_000
	
	echo AN ITERATION OF samplefile_XXX
		
		# --- FOR EACH --- adapter in a samplefile_XXX file
		for adapter in adapterfile_*
		do
		#echo adapter $adapter
		adapter_name=$(grep "=== Adapter " $adapter | sed 's/===//g')
		reads_with_this_adapter=$(sed -n -e 's/^.*Trimmed: //p' $adapter | cut -d " " -f 1)
		
		#get rid of all lines that aren't in table of number and length of adapter matches
		sed -n '/^Overview of removed sequences$/,/^$/p' $adapter > adaptertemp1.txt
		sed -n '/^Overview of removed sequences$/,/^ $/p' adaptertemp1.txt > adaptertemp2.txt
		sed '/^$/d' adaptertemp2.txt > adaptertemp1.txt
		sed '/^ $/d' adaptertemp1.txt > adapter.txt
		
		rm adaptertemp1.txt
		rm adaptertemp2.txt
		rm $adapter
		
		#echo AN ITERATION OF adapterfile_XXX
			
			# --- FOR EACH --- length of adapter trimmed in adapterfile_XXX
			if [ $reads_with_this_adapter = 0 ]
			then
			echo $run_name,$sample_name,$total_reads_processed,$reads_with_any_adapters,$total_basepairs_processed,$basepairs_in_any_adapters,$adapter_name,$reads_with_this_adapter,$length_trimmed,$times_trimmed,$times_expected >> $output_results_file
			else
			sed -e '1,/^length/d' adapter.txt | while read -r trimline;
			do
			length_trimmed=$(echo $trimline | cut -d " " -f 1)
			times_trimmed=$(echo $trimline | cut -d " " -f 2)
			times_expected=$(echo $trimline | cut -d " " -f 3)
			
			echo $run_name,$sample_name,$total_reads_processed,$reads_with_any_adapters,$total_basepairs_processed,$basepairs_in_any_adapters,$adapter_name,$reads_with_this_adapter,$length_trimmed,$times_trimmed,$times_expected >> $output_results_file
			
			done
			fi
			rm adapter.txt
		
		done
		
	done
	rm samplefile_*

rm run.txt

echo I AM ALL DONE GENERATING SUMMARY CSV OF ADAPTER REMOVAL RESULTS

wait

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
#	-a file:adapters-noRNA-mergedseqs.fa \ #fasta files with adapters to look for, the 1 best match \(if any\) is removed, highly similar adapters are merged into consensus seqs. here
#	-e 0.1 \ #error rate - 10%, default, if input sequence is 50 bp long, max of 5 base pairs are allowed to not match adapter seq. and still have adapter seq. removed
#	--overlap 4 \ #must find at least the first 4 matches of an adapter seq. to remove it, "With the default minimum overlap length of 3, only about 0.07 bases are lost per read." - Manual
#	--cores=$cpus \ #number of cores cutadapt knows it has to use
#	-o outdir/readyforrenaming_${seqfile} \ #output
#	$seqfile \ #input file
#	>> cutadapt_commonadapterremoveal_${RUN_NAME}_${DATE}.log #name of input file > file to direct stats/log info to
