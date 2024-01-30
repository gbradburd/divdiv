#!/bin/bash

#--------------- EXECUTABLE --------------------------------------------------------------

#STEP 6

#this script calculates the min, max, median, and average read length per sample, and builds a histogram like table of distribution of read lengths for each sample
#
# input: XXX.fq.gz file of nextgen genetic sequence data - already filtered to remove low quality reads and "true" adapters
#			input files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/3_adapterrmvd_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/3_adapterrmvd_reads/readyforreadlengthcheck_SRR8625100.fq.gz
#
# output: summary_of_readlengths-${RUN_NAME}.csv file of summary stats of min, max, median, mean read length per sample
#			histogram_of_readlengths-${RUN_NAME}.csv of length and count of each read (essentially a histogram in a table)
#			output files are stored in directory at <submitnode>/<user>/<submitdir>/<resultsfilesdir>
#			example: /mnt/home/rhtoczyd/processing_NCBI_data/summary_files_of_results/summary_of_readlengths-bioprj_PRJNA524160_Impatiens-capensis.csv
#					 																 /histogram_of_readlengths-bioprj_PRJNA524160_Impatiens-capensis.csv



#-----------------------------------------------------------------------------------------
#--------------- SET UP JOB ENVIRONMENT --------------------------------------------------

#define variables
DATE=`date +%m-%d-%Y_%H.%M.%S` #date/time stamp
WORKDIR=/tmp/local/$SLURM_JOB_ID #working dir on (remote) execute node
UNIQUEJOBID=${SLURM_JOBID} #define unique slurm jobid
output_results_file=summary_of_readlengths-${RUN_NAME}.csv #define output file to store results
output_histogram_file=histogram_of_readlengths-${RUN_NAME}.csv #define output file to store results

cd $WORKDIR

#so that we can check later if we processed all of the files:
#count number of input files on storage node
n_storage=($(ls $FINALREADSDIR/*.gz | wc -l))
echo there are $n_storage .gz files on storage node to process



#-----------------------------------------------------------------------------------------
#--------------- CALCULATE READ LENGTHS SUMMARY STATS ------------------------------------

#calculate mean, median, min, max read length per individual sample

#create output file
echo -e "sample_name,number_reads,mean_length,median_length,min_length,max_length" > ./$output_results_file

#for each sample .fq.gz file, count the length of each read
for samples in $FINALREADSDIR/*.gz
do {

        # borrowed from: https://unix.stackexchange.com/questions/13731/is-there-a-way-to-get-the-min-max-median-and-average-of-a-list-of-numbers-in/13779#13779
        output=$(zcat ${samples} |\
                awk '{OFS=","}{if(NR%4==2) print length($1)}' |\
                sort -n | awk '
        BEGIN {
        c = 0;
        sum = 0;
        }
        {
        a[c++] = $1;
        sum += $1;
        }
        END {
        ave = sum / c;
        if( (c % 2) == 1 ) {
        median = a[ int(c/2) ];
        } else {
        median = ( a[c/2] + a[c/2-1] ) / 2;
        }
        OFS=",";
        print c, ave, median, a[0], a[c-1];
        }')

        echo -e "$(basename $samples),$output"
		
		wait

} >> ./$output_results_file
    
    wait
        
	echo -e "done with sample $(basename $samples)"

wait 

done

echo DONE COUNTING AND CALCULATING ALL READ LENGTHS



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#check if we processed all of the files

#count number of lines in output results summary file (minus header)
n_output=$(($(less $output_results_file | wc -l)-1))
echo there are $n_output lines in $output_results_file on execute node

#compare
#requeue job in held state here if file counts don't check out
if [ $n_storage = $n_output ]
then 
	echo I counted reads for all $n_storage .gz files in $FINALREADSDIR
	#copy summary of adapter removal .csv file back to submit node
	cp $output_results_file ${SLURM_SUBMIT_DIR}/$RESULTSFILESDIR
else 
	echo -e "ERROR! - There are $n_storage .gz files in $FINALREADSDIR but $n_output lines in $output_results_file"
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi

echo DONE



#-----------------------------------------------------------------------------------------
#--------------- CALCULATE READ LENGTH HISTOGRAM TABLE -----------------------------------

#calculate number of reads of each length for each individual sample

#create output file
echo -e "sample_name,number_reads,read_length" > ./$output_histogram_file

#for each sample .fq.gz file, count the length of each read
for samples in $FINALREADSDIR/*.gz
do {

        #borrowed from: https://unix.stackexchange.com/questions/13731/is-there-a-way-to-get-the-min-max-median-and-average-of-a-list-of-numbers-in/13779#13779
        #generate a table of number of reads per each read length category
        output=$(zcat ${samples} | awk '{OFS="\t"}{if(NR%4==2) print length($1)}' | sort -n | uniq -c)
		wait

		#print each line of $output variable, separating by "," and tacking on sample name		
		while IFS= read -r output
		do
    		number_reads=$(echo $output | cut -d " " -f 1)
			read_length=$(echo $output | cut -d " " -f 2)
			echo -e "$(basename $samples),$number_reads,$read_length"
			wait
		done <<< "$output"
	
		wait

} >> ./$output_histogram_file

		echo -e "done with sample $(basename $samples)"

		wait

done

echo DONE COUNTING AND CALCULATING ALL READ LENGTHS

#copy output file back to submit dir
cp $output_histogram_file ${SLURM_SUBMIT_DIR}/$RESULTSFILESDIR
wait

#-----------------------------------------------------------------------------------------
#--------------- CLEAN UP AND FINISH JOB -------------------------------------------------

#clean up execute node
cd ../
rm -rf $WORKDIR 

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

