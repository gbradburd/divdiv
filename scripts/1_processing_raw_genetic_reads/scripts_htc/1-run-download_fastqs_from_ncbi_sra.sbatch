#!/bin/bash

#--------------- EXECUTABLE --------------------------------------------------------------

#STEP 1

#this script uses SRA ID NCBI sequence numbers listed in a text file to download the sequence files onto storage node of cluster and gzip them (to save space)
#
# input: text file list of NCBI SRA ID numbers (generated by R script prepare_to_download_sequence_data.R)
#			text file has name bioprj_<NCBIbioprj#>_<NCBIspeciesname>.txt 
#			example: - bioprj_PRJNA524160_Impatiens-capensis.txt
#						IDs in text file list are e.g. SRR8625100
#										   			   SRR8625306
#
# output: SRRXXXXXXX.XXXX.gz file of nextgen genetic sequence data 
#			downloaded files are stored in directory at <storagenode>/<user>/bioprj_<NCBIbioprj#>_<NCBIspeciesname>/1_unprocessed_reads
#			example: /mnt/scratch/rhtoczyd/bioprj_PRJNA524160_Impatiens-capensis/1_unprocessed_reads/SRR8625100.fastq.gz



#-----------------------------------------------------------------------------------------
#--------------- SET UP JOB ENVIRONMENT --------------------------------------------------

#load programs we want to use
module purge
module load SRA-Toolkit/2.10.7-centos_linux64
module list

#define variables
WORKDIR=/mnt/scratch/rhtoczyd/tempdirforfasterqdump/ #working dir on (remote) execute node
UNIQUEJOBID=${SLURM_JOBID} #define unique slurm jobid
if [ ! -d $WORKDIR ]; then mkdir $WORKDIR; fi

#copy files to execute node that we need
cp $LISTDIR/$DATASET $WORKDIR
cp $MASTERKEYDIR/$EXPECTEDREADCOUNTLIST $WORKDIR

cd $WORKDIR

#print some info to log for bookkeeping
echo These are the files currently on execute node:
ls
echo Starting to download fastq files from NCBI SRA database
echo Grabbing IDs from $DATASET
eval echo I have $(wc -l $DATASET) ids in the list
echo Saving fastq files to $RAWREADSDIR
echo I am working with ${READ_TYPE}-end reads  



#-----------------------------------------------------------------------------------------
#--------------- DOWNLOAD SEQUENCE FILES -------------------------------------------------


#download each SRA # (RUN_ACC_ID) listed in $DATASET text file
while read RUN_ACC_ID

do 


 #download the sequence file
 readstats=$(fasterq-dump $RUN_ACC_ID --split-files --threads $CPUS --mem 24000MB --bufsize 100000 --curcache 10000 --force --details --temp $WORKDIR --outdir $RAWREADSDIR 2>&1)
 echo readstats are $readstats
 emptyreads=$(echo $readstats | grep -o "reads 0-length.*" | cut -d ":" -f 2 | sed 's/,//g')
 
 #NOTE - the above fasterq-dump line uses the --split-3 option by default, which put paired end read data into 3 files 
 #     - forward and reverse reads into the _1 and _2 fastq files respectively and any singletons (unmatched reads) into the .fastq file
 #
 #	   - to put (all) reads into only two files, forward and reverse, regardless of if a read has a mate in the other file, use
 # 	   - the --split-files option, which does not put singletons in their own file. 
 # BE SURE TO CHANGE fasterq-dump COMMANDS IN ALL 3 PLACES (here and retrys below) IF YOU EDIT this command line
  
 #get the number of lines in file we just downloaded
 if [ $READ_TYPE = SINGLE ]
 then
 downloadedlinecount=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID"* | cut -d " " -f1)
 number_of_lines_per_read=4
 downloadedreadcount=$((downloadedlinecount / number_of_lines_per_read)) 
 fi
 
 if [ $READ_TYPE = PAIRED ]
 then
 downloadedlinecount1=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_1"* | cut -d " " -f1)
 downloadedlinecount2=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_2"* | cut -d " " -f1)
 number_of_lines_per_read=4
 totaldownloadedlines=$((downloadedlinecount1 + downloadedlinecount2))
 totaldownloadedreadspartial=$((totaldownloadedlines / number_of_lines_per_read))
 totaldownloadedreads=$((totaldownloadedreadspartial + emptyreads))
 paired_reads=2
 downloadedreadcount=$((totaldownloadedreads / paired_reads)) 
 fi
 
 echo I just tried to download file for SRA $RUN_ACC_ID
 echo and it has $downloadedreadcount raw reads
 
 #get expected read count according to NCBI metadata
 expectedreadcount=$(grep $RUN_ACC_ID $EXPECTEDREADCOUNTLIST | cut -d ' ' -f 2-)
 echo NCBI says it should have $expectedreadcount
 


 #check if sequence we just downloaded has the number of reads NCBI metadata says it should
 if [ $downloadedreadcount = $expectedreadcount ]
 then 
 	#print message that we have complete file and move on to next one
 	echo READ COUNTS MATCH ON TRY1, Finished downloading fastq file with NCBI run accession number $RUN_ACC_ID
 	echo  
 else 
 		#print message that we don't have complete file and try to download it again
 		echo READ COUNTS DO NOT MATCH, trying to download file again
 		readstats=$(fasterq-dump $RUN_ACC_ID --split-files --threads $CPUS --mem 24000MB --bufsize 100000 --curcache 10000 --force --details --temp $WORKDIR --outdir $RAWREADSDIR 2>&1)
 		echo readstats are $readstats
 		emptyreads=$(echo $readstats | grep -o "reads 0-length.*" | cut -d ":" -f 2 | sed 's/,//g')
 	
 		#check read counts after re-try 1
 		if [ $READ_TYPE = SINGLE ]
 		then
 		downloadedlinecount=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID"* | cut -d " " -f1)
 		number_of_lines_per_read=4
 		downloadedreadcount=$((downloadedlinecount / number_of_lines_per_read)) 
 		fi
 
 		if [ $READ_TYPE = PAIRED ]
 		then
 		downloadedlinecount1=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_1"* | cut -d " " -f1)
 		downloadedlinecount2=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_2"* | cut -d " " -f1)
 		number_of_lines_per_read=4
 		totaldownloadedlines=$((downloadedlinecount1 + downloadedlinecount2))
 		totaldownloadedreadspartial=$((totaldownloadedlines / number_of_lines_per_read))
 		totaldownloadedreads=$((totaldownloadedreadspartial + emptyreads))
 		paired_reads=2
 		downloadedreadcount=$((totaldownloadedreads / paired_reads))
 		fi
 	
 		echo I just tried to download file for SRA $RUN_ACC_ID a second time
 		echo and now it has $downloadedreadcount raw reads
 	
 		#print message about how download re-try 1 went
 		if [ $downloadedreadcount = $expectedreadcount ]
 	then 
 		#print message that we have complete file and move on to next one
 		echo READ COUNTS MATCH ON TRY2, Finished downloading fastq file with NCBI run accession number $RUN_ACC_ID
 		echo  
 	else 
 			#print message that we don't have complete file and try to download it one more time
 			echo READ COUNTS STILL DO NOT MATCH, TRY2, trying to download file one more time
 			readstats=$(fasterq-dump $RUN_ACC_ID --split-files --threads $CPUS --mem 24000MB --bufsize 100000 --curcache 10000 --force --details --temp $WORKDIR --outdir $RAWREADSDIR 2>&1)
 			echo readstats are $readstats
 			emptyreads=$(echo $readstats | grep -o "reads 0-length.*" | cut -d ":" -f 2 | sed 's/,//g')
 		
 			#check read counts after re-try 2
  			if [ $READ_TYPE = SINGLE ]
 			then
 			downloadedlinecount=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID"* | cut -d " " -f1)
 			number_of_lines_per_read=4
 			downloadedreadcount=$((downloadedlinecount / number_of_lines_per_read)) 
 			fi
 
 			if [ $READ_TYPE = PAIRED ]
 			then
 			downloadedlinecount1=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_1"* | cut -d " " -f1)
 			downloadedlinecount2=$(wc -l $RAWREADSDIR/"$RUN_ACC_ID""_2"* | cut -d " " -f1)
 			number_of_lines_per_read=4
 			totaldownloadedlines=$((downloadedlinecount1 + downloadedlinecount2))
 			totaldownloadedreadspartial=$((totaldownloadedlines / number_of_lines_per_read))
 			totaldownloadedreads=$((totaldownloadedreadspartial + emptyreads))
 			paired_reads=2
 			downloadedreadcount=$((totaldownloadedreads / paired_reads))
 			fi
 	
 			echo I just tried to download file for SRA $RUN_ACC_ID a third time
 			echo and now it has $downloadedreadcount raw reads
 		
 			#print message about how download re-try 2 went
 			if [ $downloadedreadcount = $expectedreadcount ]
 		then 
 			#print message that we have complete file and move on to next one
 			echo READ COUNTS MATCH ON TRY3, Finished downloading fastq file with NCBI run accession number $RUN_ACC_ID
 			echo   
 		else 
 			#print message that we failed to get full file after 3 tries and move on to next sample
 	 		echo READ COUNTS STILL DO NOT MATCH, TRY3
 	 		fi
 		fi
 	fi
done < $DATASET

echo  
echo FINISHED TRYING TO DOWNLOADING ALL FASTQ FILES



#-----------------------------------------------------------------------------------------
#--------------- KEEP ONLY FORWARD READS AND CHECK WE GOT ALL FILES ----------------------

#for paired-end datasets, delete the reverse read and/or unpaired read (aka files without a "_" in their name)
#(e.g. SRRXXXXXXX.fastq contains unpaired reads and SRRXXXXXXX_1.fastq and SRRXXXXXXX_2.fastq contain paired reads)
if [ $READ_TYPE = PAIRED ]
then
echo These are the unpaired reads files I am deleting:
find $RAWREADSDIR -type f ! -name '*_*' -print
find $RAWREADSDIR -type f ! -name '*_*' -exec rm -f {} \;
find $RAWREADSDIR -type f -name '*_2*' -print
find $RAWREADSDIR -type f -name '*_2*' -exec rm -f {} \;
fi



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#count how many files are in directory we just downloaded sequence data to
totalfilesdownloaded=$(ls -1 $RAWREADSDIR/ | wc -l)
echo I downloaded $totalfilesdownloaded files to $RAWREADSDIR

#count how many files we expect
#if [ $READ_TYPE = SINGLE ]; then totalfilesexpected=$(wc -l $DATASET | cut -d " " -f1); fi
#if [ $READ_TYPE = PAIRED ]; then totalfilesexpected=$(wc -l $DATASET | cut -d " " -f1); totalfilesexpected=$((totalfilesexpected * 2)); fi
totalfilesexpected=$(wc -l $DATASET | cut -d " " -f1)
echo There are $totalfilesexpected files in the list $DATASET that I asked to download
echo There are $totalfilesdownloaded downloaded files in $RAWREADSDIR

#print message about if we got all of the files or not
#requeue job in held state here if file counts don't check out
if [ $totalfilesdownloaded = $totalfilesexpected ] 
then 
	echo I DOWNLOADED THE NUMBER OF FILES IN THE LIST
else 
	echo I DID NOT DOWNLOAD THE NUMBER OF FILES IN THE LIST
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi

echo  



#-----------------------------------------------------------------------------------------
#--------------- GZIP THE FILES TO SAVE SPACE --------------------------------------------

#gzip all of the files we downloaded, so they take up less space
echo Starting to gzip all files in $RAWREADSDIR

gzip -r $RAWREADSDIR

if [ $? -ne 0 ]
then
	#error was thrown on gzip command, requeue job in held state
	echo ERROR! DURING GZIPPING FILES
  	scontrol requeuehold $UNIQUEJOBID
fi

echo Done gzipping all files in $RAWREADSDIR

echo  



#-----------------------------------------------------------------------------------------
#--------------- CHECK THAT WE PROCESSED ALL FILES ---------------------------------------

#count how many gzipped files are in directory we just downloaded sequence data to
totalfilesgzipped=$(ls $RAWREADSDIR/*.gz | wc -l)
echo I gzipped $totalfilesgzipped files in $RAWREADSDIR

#print message about if we gzipped all of the files or not
#requeue job in held state here if file counts don't check out
if [ $totalfilesgzipped = $totalfilesexpected ] 
then 
	echo I GZIPPED THE NUMBER OF FILES IN THE LIST
else 
	echo I DID NOT GZIP THE NUMBER OF FILES IN THE LIST
	echo ----------------------------------------------------------------------------------------
    echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
	(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)
	scontrol requeuehold $UNIQUEJOBID
	exit 1
fi

echo DONE



#-----------------------------------------------------------------------------------------
#--------------- CLEAN UP ----------------------------------------------------------------

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES:
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)





#graveyard --------------------------------

#download one experiment from NCBI SRA 
#./fasterq-dump SRX4883226 --progress --force --details --outdir ./tempfileshere


#file sizes don't seem to match between what NCBI says they should be and downloaded files, but read counts do, so using those to check downloads
 #get the size of the sequence file we just downloaded
 #downloadedfilesize=$(stat --format=%s "$RUN_ACC_ID"*)
 
 #get expected file size from NCBI
 #expectedfilesize=$(grep $RUN_ACC_ID $EXPECTEDREADCOUNTLIST | cut -d ' ' -f 2-)