#!/bin/bash

# The goal here is to identify if any errors were thrown and try to pin down where they happened
# This script checks to see that:
# expected logfile directories were created
# expected number of .err and .out files are present in logfiles (depends on number of samples processed)
# expected output directories are present (on storage node)
# expected output files (.tar.gz) are present in output dir/s (number of files depends on number of samples processed)
# checks to see if any errors were printed in .out and .err logfiles
# prints all info/results of this error check script to $logfile


#step1
#step2
#step3
#step4
#step5
#step6
#step7
#step8


#define variables --------------------
step=3		#variable that defines the step of pipeline we want to check for errors
list_of_datasets=testinglist.txt #name of file that contains the name of all of the datasets we want to check
# ----------------------------------------


#create output file for list of datasets we want to check:
logfile=log_summary_${list_of_datasets}_${step}.txt
if [ -f ./$logfile ]
then
	echo "$logfile already exists, removing and creating new one"
	rm ./$logfile
    touch ./$logfile
else 
    touch ./$logfile
fi
        
printf "\n################################################\n Error check for step $step \n################################################\n\n" >> ./$logfile

#count how many datasets we should have processed 
n_datasets=($(cat ../../master_keys/$list_of_datasets | sed '/^\s*$/d' | wc -l))

#define path to logfiles dir (based on step of cleaning pipeline)
if [ $step = 1 ]; then logfile_dir=./logfiles_sampledownload; fi
if [ $step = 2 ]; then logfile_dir=./logfiles_removelowqualreads; fi
if [ $step = 3 ]; then logfile_dir=./logfiles_searchforcommonadapter; fi
if [ $step = 4 ]; then logfile_dir=./logfiles_removingcommonadapter; fi
if [ $step = 5 ]; then logfile_dir=./logfiles_checkcommonadapter; fi
if [ $step = 6 ]; then logfile_dir=./logfiles_calcreadlengths; fi
if [ $step = 7 ]; then logfile_dir=./logfiles_finalreadtrim; fi

#print some stats to summary file (helpful in debugging sometimes)
echo -e "list_of_datasets: $list_of_datasets \n" >> ./$logfile
echo -e "n datasets expected: $n_datasets \n" >> ./$logfile
echo -e "step: $step \n" >> ./$logfile	
echo -e "logfile directory: $logfile_dir \n" >> ./$logfile	



#---------------- logfiles dirs ------------------

# ----- do all of the log files exist in the log directory? ----- 

printf "\n################################################\n Checking log directories\n################################################\n\n" >> ./$logfile

echo -e "##### Searching for missing .err and/or .out logfiles. ##### \n" >> ./$logfile 

# total file number should = number of samples aka $n_samps * 2 (.err and .out)
logfile_total_num=$(( $n_datasets*2 ))

#for each line aka dataset in the list_of_datasets file: 
#       check for errors
#       write out results to master summary log file
  
while read dataset
do 
        
	#define label for each dataset
    run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)
    
    #search for logfiles for each dataset
    logfiles_num=$(ls $logfile_dir | grep $run_name | grep ".err\|.out" | wc -l)
    if [ "$logfiles_num" -lt 2 ]
    then 
    	echo -e "Missing logfile/s for dataset $run_name" >> ./$logfile
    fi

done < ../../master_keys/$list_of_datasets

echo DONE >> ./$logfile


# -----  are there any error messages in log files? ----- 

echo -e '\n##### Searching for errors in .err and .out logfiles. #####' >> ./$logfile

echo -e 'Searching for non-case-sensitive: "error", "cannot stat", "warning", "cannot create", "quota exceeded", and "No such file or directory", "READ COUNTS STILL DO NOT MATCH", "DID NOT". \n' >> ./$logfile


while read dataset
do 
        
	#define label for each dataset
    run_name=$(echo $dataset | cut -d " " -f1 | cut -d "." -f1)
    echo "processesing: $run_name"
	
	#search through .err and .out files with each $run_name in the list in them (not .log files)
	#print any lines that contain the below phrases to log summary file   
	if [ $(grep -riE --include="*$run_name*" --exclude="*.log" 'error|cannot stat|warning|cannot create|quota exceeded|No such file or directory|READ COUNTS STILL DO NOT MATCH|DID NOT' $logfile_dir | grep -v "of allowed errors" | grep -v "\.out\:length" | grep -v "\.out\:WARNING\:" | wc -l) -gt 0 ]
	then
		echo -e "\nDATASET: $run_name" >> ./$logfile
    	echo "ERROR/S FOUND:" >> ./$logfile
    	grep -riE --include="*$run_name*" --exclude="*.log" 'error|cannot stat|warning|cannot create|quota exceeded|No such file or directory|READ COUNTS STILL DO NOT MATCH|DID NOT' $logfile_dir | grep -v "of allowed errors" | grep -v "\.out\:length" | grep -v "\.out\:WARNING\:" >> ./$logfile
	fi

done < ../../master_keys/$list_of_datasets

echo DONE >> ./$logfile



echo DONE