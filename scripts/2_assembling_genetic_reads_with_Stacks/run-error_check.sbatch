#!/bin/bash

# The goal here is to identify if any errors were thrown and try to pin down where they happened
# This script checks to see that:
# expected logfile directories were created
# expected number of .err and .out files are present in logfiles (depends on number of samples processed and number of Stacks param values tested)
# expected output directories are present (on storage node)
# expected output files (.tar.gz) are present in output dir/s (number of files depends on number of samples processed and number of Stacks param values tested)
# checks to see if any errors were printed in .out and .err logfiles
# prints all info/results of this error check script to $logfile

# create output file:

logfile=log_summary_${RUN_NAME}.txt

if [ -f ./$logfile ]; then
    echo "$logfile already exists, removing and creating new one"
    rm ./$logfile
    touch ./$logfile
else 
    touch ./$logfile
fi

printf "\n################################################\n Error check for ustacks through populations \n################################################\n\n" >> ./$logfile

#print some stats to summary file (helpful in debugging sometimes)
echo -e "dataset aka run_name is $RUN_NAME \n" >> ./$logfile
echo -e "$N_SAMPS samples were supposed to be processed (variable N_SAMPS) \n" >> ./$logfile
echo -e "$TOTAL_LITTLE_M_NUM values of little m set to be tested (variable TOTAL_LITTLE_M_NUM) \n" >> ./$logfile
echo -e "$TOTAL_BIG_M_NUM values of big m set to be tested (variable TOTAL_BIG_M_NUM) \n" >> ./$logfile

printf "\n################################################\n Checking log directories\n################################################\n\n" >> ./$logfile




#---------------- logfiles dirs ------------------

# do all of the log file directories exist?

echo 'Checking if logfile directories exists.' >> ./$logfile

num_logdirs=$(ls ./ | grep logfiles_ | wc -l)
if [ $num_logdirs == 8 ]
then
	echo -e 'All log directories were created. \n' >> ./$logfile
else
	echo -e "ERROR: $((9-num_logdirs)) logfile directory/s are missing! \n" >> ./$logfile
	echo -e 'The logfile directories that we do have are:' >> ./$logfile
	echo $(ls ./ | grep logfiles_) >> ./$logfile
fi


# search for errors in log files

echo 'Searching for non-case-sensitive: "error", "cannot stat", "warning", "cannot create", "quota exceeded", "unable to open", and "failed to load sample" in .err and .out log files.' >> ./$logfile

if [ $(grep -riE 'error|cannot stat|warning|cannot create|quota exceeded|unable to open|failed to load sample' ./logfiles_* | wc -l) -eq 0 ]
then
    echo -e 'No errors found in log file. \n' >> ./$logfile
else
    echo -e 'ERROR: Errors found in log files. The errors are: \n' >> ./$logfile
    grep -riE 'error|cannot stat|warning|cannot create|quota exceeded|unable to open|failed to load sample' ./logfiles_* >> ./$logfile
fi


# count the number of files in the log directories, to help detect where error/s were

echo -e '\nChecking if all log files are present. \n' >> ./$logfile


# ustacks:
# total file number should = number of samples * number of littlem values * number of bigm values * 2 (.err and .out)
ustacks_num=$(( $N_USTACKS_JOBS*$TOTAL_LITTLE_M_NUM*$TOTAL_BIG_M_NUM*2 ))
echo $ustacks_num log files expected in ustacks logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_ustacks | wc -l) == 1 ]
then
    echo 'ustacks log directory is present'>> ./$logfile
    if [ $(ls ./logfiles_ustacks | wc -l) == $ustacks_num ]
    then
        echo -e 'All ustacks log files are present. \n' >> ./$logfile
    else
        echo -e 'ERROR: Missing ustacks log files! \n' >> ./$logfile
    fi
else
    echo -e 'ERROR: Missing ustacks log directory! \n' >> ./$logfile
fi


# gathertarballs:
# total file number should = number of bigm values * number of littlem values * 2 (.err and .out)
gathertarballs_num=$(( $TOTAL_LITTLE_M_NUM*$TOTAL_BIG_M_NUM*2 ))
echo $gathertarballs_num log files expected in gathertarballs logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_gathertarballs | wc -l) == 1 ]
then
	echo 'gathertarballs log directory is present' >> ./$logfile
	if [ $(ls ./logfiles_gathertarballs | wc -l) == $gathertarballs_num ]
	then
		echo -e 'All gathertarballs log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing gathertarballs log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing gathertarballs log directory! \n' >> ./$logfile
fi


# cstacks:
# total file number should = number of bigm values * nM combos * (number of cstacks1 iterations + 1 for cstacks0) * number of littlem values * 2 (.err and .out)
cstacks_num=$(( $TOTAL_BIG_M_NUM*($NUM_CSTACKS1JOBS+1)*$TOTAL_N_NUM*$TOTAL_LITTLE_M_NUM*2 ))
echo $cstacks_num log files expected in cstacks logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_cstacks | wc -l) == 1 ]
then
	echo 'cstacks log directory is present'>> ./$logfile
	if [ $(ls ./logfiles_cstacks | wc -l) == $cstacks_num ]
	then
		echo -e 'All cstacks log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing cstacks log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing cstacks log directory! \n' >> ./$logfile
fi


# split:
# total file number should = number of bigm values * nM combos * number of littlem values * 2 (.err and .out)
split_num=$(( $TOTAL_BIG_M_NUM*$TOTAL_N_NUM*$TOTAL_LITTLE_M_NUM*2 ))
echo $split_num log files expected in splittarballs logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_splitintosampsaftercstacks | wc -l) == 1 ]
then
	echo 'splittarballs log directory is present'>> ./$logfile
	if [ $(ls ./logfiles_splitintosampsaftercstacks | wc -l) == $split_num ]
	then
		echo -e 'All splittarballs log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing splittarballs log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing splittarballs log directory! \n' >> ./$logfile
fi


# sstacks and tsv2bam:
# total file number should = number of bigm values * nM combos * number of samples * number of littlem values * 2 (.err and .out)
sstacks_num=$(( $TOTAL_BIG_M_NUM*$TOTAL_N_NUM*$N_SSTACKS_JOBS*$TOTAL_LITTLE_M_NUM*2 ))
echo $sstacks_num log files expected in sstacks logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_sstacksandtsv2bam | wc -l) == 1 ]
then
	echo 'sstacksandtsv2bam log directory is present'>> ./$logfile
	if [ $(ls ./logfiles_sstacksandtsv2bam | wc -l) == $sstacks_num ]
	then
		echo -e 'All sstacksandtsv2bam log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing sstacksandtsv2bam log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing sstacksandtsv2bam log directory! \n' >> ./$logfile
fi


# gstacks:
# total file number should = number of bigm values * nM combos * number of littlem values * 2 (.err and .out)
gstacks_num=$(( $TOTAL_BIG_M_NUM*$TOTAL_N_NUM*$TOTAL_LITTLE_M_NUM*2 ))
echo $gstacks_num log files expected in gstacks logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_gstacks | wc -l) == 1 ]
then
	echo 'gstacks log directory is present'>> ./$logfile
	if [ $(ls ./logfiles_gstacks | wc -l) == $gstacks_num ]
	then
		echo -e 'All gstacks log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing gstacks log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing gstacks log directory! \n' >> ./$logfile
fi


# populations:
# number of bigm values * nM combos * number of littlem values * 2 (.err and .out)
populations_num=$(( $TOTAL_BIG_M_NUM*$TOTAL_N_NUM*$TOTAL_LITTLE_M_NUM*2 ))
echo $populations_num log files expected in populations logfile dir >> ./$logfile

if [ $(ls ./ | grep logfiles_populations | wc -l) == 1 ]
then
	echo 'populations log directory is present'>> ./$logfile
	if [ $(ls ./logfiles_populations | wc -l) == $populations_num ]
	then
		echo -e 'All populations log files are present. \n' >> ./$logfile
	else
		echo -e 'ERROR: Missing populations log files! \n' >> ./$logfile
	fi
else
	echo -e 'ERROR: Missing populations log directory! \n' >> ./$logfile
fi




#---------------- output files ------------------

printf "\n################################################\n Checking for presence of output files\n################################################\n\n" >> ./$logfile


# check for the presence of output directories

echo 'Checking if output directories exist.' >> ./$logfile

present_outdirs=$(ls -l $STORAGENODE/$RUN_NAME/ | grep stacks_littlem_ | grep -v .tar.gz | wc -l)
expected_outdirs=$(( $TOTAL_BIG_M_NUM*$TOTAL_LITTLE_M_NUM ))

if [ $present_outdirs == $expected_outdirs ]
then
	echo -e 'All intermediate directories where output tarballs are stored exist. \n' >> ./$logfile
else
	echo -e "ERROR: Missing $(($expected_outdirs-$present_outdirs)) directory/s where output tarballs are stored! Does the output path exist? Are your quotas (file number and/or disk space) high enough? \n" >> ./$logfile
fi



# check for the presence of final output directory (/stacks_final_output_from_populations)

echo 'Checking if final output directory - /stacks_final_output_from_populations - exists.' >> ./$logfile

if [ -d $STORAGENODE/$RUN_NAME/stacks_final_output_from_populations ]
then 
	echo -e 'Final output directory /stacks_final_output_from_populations exists. \n' >> ./$logfile
else
	echo -e "ERROR: Missing /stacks_final_output_from_populations directory! Does the output path exist? Are your quotas (file number and/or disk space) high enough? \n" >> ./$logfile
fi



# check for presence of output files in final output dir (/stacks_final_output_from_populations)

echo 'Checking /stacks_final_output_from_populations for missing files.' >> ./$logfile

# total file number should = number of bigm values * nM combos * number of littlem values
finaltarball_num=$(( $TOTAL_BIG_M_NUM*$TOTAL_N_NUM*$TOTAL_LITTLE_M_NUM ))
echo $finaltarball_num output tarballs expected in stacks_final_output_from_populations dir >> ./$logfile
present_tarballs=$(ls $STORAGENODE/$RUN_NAME/stacks_final_output_from_populations | wc -l)
if [ $present_tarballs == $finaltarball_num ]
then
	echo -e "All final output tarballs (from populations) are present in /stacks_final_output_from_populations. \n" >> ./$logfile
else
    echo -e "ERROR: Missing $(($finaltarball_num-$present_tarballs)) output tarballs in /stacks_final_output_from_populations. \n" >> ./$logfile
fi



# check for the presence of output files in intermediate dirs

cat_num=$(( ($N_SAMPS + 1)* $TOTAL_N_NUM ))

total_num=$(( $N_SAMPS + 1 + $TOTAL_N_NUM + $cat_num ))


for subdir in $(ls -d $STORAGENODE/$RUN_NAME/stacks_littlem_*)
do
    echo "Checking $subdir for missing files" >> ./$logfile    

    if [ $(find $subdir -type f | wc -l) == $total_num ]
    then
        echo -e "All output tarballs are present in $subdir, no problems here. \n" >> ./$logfile
    else
        echo "ERROR: Missing tarball/s in output dir: $subdir" >> ./$logfile
		
		# check for ustacks outputs
        present_tarballs=$(ls $subdir | grep "${RUN_NAME}_stacks_littlem_[0-9]*_bigm_[0-9]*_sample" | wc -l)
        if [ $present_tarballs == $N_SAMPS ]
        then
            echo "All output tarballs from ustacks are present." >> ./$logfile
        else
            echo "ERROR: Missing $(($N_SAMPS-$present_tarballs)) output tarballs from ustacks. Problem likely occurred at this step." >> ./$logfile
        fi

		# check for gather outputs
        if [ $(ls $subdir | grep "${RUN_NAME}_stacks_littlem_[0-9]*_bigm_[0-9]*.tar.gz" | wc -l) == 1 ]
        then
            echo "The 1 tarball from gathertarballs step is present." >> ./$logfile
        else
            echo "ERROR: Missing the output tarball from gathertarballs step. Problem likely occurred here." >> ./$logfile
        fi

        # check for cstacks0/cstacks1/gstacks/populations output files (all of these steps have same output name so file is written over each time)
        # in each subdirectory, should be equal to subname value length
        present_tarballs=$(ls $subdir | grep "${RUN_NAME}_stacks_littlem_[0-9]*_bigm_[0-9]*_n[0-9]*_" | grep -v "sample" | grep -v "catalog" | wc -l)
        if [ $present_tarballs == $TOTAL_N_NUM ]
        then
            echo "All tarballs from cstacks0 and/or cstacks1 and/or gstacks and/or populations are present (output file names are the same for each/all of these steps)." >> ./$logfile
        else
            echo "ERROR: Missing $(($TOTAL_N_NUM-$present_tarballs)) output tarballs from cstacks and/or gstacks and/or populations. Problem likely occurred at one or multiple of these steps (the output file name is the same for each/all of these steps)." >> ./$logfile
        fi

        # check for split/sstacks/tsv2bam output files (all of these steps have same output name so file is written over each time)
        # in each subdirectory, should be equal to  (1 file per sample + catalog) * number of nM comobs, which is cat_num here
        present_tarballs=$(ls $subdir | grep "${RUN_NAME}_stacks_littlem_[0-9]*_bigm_[0-9]*_n[0-9]*_" | grep "sample\|catalog" | wc -l)
        if [ $present_tarballs == $cat_num ]
        then
            echo -e "All tarballs from split and/or sstacks and/or tsv2bam are present (output file names are the same for each/all of these steps). \n" >> ./$logfile
        else
            echo -e "ERROR: Missing $(($cat_num-$present_tarballs)) output tarballs from split and/or sstacks and/or tsv2bam. Problem likely occurred at one or multiple of these steps (the output file name is the same for each/all of these steps). \n" >> ./$logfile
        fi
    fi
done

echo >> ./$logfile
echo "----- Error check is finished -----" >> ./$logfile

#print some environment variables to stdout for records
echo ----------------------------------------------------------------------------------------
echo PRINTING SUBSET OF ENVIRONMENT VARIABLES: 
(set -o posix ; set | grep -v ^_ | grep -v ^EB | grep -v ^BASH | grep -v PATH | grep -v LS_COLORS)

