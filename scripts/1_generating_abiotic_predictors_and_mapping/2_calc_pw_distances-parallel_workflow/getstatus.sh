#!/bin/bash

# pick logfiles dir to search for errors in

logdirtosearch=logfiles_marmapdists




# create output file:

logfile=status_summary-${logdirtosearch}.txt

if [ -f ./$logfile ]; then
    echo "$logfile already exists, removing and creating new one"
    rm ./$logfile
    touch ./$logfile
else 
    touch ./$logfile
fi




# search for last landpoint processed in .out files

echo 'Searching for lines that start with: "done with landpoint" in .out files' >> ./$logfile

for file in ./$logdirtosearch/*.out
do
	echo starting file $file
	prefix=$(echo $file | rev | cut -d/ -f 1 | rev)
	echo $prefix > out.txt
	grep -riE 'done with landpoint' $file > temp.txt
	lastpoint=$(less temp.txt | tail -n1 | cut -d ' '  -f 5)
	totalpoints=$(less temp.txt | tail -n1 | cut -d ' '  -f 7 | cut -d '"' -f 1)
	if [ $lastpoint = $totalpoints ]
	then 
		echo match
	else 
		echo mismatch
		less temp.txt | tail -n1 >> out.txt
		cat out.txt >> ./$logfile
	fi
	rm temp.txt
	rm out.txt
done
