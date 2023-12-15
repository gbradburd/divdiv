#!/bin/bash

for file in lat_long_tables_per_dataset/lat_long_table-*
do

echo file is $file

run_name=$(echo $file | cut -d "/" -f2 | sed -r 's/lat_long_table-//g' | sed -r 's/.txt//g')
echo run_name is $run_name

endpointdir=/mnt/research/bradburd_lab/divdiv_working_popgen/$run_name
if [ ! -d $endpointdir ]; then mkdir $endpointdir; fi

cp $file $endpointdir

done

echo ALL DONE