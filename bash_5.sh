#!/bin/bash

echo > Problem_5_output.txt
for file in $(find ../../gene-finder-tool/14_files/ -name 'GCA*fna')
do
    echo $file
    echo $file >> Problem_5_output.txt
    python Problem_5.py $file >> Problem_5_output.txt
done
