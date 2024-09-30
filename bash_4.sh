#!/bin/bash

echo > Problem_4_output.txt
for file in $(find ../../gene-finder-tool/14_files/ -name 'GCA*fna')
do
    echo $file
    echo $file >> Problem_4_output.txt
    python Problem_4.py $file >> Problem_4_output.txt
done
