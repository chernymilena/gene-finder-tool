#!/bin/bash

echo > Problem_6_output.txt
for file in $(find ../../gene-finder-tool/14_files/ -name 'GCA*fna')
do
    echo $file
    echo $file >> Problem_6_output.txt
    python Problem_6.py $file >> Problem_6_output.txt
done
