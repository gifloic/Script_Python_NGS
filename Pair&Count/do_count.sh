#!/bin/bash

## Example of required fastq name format  : Sample-02_S2_L001_1.fastq.gz (read 1)
#                                           Sample-02_S2_L001_2.fastq.gz (read 2)

###   Definitions section   ###
data_dir="Data/gz/Assemblies"
i=0
declare -a Motifs=("TAATGAGCCCTTA" "TAATGAGCCCGTC" "TAATGAGCCCGAA" "TAATGAGCCCGCG" "TAATGAGCCCATT")     # Defining group of samples

### Preparing the run  ###
echo "Starting run preliminary jobs:"
date

### Count motifs in fastq.gz samples ###
date
for entry_1 in `ls $data_dir`; do
    for motif in ${Motifs[@]}; do
        echo "Analyzing $data_dir/$entry_1 for motif $motif "
        gunzip -cd $data_dir/$entry_1 | grep -c $motif >> Results.txt
    done
done

echo -e "\nJust finished all jobs !!! \n***************************  "
date
