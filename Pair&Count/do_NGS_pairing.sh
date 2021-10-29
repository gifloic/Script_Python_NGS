#!/bin/bash

## Example of required fastq name format  : Sample-02_S2_L001_1.fastq.gz (read 1)
#                                           Sample-02_S2_L001_2.fastq.gz (read 2)

###   Definitions section   ###
data_dir="Data"
i=0

NbrOfFiles=$(( ${#Gr_eB2H[@]} + ${#Gr_Bevo_v2[@]} + ${#Gr_Bevo_v3[@]} ))
#NbrOfCores="$(grep -c ^processor /proc/cpuinfo -1)"
#HalfNbrOfCores=$(( ${NbrOfCores} / 2 ))
NbrOfCores=50
 
### Preparing the run  ###
echo "Starting run preliminary jobs:"
date
rm -fR $data_dir/gz/Assemblies

# Pair corresponding files
date
echo "Pairing read_1 and read_2 ..."
for entry_1 in `ls $data_dir/gz`; do
    for entry_2 in `ls $data_dir/gz`; do
	sample_1=$(echo $entry_1  | cut -d"." -f1)
	sample_2=$(echo $entry_2  | cut -d"." -f1)

	sample_1_prefix=$(cut -d '_' -f 1 <<< "$sample_1")
	sample_2_prefix=$(cut -d '_' -f 1 <<< "$sample_2")

	sample_1_sufix=${sample_1##*_}
	sample_2_sufix=${sample_2##*_}


#	{
		mkdir $data_dir/gz/Assemblies  > /dev/null 2>&1 # no output
#		mkdir $data_dir/gz/Assemblies/logs > /dev/null 2>&1      # no output; uncomment if log output is required  
#	}
        if [ $sample_1_prefix == $sample_2_prefix ]  && [ $sample_1_sufix -lt $sample_2_sufix ] 
	then
	    echo "Now $entry_1 will be paired with $entry_2"
	    date
#	    time NGmerge -n ${NbrOfCores} -l -v -z $data_dir/gz/Assemblies/logs/$sample_1_prefix.log  -1 $data_dir/gz/$entry_1 -2 $data_dir/gz/$entry_2 -o $data_dir/gz/Assemblies/$sample_1_prefix\_assembled.gz   # uncomment if log output is required    
	    time NGmerge -n ${NbrOfCores} -v -z -1 $data_dir/gz/$entry_1 -2 $data_dir/gz/$entry_2 -o $data_dir/gz/Assemblies/$sample_1_prefix\_assembled.gz  # uncomment if log output is not required
#	    time pear -j ${NbrOfCores} -f $data_dir/gz/$entry_1 -r $data_dir/gz/$entry_2 -o $data_dir/gz/Assemblies/$sample_1_prefix
	    if [ $i -eq 0 ]
	    then
		LIST=${sample_1_prefix}_assembled.gz
		i="-1"
	    else
		LIST=$LIST,${sample_1_prefix}_assembled.gz
	    fi
	fi
    done
done

echo -e "\nJust finished all jobs !!! \n***************************  "
date
