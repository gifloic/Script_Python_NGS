#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import collections
import pandas as pd

# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='Output the table of number of ZMW corresponding to each observed number of passes (Pacbio subreads as input).', prefix_chars='-+')
parser.add_argument("-i", help="the fasta file path + name")
parser.add_argument("-o", help="the output file path + name")
parser.add_argument("-it", help="the input type: fasta or fastq (default = fasta)")

#actually parse command line options
args = parser.parse_args()


if args.i is None: # for command line options
        print('Please provide an option file (-i)')
        quit()
else:
        InFN = args.i

if args.o is None: # for command line options
        OutFN = 'Passes_vs_NbrZMW.tsv'
        print('No output file name provided. Will output as "%s".' % (OutFN))
else:
        OutFN = args.o

if args.it is None: # for command line options
        InType = 'FASTA'
        print('No input type provided. Will assume the input in %s format.' % (InType))
else:
        InType = args.it.upper()
        
# transfer the content of a file to the RAM  
def get_file_content(filename):
        with open(filename) as f:
                data = f.read()
                my_splitlines = data.splitlines()
                f.close()
        return my_splitlines
        
# transfer the content of a file to the RAM  
def get_gzfile_content(filename):

	try:
		with gzip.open(filename, mode='rb') as f:
			data = f.read().decode("utf-8")
			my_splitlines = data.splitlines()
			f.close()
			return my_splitlines
	except FileNotFoundError:
		print("The file %s was not found !" % (filename))
        
### MAIN SECTION ###
# read input file
extension = InFN.split(".")[-1]
if  extension == "gz" or extension == "gzip":
	FastaFile = get_gzfile_content(InFN)
else:
	FastaFile = get_file_content(InFN)

# create a dictionary of number of passes / ZMW
ZMWdict = {}
for line in FastaFile:
        if InType == "FASTQ":
                if line[0] == '@':
                        ZMWName = line.split('/')[1]
                        if ZMWName in ZMWdict.keys():
                                ZMWdict[ZMWName] += 1
                        else:
                                ZMWdict[ZMWName] = 1
        else:
                if line[0] == '>':
                        ZMWName = line.split('/')[1]
                        if ZMWName in ZMWdict.keys():
                                ZMWdict[ZMWName] += 1
                        else:
                                ZMWdict[ZMWName] = 1

# order dictionary based on the number of passes
ZMWdict=collections.OrderedDict(sorted(ZMWdict.items(), key=lambda t: t[1], reverse=False))

# create a dictionary of number of ZMW / number of passes
NbrOfPassesDict = {}
for ZMWName in ZMWdict:
        if ZMWdict[ZMWName] in NbrOfPassesDict:
                NbrOfPassesDict[ZMWdict[ZMWName]] += 1
        else:
                NbrOfPassesDict[ZMWdict[ZMWName]] = 1

#order the dictionary based on number of passes
NbrOfPassesDict=collections.OrderedDict(sorted(NbrOfPassesDict.items(), key=lambda t: t[1], reverse=True))

# convert dictionary to pandas DF and output to file
NbrOfPassesDF = pd.DataFrame(NbrOfPassesDict.items(), columns=['Passes', 'NbrOfZMW'])
NbrOfPassesDF.to_csv(OutFN, index = False, sep="\t", mode='w')

#print the final DF
print(NbrOfPassesDF)

#Exit 
quit()
                

        
