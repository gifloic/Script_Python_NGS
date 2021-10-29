#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='Extract BC pairs from illumina fastQ. ', prefix_chars='-+')
parser.add_argument("-i", help="the fasta file path + name")
parser.add_argument("-o", help="the output file path + name")

#actually parse command line options
args = parser.parse_args()

if args.i is None: # for command line options
        print('Please provide an option file (-i)')
        quit()

if args.o is None: # for command line options
        print('No output file name provided. Will output as "OutputFasta.fa".')
        OutFFN = 'OutputFasta.fa'
else:
        OutFFN = args.o
        
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
		
# Load FastaDict from fasta file	
def getFastaData (FastaFileName):
	extension = FastaFileName.split(".")[-1]
	if  extension == "gz" or extension == "gzip":
		fasta = get_gzfile_content(FastaFileName)
	else:
		fasta = get_file_content(FastaFileName)
	a = 0
	TotalLines = len(fasta)
	PercLast = 0
	header = ""
	seq = ""
	FastaDict = {}
	FirstSeqLine = True
	for line in fasta:
		if line[:1] == ">":
			header = line.replace(">", "").strip().replace("\n", "")
			FirstSeqLine = True
			a += 1
		else:
			if FirstSeqLine == True:
				seq = line.upper().strip().replace("\n","")
				FirstSeqLine = False
			else:
				seq = "%s%s" % (seq, line.upper().replace("\n","")) 
		FastaDict[header] = seq
		a += 1
		PercNow = int((a * 100)/TotalLines)
		if PercLast != PercNow:
			print("  Reading sequences from file %s ...                    " % (FastaFileName), end="\r", flush=True)
			PercLast = PercNow
	print("  Reading sequences from file %s ... Done!                    " % (FastaFileName))
		
	return FastaDict

#remove non allowed or degenerated characters from sequences
def removeNonAllowDeg(seq):
        seq = seq.upper()
        SeqList = list(seq)
        NewSeqList = []
        for n in range(0, len(SeqList)):
                if SeqList[n] in ['A','C','G','T']:
                        NewSeqList.append(SeqList[n])
        NewSeq = "".join(NewSeqList)

        return NewSeq

# write fasta file from FastaDict
def WriteFastaDict(FastaDict):
        FastaFile=open(OutFFN, "w")
        for header in FastaDict:
                FastaFile.write('>%s\n%s\n' % (header, FastaDict[header]))
        
### MAIN SECTION ###
MyFastaDict = getFastaData(args.i)

for seqname in MyFastaDict.keys():
        NewSeq = removeNonAllowDeg(MyFastaDict[seqname])
        MyFastaDict[seqname] = NewSeq
        
WriteFastaDict(MyFastaDict)

        
