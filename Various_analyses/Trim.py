#!/usr/bin/env python

#
#  Trim the readsin order to comprise only the region of interest
#
# Oscar Ramos (14/11/2018)
#

#Import required module(s)
import argparse
import datetime
import locale
from os import listdir
from os.path import isfile, join

# American global settings
locale.setlocale(locale.LC_ALL, 'C')

# Global variables section
mypath = None
ANCHOR5 = None
LENGTH = None
ANCHOR3 = None
hscale = None
rc_boolean = None
expected = None
mincounts = None
EnrichCompList = []

# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='Clean, count, rank and calculate the hydrophobicity score of sequences using a fastq file as input. ', prefix_chars='-+')
parser.add_argument("-p", help="the path to the directory containing the fastq files that should be analyzed.")
parser.add_argument("-l", help="the anchor sequence at 5' side (left)")
parser.add_argument("-s", help="the size of the region between anchors to extract (number of consecutive bases)")
parser.add_argument("-r", help="the anchor sequence at 3' side (right)")
parser.add_argument("-rc", help="Use the reverse complement. 0: no; 1: yes.")
parser.add_argument("-op", help="the path to the options file.")

#actually parse command line options
args = parser.parse_args()

if args.op is None: # for command line options
	mypath = args.p
	ANCHOR5 = args.l
	LENGTH = int(args.s)
	ANCHOR3 = args.r
	rc_boolean = int(args.rc)

else: # for options file
	optionsfile=logfile=open(args.op, "r")
	optionlinenbr = 0
	for line in optionsfile :
		if line[:1] != "@" and line!= "\n":
			if optionlinenbr == 1:
				mypath = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 2:
				rc_boolean = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 3:
				ANCHOR5 = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 4:
				LENGTH = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 5:
				ANCHOR3 = (line.rsplit( ":")[ 1 ]).replace("\n", "")
		optionlinenbr += 1

#Dictionaries for the messages
rc_boolean_dic = {0:'No',1:'Yes'}

if rc_boolean != 1:
	rc_boolean = 0

# Begin to log job events
logfile=open("batch_job.log", "w")
logfile.write("***\nUSED OPTIONS\n")
print("***\nUSED OPTIONS")
logfile.write("Directory containing the files to be analyzed : " + mypath + "\n")
print("Directory containing the files to be analyzed : " + mypath)
logfile.write('Create reverse complement : %s \n' % rc_boolean_dic[rc_boolean])
print('Create reverse complement : %s' % rc_boolean_dic[rc_boolean])
logfile.write("5' anchor sequence : " + ANCHOR5 + "\n")
print("5' anchor sequence : " + ANCHOR5)
logfile.write("Size of the region of interest (nbr of bases) : " + str(LENGTH) + "\n")
print("Size of the region of interest (nbr of bases) : " + str(LENGTH))
logfile.write("3' anchor sequence : " + ANCHOR3 + "\n")
print("3' anchor sequence : " + ANCHOR3)

StartTime = datetime.datetime.now()
logfile.write("Job started : " + str(StartTime) + "\n")
print("Job started : " + str(StartTime))
logfile.close()

# transfer the content of a file to the RAM  
def get_file_content(filename):
    with open(filename) as f:
	data = f.read()
	my_splitlines = data.splitlines()
	f.close()
	return my_splitlines

# isolate the dna region of interest
def find_between( text, left, right ):
    try:
	if right != "" and left == "":
		start = 0
        	end = text.index( right, start )
		my_text = text[start:end]
        	return my_text,start,end
	if right == "":
		start = text.index( left ) + len( left )
		end = len(text)
		my_text = text[start:end]
        	return my_text,start,end
	if right != "" and left != "":
		start = text.index( left ) + len( left )
		end = text.index( right, start )
		my_text = text[start:end]
		return my_text,start,end
    except ValueError:
        return "", "", ""

# get the reverse complement of a DNA sequence
def rc_dna(sequence):
	
	rc_sequence = ''
	reverse = my_seq = line[::-1]
	
	# complementary bases dictionary
	complement_nt = { 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}

	for n in range (0, len(reverse)):
		if reverse[n] in complement_nt.keys():
			rc_sequence += complement_nt.get(reverse[n], "N")
	sequence = ''
	return rc_sequence

# Do the job. Get fastq file list and for each fastq file, get sequences, isolate the region of interest between 5' and 3' anchors, 
# translate, clean non-expected sequences and count their occurences. If a hydrophobicity scale is choosen then a score is calculated
# for each sequence using it.
# The output is a fasta file corresponding to the isolated region after translation, an weblogo file and the ranked sequences in csv format 

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for in_file_name in onlyfiles:
	myseq_flag = 0
	i = 0
	header = None
	my_start = None
	my_end = None
	IsMySeqGood = None
	myAllowedSeqs_list = []
	myfastqseqs_list = []
	myallowedseqs_list = []
	myUMIComb_dict = {}
	FileTime = datetime.datetime.now()
	logfile = open("batch_job.log", "a")
	logfile.write("Reading FASTQ file " + in_file_name + " : " + str(datetime.datetime.now()) + "\n")
	print("Reading FASTQ file " + in_file_name + " : " + str(datetime.datetime.now()))
	f1_lines = get_file_content(mypath+in_file_name)
	logfile.write("    Finished reading FASTQ file : " + str(datetime.datetime.now()) + "\n")
	print("    Finished reading FASTQ file : " + str(datetime.datetime.now()))

	#Parse fastq file to retrieve only the region of interest 
	for line in f1_lines:
		if myseq_flag == 3:
			if my_start != "" and my_end != "" and line != "" and IsMySeqGood == 1 :
				my_quality = (line[my_start:my_end] + "\n")
	#			print ("my_start = " + str(my_start) + "; my_end = " + str(my_end) + "\n" + my_quality)
				current_read = (header + "\n" + my_seq + "\n" + my_plus + my_quality)
				myfastqseqs_list.append(current_read)
			my_start = None
			my_end = None
			myseq_flag = 0
			IsMySeqGood = 0
		if myseq_flag == 2:
			my_plus = (line + "\n")
			myseq_flag = 3
		if myseq_flag == 1:
			if rc_boolean == 1:
				my_seq, my_start, my_end = find_between( rc_dna(line), ANCHOR5, ANCHOR3 )
				if len(my_seq) == LENGTH:
					myAllowedSeqs_list.append(my_seq)
				#	my_seq = (header + " \n" + my_seq + "\n")
				#	myfastaseqs_list.append(my_seq)
					IsMySeqGood = 1
			else:
				my_seq, my_start, my_end = find_between( line, ANCHOR5, ANCHOR3 )
				if len(my_seq) == LENGTH:
					myAllowedSeqs_list.append(my_seq)
				#	my_seq = (header + " \n" + my_seq + "\n")
				#	myfastaseqs_list.append(my_seq)
					IsMySeqGood = 1
			myseq_flag = 2
		if myseq_flag == 0 and line[:1] == "@":
				header=line
				myseq_flag = 1
	
# Write fasta file with all allowed sequences
	fastq_out_file = in_file_name.rsplit( ".", 1 )[ 0 ] + "_trimmed.fastq"
	fastqfile = open(fastq_out_file, "w")
	for read in myfastqseqs_list:	
		fastqfile.write(read)
	fastqfile.close()

	logfile.write("    Finished writing " + str(i) + " allowed sequences to fasta file " + fastq_out_file + " : " + str(datetime.datetime.now()) + "\n")
	print("    Finished writing " + str(i) + " allowed sequences to fasta file " + fastq_out_file + " : " + str(datetime.datetime.now()))

# Last log statements
EndTime = datetime.datetime.now()
logfile.write("Batch job ended : " + str(EndTime) + "\n")
print("Batch job ended : " + str(EndTime))

logfile.write("Duration : " + str(EndTime - StartTime) + "\n")
print("Duration : " + str(EndTime - StartTime))

logfile.close()
