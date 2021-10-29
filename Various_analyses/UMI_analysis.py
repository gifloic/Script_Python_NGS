#!/usr/bin/env python

#
# Analyze UMI. 
#
# User can specify the parameters that will used for the run in an option file   
# Here an example that can be used as template after deleting # at the begining of each line
#
#@parameters
#path to the fastq files:../../
#reverse complement[0=no,1=yes]:1
#5'anchor sequence:GTGCA
#size of the region of interest:45
#3'anchor sequence:GGCTCAGAAGCAAAATGGGC
#hydrophobicity scale[0=no hydrophobicity calculation, 1 = Kyte-Doolittle, 2 = Miyazawa-Jernigen]:2
#Allowed protein sequences corresponding to the region of interest (use regular expression):[AS]{2}[QEK]{2}[ILMFVAST]{2}[QEK]{2}[ILMFVAST][ASTGR]{6}
#Minimal counts for each allowed sequence to be considered:2
#@enrichement calculations (enrichement = file_2 frequencies / file_1 frequencies) => specify pairs of files for enrichement calculation => file_1 : file_2 
#01_R0PCR_S10_all_R1_001_cutadapt.fastq:02_R0K7PCR_S14_all_R1_001_cutadapt.fastq
#04_R2PCR_S12_all_R1_001_cutadapt.fastq:05_R2K7PCR_S16_all_R1_001_cutadapt.fastq
#06_R3PCR_S13_all_R1_001_cutadapt.fastq:07_R3K2PCR_S15_all_R1_001_cutadapt.fastq
#06_R3PCR_S13_all_R1_001_cutadapt.fastq:08_R3K7PCR_S17_all_R1_001_cutadapt.fastq
#07_R3K2PCR_S15_all_R1_001_cutadapt.fastq:08_R3K7PCR_S17_all_R1_001_cutadapt.fastq 
#
#
# Oscar Ramos (mer. 16 janv. 2019 09:36:12 CET )
#

#Import required module(s)
import collections
import argparse
import re
import datetime
import locale
import subprocess
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
parser.add_argument("-hs", help="the input hydrophobicity scale to be used. 0: no hydrophobicity calculation; 1: Kyte-Doolittle. J. Mol. Biol. 157:105-132(1982); 2: Miyazawa-Jernigen. Macromolecules 18:534-552(1985)")
parser.add_argument("-rc", help="Use the reverse complement. 0: no; 1: yes.")
parser.add_argument("-allowed", help="Allowed protein sequences corresponding to the region of interest (use regular expression)")
parser.add_argument("-mc", help="Minimal counts for each allowed sequence to be considered")
parser.add_argument("-op", help="the path to the options file.")

#actually parse command line options
args = parser.parse_args()

if args.op is None: # for command line options
	mypath = args.p
	ANCHOR5 = args.l
	LENGTH = int(args.s)
	ANCHOR3 = args.r
	hscale = int(args.hs)
	rc_boolean = int(args.rc)
	expected = args.allowed
	mincounts = args.mc

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
			if optionlinenbr == 7:
				expected = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 8:
				mincounts = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr >= 10:
				files = (line.split(":"))
				EnrichCompList.append([files[0].replace(".fastq", ".csv"), files[1].replace(".fastq\n", ".csv")])
		optionlinenbr += 1

#Dictionaries for the messages
rc_boolean_dic = {0:'No',1:'Yes'}
hscale_dic = {0:'None, all values will be set to 0',1:'Kyte-Doolittle',2:'Miyazawa-Jernigen'}

if rc_boolean != 1:
	rc_boolean = 0

if hscale is None or hscale >= 3:
	print (" \n'Warning: No valid hydrophobic scale was choosen! No hydrophobicity calculation will be made and all values will be set to 0.'\n")
	hscale = 0

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
logfile.write('Hydrophobicity scale : %s \n' % hscale_dic[hscale])
print('Hydrophobicity scale : %s' % hscale_dic[hscale])
print('Mininal number of counts : %i' % mincounts)
if EnrichCompList:
	for pair in EnrichCompList:
		logfile.write("The enrichment will be calculated for: " + pair[0].replace(".csv", ".fastq") + " vs " + pair[1].replace(".csv", ".fastq") + "\n")
		print ("The enrichment will be calculated for: " + pair[0].replace(".csv", ".fastq") + " vs " + pair[1].replace(".csv", ".fastq"))
	logfile.write("***\n")
	print("***")
else:
	logfile.write("Enrichement will not be calculated\n***\n")
	print ("Enrichement will not be calculated\n***")

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
        start = text.index( left ) + len( left )
        end = text.index( right, start )
        return text[start:end]
    except ValueError:
        return ""

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

# translate a DNA sequence to protein sequence
def translate_dna(sequence):
	proteinsequence = ''

	codontable = {
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	    }
	
	for n in range(0,len(sequence),3):
		if sequence[n:n+3] in codontable.keys():
		    proteinsequence += codontable[sequence[n:n+3]]
		else:
		    proteinsequence += "#" # non-expected codon (maybe codon composed of non standard DNA bases A, C, G, T or incomplete)
	sequence = ''
	return proteinsequence

# calculate hydrophobicity using the choosen scale  
def get_hydrophob(seq):
	hydrophob = 0
	hscale_1 = {
		    'A':'1.8', 'C':'2.5', 'D':'-3.5', 'E':'-3.5',
		    'F':'2.8', 'G':'-0.4', 'H':'-3.2', 'I':'4.5',
		    'K':'-3.9', 'L':'3.8', 'M':'1.9', 'N':'-3.5',
		    'P':'-1.6', 'Q':'-3.5', 'R':'-4.5', 'S':'-0.8',
		    'T':'-0.7', 'V':'4.2', 'W':'-0.9', 'Y':'-1.3', '*':'0', '#':'0'
		    }

	hscale_2 = {
		    'A':'5.330', 'C':'7.930', 'D':'3.590', 'E':'3.650',
		    'F':'9.030', 'G':'4.480', 'H':'5.100', 'I':'8.830',
		    'K':'2.950', 'L':'8.470', 'M':'8.950', 'N':'3.710',
		    'P':'3.870', 'Q':'3.870', 'R':'4.180', 'S':'4.090',
		    'T':'4.490', 'V':'7.660', 'W':'5.890', 'Y':'7.630', '*':'0', '#':'0'		    
		    }

	if hscale == 1:
		for char in seq:   
			hydrophob = hydrophob + float(hscale_1[char])
		return hydrophob

	elif hscale == 2:
		for char in seq:   
			hydrophob = hydrophob + float(hscale_2[char])
		return hydrophob

# Do the job. Get fastq file list and for each fastq file, get sequences, isolate the region of interest between 5' and 3' anchors, 
# translate, clean non-expected sequences and count their occurences. If a hydrophobicity scale is choosen then a score is calculated
# for each sequence using it.
# The output is a fasta file corresponding to the isolated region after translation, an weblogo file and the ranked sequences in csv format 

expected_aas = re.compile(r"%s" % expected);
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for in_file_name in onlyfiles:
	myseq_flag = 0
	i = 0
	myAllowedSeqs_list = []
	myfastaseqs_list = []
	myallowedseqs_list = []
	myUMIComb_dict = {}
	FileTime = datetime.datetime.now()
	logfile = open("batch_job.log", "a")
	logfile.write("Reading FASTQ file " + in_file_name + " : " + str(datetime.datetime.now()) + "\n")
	print("Reading FASTQ file " + in_file_name + " : " + str(datetime.datetime.now()))
	f1_lines = get_file_content(mypath+in_file_name)
	logfile.write("    Finished reading FASTQ file : " + str(datetime.datetime.now()) + "\n")
	print("    Finished reading FASTQ file : " + str(datetime.datetime.now()))

	#Parse fastq file to recover only the sequences
	for line in f1_lines:
		if myseq_flag == 1:
			if rc_boolean == 1:
				my_seq = find_between( rc_dna(line), ANCHOR5, ANCHOR3 )
				if len(my_seq) == LENGTH:
					i += 1
					myAllowedSeqs_list.append(my_seq)
					my_seq = (">Allowed read_" + str(i) + " \n" + my_seq + "\n")
					myfastaseqs_list.append(my_seq)
			else:
				my_seq = find_between( line, ANCHOR5, ANCHOR3 )
				if len(my_seq) == LENGTH:
					i += 1
					myAllowedSeqs_list.append(my_seq)
					my_seq = (">Allowed read_" + str(i) + " \n" + my_seq + "\n")
					myfastaseqs_list.append(my_seq)
			myseq_flag = 0
		if myseq_flag == 0 and line[:1] == "@":
			myseq_flag = 1
	
# Write fasta file with all allowed sequences
	fasta_out_file = in_file_name.rsplit( ".", 1 )[ 0 ] + ".fa"
	fastafile = open(fasta_out_file, "w")
	for read in myfastaseqs_list:	
		fastafile.write(read)
	fastafile.close()

	logfile.write("    Finished writing " + str(i) + " allowed sequences to fasta file " + fasta_out_file + " : " + str(datetime.datetime.now()) + "\n")
	print("    Finished writing " + str(i) + " allowed sequences to fasta file " + fasta_out_file + " : " + str(datetime.datetime.now()))

	for read in myAllowedSeqs_list:
		UMIcomb= read[:6] + read[234:240]
		if UMIcomb in myUMIComb_dict.keys():
			myUMIComb_dict[UMIcomb][0] = myUMIComb_dict[UMIcomb][0] + 1
		else:
			myUMIComb_dict[UMIcomb] = ([1]) # key = UMIcomb; values = [counts]
	
	#Delete a sequence if it has less than the minimal number of reads
	myUMIComb_dict = {k: v for k, v in myUMIComb_dict.iteritems() if v[0] >= mincounts}	

	# Rank the sequences based on counts (t[1][1] = second value in the dictionary)
	myUMIComb_dict = collections.OrderedDict(sorted(myUMIComb_dict.items(), key=lambda t: t[1], reverse=True))

	out_file_name = in_file_name.replace(".fastq", ".csv") 
	f2=open(out_file_name, "wt")
	f2.write("counts,UMI_1,UMI_2\n")

	for key, value in myUMIComb_dict.iteritems() :
		f2.write ('%d,%s,%s\n' % (value[0], key[:6], key[6:12]))
	
	f2.close()

	logfile.write("    Finished writing UMI combinations analysis to file " + fasta_out_file + " : " + str(datetime.datetime.now()) + "\n")
	print("    Finished writing UMI combinations analysis to file  " + fasta_out_file + " : " + str(datetime.datetime.now()))


# Last log statements
EndTime = datetime.datetime.now()
logfile.write("Batch job ended : " + str(EndTime) + "\n")
print("Batch job ended : " + str(EndTime))

logfile.write("Duration : " + str(EndTime - StartTime) + "\n")
print("Duration : " + str(EndTime - StartTime))

logfile.close()
