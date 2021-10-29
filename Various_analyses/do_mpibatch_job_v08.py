#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Clean/demultiplex, count, calculate the hydrophobicity score and frequency
# of sequences and rank using a fastq file as input. I also create a fasta file
# comprising the corresponding protein sequences and generates a weblogo. 
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
# This algorithm takes only is ~2.65% of the time of the previous version concerning the enrichment calculation
#
# Oscar Ramos (21/03/2018)
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

# MPI=> module import
from mpi4py import MPI

# Use default communicator. No need to complicate things.
COMM = MPI.COMM_WORLD   # get MPI communicator object

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
	optionsfile=open(args.op, "r")
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
			if optionlinenbr == 6:
				hscale = int(line.rsplit( ":")[ 1 ])
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
logfile.write('Allowed protein sequence pattern : %s\n' % expected)
print('Allowed protein sequence pattern : %s' % expected)
logfile.write('Mininal number of counts : %i\n' % mincounts)
print('Mininal number of counts : %i' % mincounts)
if EnrichCompList:
	for pair in EnrichCompList:
		logfile.write("The enrichment will be calculated for: " + pair[0].replace(".csv", ".fastq") + " => " + pair[1].replace(".csv", ".fastq") + "\n")
		print ("The enrichment will be calculated for: " + pair[0].replace(".csv", ".fastq") + " => " + pair[1].replace(".csv", ".fastq"))
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

# MPI=>function required to split jobs
def split(container, count):
    """
    Simple function splitting a container into equal length chunks.

    Order is not preserved but this is potentially an advantage depending on
    the use case.
    """
    return [container[_i::count] for _i in range(count)]



# Do the job. Get fastq file list and for each fastq file, get sequences, isolate the region of interest between 5' and 3' anchors, 
# translate, clean non-expected sequences and count their occurences. If a hydrophobicity scale is choosen then a score is calculated
# for each sequence using it.
# The output is a fasta file corresponding to the isolated region after translation, an weblogo file and the ranked sequences in csv format 

expected_aas = re.compile(r"%s" % expected);
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

## MPI related block starts here

# MPI=>Master define chuncks
if COMM.rank == 0:
	# Split into however many cores are available.
	jobs = split(onlyfiles, COMM.size)
else:
	jobs = None

# MPI=> Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)
results = []

for job in jobs:
	in_file_name = job
	ranking = {}
	myseq_flag = 0
	i = 0
	myseqs_list = []
	myallowedseqs_list = []
	FileTime = datetime.datetime.now()
	logfile = open("batch_job.log", "a")
	logfile.write("Reading FASTQ file " + in_file_name + " from " + MPI.Get_processor_name() + ", task " + str(COMM.rank+1) + ". Current date and time is : " + str(datetime.datetime.now()) + "\n")
	print("Reading FASTQ file " + in_file_name + " from " + MPI.Get_processor_name() + ", task " + str(COMM.rank+1) + ". Current date and time is : " + str(datetime.datetime.now()))
	f1_lines = get_file_content(mypath+in_file_name)
	logfile.write("    Finished reading FASTQ file : " + str(datetime.datetime.now()) + "\n")
	print("    Finished reading FASTQ file : " + str(datetime.datetime.now()))

	#Parse fastq file to recover only the sequences
	for line in f1_lines:
		if myseq_flag == 1:
			if rc_boolean == 1:
				my_seq = rc_dna(line)
			else:
				my_seq = line
			myseqs_list.append(my_seq)
			myseq_flag = 0
		if myseq_flag == 0 and line[:1] == "@":
			myseq_flag = 1

	del f1_lines # release fastq file from memory

	global_counts = len(myseqs_list)

	logfile.write("    Finished parsing FASTQ sequences : " + str(datetime.datetime.now()) + "\n")
	print("    Finished parsing FASTQ sequences : " + str(datetime.datetime.now()))

	# Count reads
	for read in myseqs_list:
		my_seq = find_between( read, ANCHOR5, ANCHOR3 )
		if len(my_seq) == LENGTH:
			my_protseq = translate_dna(my_seq)
			if expected_aas.search(my_protseq):
				i += 1
				myread = ("> Allowed_read_" + str(i) + " \n" + my_protseq + "\n")
				myallowedseqs_list.append(myread)
				if my_seq in ranking.keys():
					ranking[my_seq][1] = ranking[my_seq][1] + 1
				else:
					ranking[my_seq] = ([my_protseq, 1, 0, 0]) # key = my_seq; values = [my_protseq, counts, frequency, hydrophobicity]
	
	del myseqs_list

	# Write fasta file with all allowed sequences
	fasta_out_file = in_file_name.rsplit( ".", 1 )[ 0 ] + ".fa"
	fastafile = open(fasta_out_file, "w")

	for reads in myallowedseqs_list:
		fastafile.write(reads)

	# Count the total number of allowed reads
	allowedseqs_counts = sum(x[1] for x in ranking.values())

	logfile.write("    Finished counting, translation and analysis of %i reads (%i retained for further comparisons) : " % (global_counts, allowedseqs_counts) + str(datetime.datetime.now()) + "\n")		
	print("    Finished counting, translation and analysis of %i reads (%i retained for further comparisons) : " % (global_counts, allowedseqs_counts) + str(datetime.datetime.now()))

	fastafile.close()

	# Call weblogo to produce the output file
	subprocess.call(["weblogo", "-f", fasta_out_file, "-D", "fasta", "-A", "protein", "-F", "pdf", "-U", "probability", "-o", fasta_out_file.replace(".fa", ".pdf")])

	logfile.write("    Finished the creation of a weblogo : " + str(datetime.datetime.now()) + "\n")
	print("    Finished the creation of a weblogo : " + str(datetime.datetime.now()))

	#Delete a sequence if it has less than the minimal number of reads
	ranking = {k: v for k, v in ranking.iteritems() if v[1] >= mincounts}	

	# Rank the sequences based on counts (t[1][1] = second value in the dictionary)
	ranking=collections.OrderedDict(sorted(ranking.items(), key=lambda t: t[1][1], reverse=True))				

	logfile.write("    Finished ranking %i allowed reads (%i non-redundant sequences) : " % (allowedseqs_counts, len(ranking)) + str(datetime.datetime.now()) + "\n")
	print("    Finished ranking %i allowed reads (%i non-redundant sequences) : " % (allowedseqs_counts, len(ranking)) + str(datetime.datetime.now()))

	# open the output file, write the heading + ranked sequences
	out_file_name = in_file_name.replace(".fastq", ".csv") 
	f2=open(out_file_name, "wt")
	f2.write("frequency,counts,aa_sequence,nt_sequence,hydrophobicity\n")

	if hscale == 0:
		for key, value in ranking.iteritems() :
			ranking[key][2] = value[1]/float(global_counts)
	    		f2.write ('%0.8f,%d,%s,%s,%0.2f\n' % (value[2], value[1], value[0], key, 0))

	elif hscale > 0 and hscale < 3:
		for key, value in ranking.iteritems() :
			ranking[key][2] = value[1]/float(global_counts)
			ranking[key][3] = get_hydrophob(str(value[0]))
			f2.write ('%0.8f,%d,%s,%s,%0.2f\n' % (value[2], value[1], value[0], key, value[3]))
	
	f2.write("#, %i, allowed reads, %i, total reads" % (allowedseqs_counts, global_counts))

	f2.close()

	logfile.write("    Finished writing results : " + str(datetime.datetime.now()) + "\n")
	print("    Finished writing results : " + str(datetime.datetime.now()))
	
	logfile.write("    Time to process the file : " + str(datetime.datetime.now() - FileTime) + "\n")
	print("    Time to process the file : " + str(datetime.datetime.now() - FileTime))

	del ranking

	logfile.close()

	## MPI related block ends here

## MPI=>Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)


## MPI related block restarts here to calculate enrichment between pair of files 

## MPI=>Master define chuncks
if COMM.rank == 0:
	# Split into however many cores are available.
	jobs = split(EnrichCompList, COMM.size)
else:
	jobs = None

## MPI=> Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)
results = []

print("the jobs are: %s" % jobs)
# Calculate the enrichement based on two the csv files
for job in jobs:
	# Get file pairs to compare
	enrichfilename = job[0].replace(".csv", "") + "_vs_" + job[1]

	logfile = open("batch_job.log", "a")
	logfile.write("Starting enrichment analysis for %s from %s, task %d : " % (enrichfilename, MPI.Get_processor_name(), COMM.rank+1) + str(datetime.datetime.now()) + "\n")
	print("Starting enrichment analysis for %s from %s, task %d : " % (enrichfilename, MPI.Get_processor_name(), COMM.rank+1) + str(datetime.datetime.now()) + "\n")
	csv_file_1 = get_file_content(pair[0]) # Read 1st .csv file
	csv_file_2 = get_file_content(pair[1]) # Read 2nd .csv file
	file1_fields = {} # allocate dictionary for the 1st file fields
	file2_fields = {} # allocate dictionary for the 2nd file fields

	# split lines in fields for he 1st file
	for line1 in csv_file_1:
		if line1[:9] != "frequency" and line1[:1] != "#" and line1!= "\n":
			line1_fields = line1.split(",")
			file1_fields[line1_fields[3]] = ([line1_fields[0], line1_fields[1], line1_fields[2], line1_fields[4]]) # key = line1_fields[3]; values = [frequency, counts, aa_seq, hydrophobicity] 

	del csv_file_1 # release memory
	 	
	for line2 in csv_file_2:
		if line2[:9] != "frequency" and line2[:1] != "#" and line2!= "\n":
			line2_fields = line2.split(",")
			file2_fields[line2_fields[3]] = ([line2_fields[0], line2_fields[1], line2_fields[2], line2_fields[4]]) # key = line1_fields[3]; values = [frequency, counts, aa_seq, hydrophobicity] 

	del csv_file_2 # release memory

	enrichment = {} # create enrichement dictionary

	# Calculate enrichment for each nucleotide sequence 
	for nucleoseq in file2_fields:
		if nucleoseq in file1_fields.keys():
			enrichment[nucleoseq] = ([float(file2_fields[nucleoseq][0]) / float(file1_fields[nucleoseq][0]), file1_fields[nucleoseq][0], file2_fields[nucleoseq][0], file1_fields[nucleoseq][1], file2_fields[nucleoseq][1], file2_fields[nucleoseq][2], file1_fields[nucleoseq][3]])

	del file1_fields # release memory
	del file2_fields # release memory	
	

	# Rank the sequences based on their enrichment (t[1][0] = first value in the dictionary)
	enrichment=collections.OrderedDict(sorted(enrichment.items(), key=lambda t: t[1][0], reverse=True))

	# Write the enrichement file
	
	enrichfile=open(enrichfilename, "w")
	enrichfile.write("enrichment,freq1,freq2,count1,count2,aa_sequence,nt_sequence,hydrophobicity\n")
	for key, value in enrichment.iteritems() :
		enrichfile.write ('%0.3f,%s,%s,%s,%s,%s,%s,%s\n' % (value[0], value[1], value[2], value[3], value[4], value[5], key, value[6]))
	enrichfile.close()

	logfile.write("Finished writing enrichement file %s from %s, task %d : " % (enrichfilename, MPI.Get_processor_name(), COMM.rank+1) + str(datetime.datetime.now()) + "\n")
	logfile.close()
	print("Finished writing enrichement file %s from %s, task %d : " % (enrichfilename, MPI.Get_processor_name(), COMM.rank+1) + str(datetime.datetime.now()))		

	# Release memory	
	del enrichment

## MPI=>Gather results on rank 0.
results = MPI.COMM_WORLD.gather(results, root=0)

# If all results are finished and gathered, the server writes the final statements in the log file
if COMM.rank == 0:
	# Last log statements
	EndTime = datetime.datetime.now()
	logfile = open("batch_job.log", "a")
	logfile.write("Batch job successfully ended : " + str(EndTime) + "\n")
	print("Batch job successfully ended : " + str(EndTime))

	logfile.write("Duration : " + str(EndTime - StartTime) + "\n")
	print("Duration : " + str(EndTime - StartTime))

	logfile.close()
	quit()
