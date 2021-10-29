#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

#
# Treat .fastq files using the notion of UMIs
# 
# User can specify the parameters that will used for the run in an option file   
# Here an example that can be used as template after deleting # at the begining of each line
#
#@parameters section
#path to the fastq files:../
#reverse complement[0=no,1=yes]:0
#5'anchor sequence:ATTCCGGAAACCGAC
#size of the region of interest:201
#3'anchor sequence:GGATCTTGATAAGAG
#hydrophobicity scale[0 = none, 1 = Kyte-Doolittle, 2 = Miyazawa-Jernigen]:2
#Allowed protein sequences corresponding to the region of interest (use regular expression):[KNTIMQHRLEDAGVYWFCPS*#]{67}
#Allowed UMI at 5' (use regular expression):(^(CGC|CT|A)[ACT][ACT][ACGT][ACT][ACT][ACGT][ACT])
#Allowed UMI at 3' (use regular expression):([AGT][ACGT][AGT][AGT][ACGT][AGT][AGT](GCG|AG|T)$)
#Minimal counts for each merged UMI + Variant to be considered:2
#Minimal ratio between the master variant / second best variant per merged UMI:2
#Minimal number molecules relaible molecules to be considered:2
#Wild type nucleotide sequence:GCAGTTGGTGTCACGGTGGTTCTGATTACCTGCACGTATCACGGCCAGGAATTCATCCGCGTGGGTTATTACGTTAACAATGAATACCTGAACCCGGAACTGCGTGAAAATCCGCCGATGAAACCGGATTTTAGTCAGCTGCAACGCAACATTCTGGCATCCAATCCGCGTGTTACCCGCTTCCATATCAATTGGGATAAC
#Log file name:batch_job.log
#@enrichement calculations section (enrichement = file_2 frequencies / file_1 frequencies) => specify pairs of files for enrichement calculation => file_1 : file_2 
#Sample-01.assembled.fastq:Sample-02.assembled.fastq
#Sample-03.assembled.fastq:Sample-04.assembled.fastq
#
#
#  to get the time of MPI runs use the command: time mpirun -np 4 --mca btl ^openib do_MpiBatchJob_UMIs_v01.py -op options.txt &
#
# Oscar Ramos (20/12/2019)
#

#Import required module(s)
import gzip
import collections
import argparse
import re
import datetime
import locale
import subprocess
import sys
import socket
from os import listdir, path, getcwd
from os.path import isfile, join
from Bio import pairwise2
from Bio.Seq import Seq 

import timeit

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
umi1 = None
umi2 = None
mc = None
mr = None
mm = None
wt_ntseq = None
o = None
EnrichCompList = []

# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='Clean, count, rank and calculate the hydrophobicity score of sequences using a fastq file as input. ', prefix_chars='-+')
parser.add_argument("-imp", help="Import previously calculated UMIsmerge data? (Yes/No) (yes/no)")
parser.add_argument("-p", help="the path to the file that should be analyzed.")
parser.add_argument("-l", help="the anchor sequence at 5' side (left)")
parser.add_argument("-s", help="the size of the region between anchors to extract (number of consecutive bases)")
parser.add_argument("-r", help="the anchor sequence at 3' side (right)")
parser.add_argument("-hs", help="the input hydrophobicity scale to be used. 0: no hydrophobicity calculation; 1: Kyte-Doolittle. J. Mol. Biol. 157:105-132(1982); 2: Miyazawa-Jernigen. Macromolecules 18:534-552(1985)")
parser.add_argument("-rc", help="Use the reverse complement. 0: no; 1: yes.")
parser.add_argument("-allowed", help="Allowed protein sequences corresponding to the region of interest (use regular expression)")
parser.add_argument("-umi1", help="Allowed UMI at 5' (use regular expression)")
parser.add_argument("-umi2", help="Allowed UMI at 3' (use regular expression)")
parser.add_argument("-mc", help="Minimal counts for each allowed sequence to be considered for a merged UMI")
parser.add_argument("-mr", help="Minimal ratio of the best variant counts / second best variant counts per merged UMI")
parser.add_argument("-mm", help="Minimal number of master variant / merged UMI pairs to be retained")
parser.add_argument("-wnt", help="Wild type nucleotide sequence")
parser.add_argument("-o", help="Name of the log file (output).")
parser.add_argument("-op", help="the path to the options file.")

#actually parse command line options
args = parser.parse_args()

if args.op is None: # for command line options
	ImpUMIData = args.imp	
	mypath = args.p
	ANCHOR5 = args.l
	LENGTH = int(args.s)
	ANCHOR3 = args.r
	hscale = int(args.hs)
	rc_boolean = int(args.rc)
	expected = args.allowed
	UMI1_regex = args.umi1
	UMI2_regex = args.umi2
	mincounts = args.mc
	minratio = args.mr
	minmols = args.mm
	wt_ntseq = args.wnt	
	o = args.o

else: # for options file
	optionsfile=open(args.op, "r")
	optionlinenbr = 0
	for line in optionsfile :
		if line[:1] != "@" and line!= "\n":
			if optionlinenbr == 1:
				ImpUMIData = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 2:
				mypath = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 3:
				rc_boolean = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 4:
				ANCHOR5 = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 5:
				LENGTH = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 6:
				ANCHOR3 = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 7:
				hscale = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 8:
				expected = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 9:
				UMI1_regex = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 10:
				UMI2_regex = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 11:
				mincounts = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 12:
				minratio = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 13:
				minmols = int(line.rsplit( ":")[ 1 ])
			if optionlinenbr == 14:
				wt_ntseq = line.rsplit( ":")[ 1 ].replace("\n", "")
			if optionlinenbr == 15:
				o = line.rsplit( ":")[ 1 ].replace("\n", "")
			if optionlinenbr >= 16:
				files = (line.split(":"))
				files[-1]=files[-1].replace("\n", "")
				EnrichCompList.append(files)
		optionlinenbr += 1

#Dictionaries for the messages
rc_boolean_dic = {0:'No',1:'Yes'}
hscale_dic = {0:'None, all values will be set to 0',1:'Kyte-Doolittle',2:'Miyazawa-Jernigen'}

if rc_boolean != 1:
	rc_boolean = 0

if hscale is None or hscale >= 3:
	print (" \n'Warning: No valid hydrophobic scale was choosen! No hydrophobicity calculation will be made and all values will be set to 0.'\n")
	hscale = 0

if COMM.rank == 0:
	# Begin to log job events
	logfile=open(o, "w")
	print("***\nHostname : %s \nUSED OPTIONS" % (socket.gethostname()))
	logfile.write("***\nHostname : %s \nUSED OPTIONS\n" % (socket.gethostname()))
	logfile.write("Import UMIDATA? : " + ImpUMIData + "\n")
	print("Import UMIDATA? : " + ImpUMIData)
	logfile.write("File to be analyzed : " + mypath + "\n")
	print("File to be analyzed : " + mypath)
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
	logfile.write("Allowed 5'UMI pattern : %s\n" % UMI1_regex)
	print("Allowed 5'UMI pattern : %s" % UMI1_regex)
	logfile.write("Allowed 3'UMI pattern : %s\n" % UMI2_regex)
	print("Allowed 3'UMI pattern : %s" % UMI2_regex)
	logfile.write('Mininal number of counts : %i\n' % mincounts)
	print('Mininal number of counts : %i' % mincounts)
	logfile.write('Mininal best/second ratio : %i\n' % minratio)
	print('Mininal best/second ratio : %i' % minratio)
	logfile.write('Minimal number of master variant / merged UMI pairs to be retained : %i\n' % minmols)
	print('Minimal number of master variant / merged UMI pairs to be retained : %i' % minmols)
	logfile.write('Wild type nucleotide sequence : %s\n' % wt_ntseq)
	print('Wild type nucleotide sequence : %s' % wt_ntseq)
	logfile.write('Output log file : %s\n' % o)
	print('Output log file : %s' % o)
	if EnrichCompList:
		for pair in EnrichCompList:
			logfile.write("The enrichment will be calculated for: %s => %s\n" % (pair[0], pair[1]))
			print ("The enrichment will be calculated for: %s => %s" % (pair[0], pair[1]))
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
	try:
#		with open(path.dirname(__file__) + '/' + filename) as f:
		with open(filename) as f:
			data = f.read()
			my_splitlines = data.splitlines()
			f.close()
			return my_splitlines
	except FileNotFoundError:
		print("The file %s was not found !" % (filename))

# transfer the content of a file to the RAM  
def get_gzfile_content(filename):
	try:
		with gzip.open(filename, 'rb') as f:
			data = f.read().decode("utf-8")
			my_splitlines = data.splitlines()
			f.close()
			return my_splitlines
	except FileNotFoundError:
			print("The file %s was not found !" % (filename))

# Find mutations and return a list in which element is composed by position in sequence (pos), residue 1 (res1), residue 2 (res2) 
def find_mutations(seq1, seq2):
	mutation = []
	mutations = []
	try:
		position = None
		res1 = None
		res2 = None
		if len(seq1) == len(seq1):
			for pos in range(0,len(seq1)):
				res1 = seq1[pos]
				res2 = seq2[pos]
				if res1 != res2:
					mutation = (pos, res1, res2)
					mutations.append(mutation)
		del mutation
		if mutations != []:
			return mutations
		else:
			return [["WT"]]
	except ValueError:
        	return "error while trying to find mutations"

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

# Divide list in chunks
def make_chunks(f1_lines, chunks_nbr):
	chunk_size = (int((len(f1_lines) / chunks_nbr)/4) * 4) + 4  
	Mychunks = [f1_lines[x:x+chunk_size] for x in range(0, len(f1_lines), chunk_size)]
	return Mychunks

class CounterContainer:
	def __init__(self):
		self.error_counter = 0
		self.base_counter = 0
		self.deletion_counter = 0
		self.insertion_counter = 0
		self.mutation_counter = 0

#expected_aas = re.compile(r"%s" % expected);


## MPI related block starts here
# MPI=>Master define chuncks
if COMM.rank == 0:
	onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	onlyfiles.sort()
	logfile.close()

	# MPI=> Assign jobs for workers
	jobs = split(onlyfiles, COMM.size)
else:
	jobs = None
	
COMM.barrier()

# MPI=> Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)
cc1_list = []
#print("Hello from %s, task %s, size of Mychunks = %d. Current date and time is : %s \n %s" % (MPI.Get_processor_name(), str(COMM.rank+1), len(Mychunks), str(datetime.datetime.now()), Mychunks [0:4]))

for in_file_name in jobs:
	myseq_flag = 0
	myseqs_list = []
	myallowedseqs_list = []	
	alignments = []
	cc1 = CounterContainer()
	
	logfile = open(o, "a")
#	logfile.write("%s, %s will treated by worker %s : %s \n" % (MPI.Get_processor_name(), in_file_name, str(COMM.rank+1), str(datetime.datetime.now())))
	print("%s will treated by worker %s : %s " % (in_file_name, str(COMM.rank+1), str(datetime.datetime.now())))
	splited_filename = in_file_name.rsplit( ".")
	if splited_filename[len(splited_filename) - 1] != "gz":
		f1_lines = get_file_content(mypath+in_file_name)
		gz = 'FALSE' 
	else:
		f1_lines = get_gzfile_content(mypath+in_file_name)
		gz = 'TRUE'		
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

	del f1_lines # release memory


	Striped_Seqs_List = []

	for Seq in myseqs_list:
		my_ntseq = find_between( Seq , ANCHOR5, ANCHOR3 )
		if len(my_ntseq) > 50:
			Striped_Seqs_List.append(my_ntseq)

	del myseqs_list # release memory

	for CurrSeq in Striped_Seqs_List:
		if CurrSeq != wt_ntseq:
			alignment = pairwise2.align.globalxs(wt_ntseq, CurrSeq, -2, -2)
			errors=find_mutations(alignment[0][0] , alignment[0][1])
			cc1.base_counter += len(CurrSeq)
			cc1.error_counter += len(errors)
			for error in errors:
				if error[1] != "-" and error[2] != "-":
					cc1.mutation_counter += 1
				elif error[1] == "-":
					cc1.insertion_counter += 1
				elif error[2] == "-":
					cc1.deletion_counter += 1
		else:
			cc1.base_counter += len(CurrSeq)

				

	del Striped_Seqs_List
	logfile = open(o, "a")
	logfile.write("\nWorker %s finished analysis of file %s at %s \n" % (str(COMM.rank+1), in_file_name, str(datetime.datetime.now())))
	print("\nWorker %s finished analysis of file %s at %s" % (str(COMM.rank+1), in_file_name, str(datetime.datetime.now())))
	logfile.write("Results for the chunk \nBases = %d\n  Errors = %d, Error_Freq = %0.8f\n  Deletions = %d, Del_Freq = %0.8f\n  Insertions = %d, Ins_Freq = %0.8f\n  Mutations = %d, Mut_Freq = %0.8f\n" % (cc1.base_counter, cc1.error_counter, float(cc1.error_counter / cc1.base_counter), cc1.deletion_counter, float(cc1.deletion_counter / cc1.base_counter), cc1.insertion_counter, float(cc1.insertion_counter / cc1.base_counter), cc1.mutation_counter, float(cc1.mutation_counter / cc1.base_counter)))
	print ("Results for the chunk \nBases = %d\n  Errors = %d, Error_Freq = %0.8f\n  Deletions = %d, Del_Freq = %0.8f\n  Insertions = %d, Ins_Freq = %0.8f\n  Mutations = %d, Mut_Freq = %0.8f" % (cc1.base_counter, cc1.error_counter, float(cc1.error_counter / cc1.base_counter), cc1.deletion_counter, float(cc1.deletion_counter / cc1.base_counter), cc1.insertion_counter, float(cc1.insertion_counter / cc1.base_counter), cc1.mutation_counter, float(cc1.mutation_counter / cc1.base_counter)))
	logfile.close()
	cc1_list.append(cc1)


# MPI=> Wait all workers and gather results on rank 0.

COMM.barrier()
cc1_list =  MPI.COMM_WORLD.gather(cc1_list, root =0)

if COMM.rank == 0:
	cc2_list = []
	for n in range(0, len(cc1_list)):
	    cc2_list = cc2_list + cc1_list[n]
	del cc1_list
	FinalDelCounter = 0
	FinalErrorCounter = 0
	FinalErrorCounter2 = 0
	FinalBaseCounter = 0
	FinalInsCounter = 0
	FinalMutCounter = 0
	for n in range(len(cc2_list)):
		FinalBaseCounter += cc2_list[n].base_counter
		FinalErrorCounter += cc2_list[n].error_counter
		FinalMutCounter += cc2_list[n].mutation_counter
		FinalDelCounter += cc2_list[n].deletion_counter
		FinalInsCounter += cc2_list[n].insertion_counter

	logfile = open(o, "a")
	print ("Final Results\nBases = %d\n  Errors = %d, Error_Freq = %0.8f\n  Deletions = %d, Del_Freq = %0.8f\n  Insertions = %d, Ins_Freq = %0.8f\n  Mutations = %d, Mut_Freq = %0.8f\n" % (FinalBaseCounter, FinalErrorCounter, float(FinalErrorCounter / FinalBaseCounter), FinalDelCounter, float(FinalDelCounter / FinalBaseCounter), FinalInsCounter, float(FinalInsCounter / FinalBaseCounter), FinalMutCounter, float(FinalMutCounter / FinalBaseCounter)))
	logfile.write("\nFinal Results\n  Bases = %d\n  Errors = %d, Error_Freq = %0.8f\n  Deletions = %d, Del_Freq = %0.8f\n  Insertions = %d, Ins_Freq = %0.8f\n  Mutations = %d, Mut_Freq = %0.8f\n" % (FinalBaseCounter, FinalErrorCounter, float(FinalErrorCounter / FinalBaseCounter), FinalDelCounter, float(FinalDelCounter / FinalBaseCounter), FinalInsCounter, float(FinalInsCounter / FinalBaseCounter), FinalMutCounter, float(FinalMutCounter / FinalBaseCounter), ))
	logfile.write("\nFinished : %s \n" % (str(datetime.datetime.now())))
	print("Finished : %s \n" % (str(datetime.datetime.now())))
	logfile.close()
else:
	pass


quit()

