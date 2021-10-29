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
parser.add_argument("-p", help="the path to the directory containing the fastq files that should be analyzed.")
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

# Begin to log job events
logfile=open(o, "w")
print("***\nHostname : %s \nUSED OPTIONS" % (socket.gethostname()))
logfile.write("***\nHostname : %s \nUSED OPTIONS\n" % (socket.gethostname()))
logfile.write("Import UMIDATA? : " + ImpUMIData + "\n")
print("Import UMIDATA? : " + ImpUMIData)
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
		with open(filename) as f:
			data = f.read()
			my_splitlines = data.splitlines()
			f.close()
			return my_splitlines
	except FileNotFoundError:
		print("The file %s was not found !" % (filename))

# transfer the content of a file to the RAM  
def get_gzfile_content(filename):
	with gzip.open(filename, 'rb') as f:
		data = f.read().decode("utf-8")
		my_splitlines = data.splitlines()
		f.close()
		return my_splitlines

# Fill a set container with NR_UMIs and return it
def get_NRUMIs_set(f1_lines):
	NRUMIs = set()
	myseq_flag = 0
	seqnbr = 0
	for line in f1_lines:
		if myseq_flag == 1:
			if rc_boolean == 1:
				my_seq = rc_dna(line)
			else:
				my_seq = line
			if expected_UMI01.search(line):
				UMI5p = expected_UMI01.findall(line)[0][0]
				if expected_UMI02.search(line):
					UMI3p = expected_UMI02.findall(line)[0][0]
					my_ntseq = find_between( my_seq , ANCHOR5, ANCHOR3 )
					if len(my_ntseq) == LENGTH:
						my_protseq = translate_dna(my_ntseq)
						if expected_aas.search(my_protseq):
							UMIsmerge = (("%s|%s") % (UMI5p,UMI3p))
							NRUMIs.add(UMIsmerge)
							del UMIsmerge
			myseq_flag = 0
			seqnbr += 1
#			sys.stdout.write("\r      Getting NR merged UMIs. Now at sequence %d of %d : %d NR merged UMIs found!         \r" % (seqnbr, (len(f1_lines)/4), len(NRUMIs)))
#			sys.stdout.flush()
		if myseq_flag == 0 and line[:1] == "@":
			myseq_flag = 1
	return NRUMIs
	del NRUMIs

# Construct a set for NR NTs from a dictionary
def get_NRNTs_set(UMIDATA):
	NRNTs = set()
	for Umerge in UMIDATA_0:
		for varnbr in range(0, len(UMIDATA_0[Umerge])):
			if UMIDATA_0[Umerge][varnbr][4] != 0:
				myfactor = UMIDATA_0[Umerge][varnbr][4]
			else:
				myfactor = 0.1
			if UMIDATA_0[Umerge][varnbr][5] == "Yes" and UMIDATA_0[Umerge][varnbr][2] >= mincounts and (UMIDATA_0[Umerge][varnbr][3] / myfactor) >= minratio:
				NRNTs.add(UMIDATA_0[Umerge][varnbr][0])
	return NRNTs
	del NRNTs

# For each NRUMI get the associated variants and return the corresponding dictionary 
def get_Var4NRUMIs(f1_lines, NRUMIs):
	UMIDATA = {}
	UMIDATA = dict.fromkeys(NRUMIs)
	seqnbr = 0
	myseq_flag = 0
	for line in f1_lines:
		if myseq_flag == 1:
			if rc_boolean == 1:
				my_seq = rc_dna(line)
			else:
				my_seq = line
			if expected_UMI01.search(line):
				UMI5p = expected_UMI01.findall(line)[0][0]
				if expected_UMI02.search(line):
					UMI3p = expected_UMI02.findall(line)[0][0]
					my_ntseq = find_between( my_seq , ANCHOR5, ANCHOR3 )
					if len(my_ntseq) == LENGTH:
						my_protseq = translate_dna(my_ntseq)
						if expected_aas.search(my_protseq):
							Umerge = (("%s|%s") % (UMI5p,UMI3p))
							seqnbr += 1
							if UMIDATA[Umerge] != None:
								found_flag = 0
								for varnbr in range(0, len(UMIDATA[Umerge])):
									if my_ntseq == UMIDATA[Umerge][varnbr][0]:
										UMIDATA[Umerge][varnbr][2] += 1
										found_flag = 1
										break
								if found_flag == 0:
									UMIDATA[Umerge].append([my_ntseq,my_protseq,1,0,0,"No"]) # my_ntseq (0), my_protseq (1), counts(2), mastercounts(3), secondcounts(4), ismaster?(5)	
							else:
								UMIDATA[Umerge] = [[my_ntseq,my_protseq,1,0,0,"No"]] # my_ntseq (0), my_protseq (1), counts(2), mastercounts(3), secondcounts(4), ismaster?(5)
							del Umerge
			myseq_flag = 0
#			sys.stdout.write("\r      Getting variants for NR UMIs. Now at accepted sequence %d of %d max.         \r" % (seqnbr, (len(f1_lines)/4)))
#			sys.stdout.flush()
		if myseq_flag == 0 and line[:1] == "@":
			myseq_flag = 1
		del line
	return UMIDATA
	del UMIDATA

def find_master(in_file_name, UMIDATA):
	# Find master and second variant and write results to file
	UMIsmergeFN = in_file_name.rsplit(".")[ 0 ] + "_UMIsmerge.csv.gz"
	print("      Now will going to write the csv file : %s " % (UMIsmergeFN))
	logfile.write("      Now will going to write the csv file : %s \n" % (UMIsmergeFN))
	UMIcouplesFile = gzip.open(UMIsmergeFN, 'wt')
	UMIcouplesFile.write("merged UMIs(5'|3'),nt_sequence,aa_sequence,counts,UMI_counts,master_counts, second_counts,is_master?\n")			
	NowAtSeq = 0
	UMIcounts = 0
	for Umerge in UMIDATA:
		best_count = 0
		second_best_count = 0
		equals_best = 0
		# Find the counts value of the master and the second most prevalent variant for the UMI
		for varnbr in range(0, len(UMIDATA[Umerge])):
			UMIcounts += UMIDATA[Umerge][varnbr][2]
			if UMIDATA[Umerge][varnbr][2] > best_count:
				second_best_count = best_count
				best_count = UMIDATA[Umerge][varnbr][2]
			elif UMIDATA[Umerge][varnbr][2] == best_count:
			 	second_best_count = best_count
			elif UMIDATA[Umerge][varnbr][2] < best_count and UMIDATA[Umerge][varnbr][2] > second_best_count:
			 	second_best_count = UMIDATA[Umerge][varnbr][2]
			else:
				pass
		if best_count > second_best_count:
			hasmaster = "Yes"
		else:
			hasmaster = "No"
		# write results
		for varnbr in range(0, len(UMIDATA[Umerge])):
			NowAtSeq += 1
			UMIcouplesFile.write("%s,%s,%s,%d,%d,%d,%d," % (Umerge, UMIDATA[Umerge][varnbr][0], UMIDATA[Umerge][varnbr][1], UMIDATA[Umerge][varnbr][2], UMIcounts, best_count, second_best_count))
			if hasmaster == "Yes" and UMIDATA[Umerge][varnbr][2] == best_count:
				UMIDATA[Umerge][varnbr][3] = best_count
				UMIDATA[Umerge][varnbr][4] = second_best_count
				UMIDATA[Umerge][varnbr][5] = "Yes"
				UMIcouplesFile.write("%s\n" % ("Yes"))
			else:
				UMIcouplesFile.write("%s\n" % ("No"))
#			sys.stdout.write("      Now at the retained NR merged UMI + variant %d. %d sequences with the merged UMI %s comprising %d variants      \r" % (NowAtSeq, UMIcounts, Umerge, len(UMIDATA[Umerge])))
#			sys.stdout.flush()
		UMIcounts = 0	
	UMIcouplesFile.close()
	return UMIDATA
	del UMIDATA

# Import previously collected NR merged UMIs data
def ImpUMIDATA(UMIsmergeFN):
	UMISmerge_lines = get_file_content(UMIsmergeFN)
	del UMISmerge_lines[0]
	NRUMIs = set()
	for line in UMISmerge_lines:
		splited_line = line.rsplit(",")
		NRUMIs.add(splited_line[0])
	UMIDATA = dict.fromkeys(NRUMIs)
	for line in UMISmerge_lines:
		splited_line = line.rsplit(",")
		if UMIDATA[splited_line[0]] != None:
			UMIDATA[splited_line[0]].append([splited_line[1], splited_line[2], int(splited_line[3]), int(splited_line[5]), int(splited_line[6]), splited_line[7]])
		else:
			UMIDATA[splited_line[0]] = [[splited_line[1], splited_line[2], int(splited_line[3]), int(splited_line[5]), int(splited_line[6]), splited_line[7]]]
	return UMIDATA

# Get masters that correspond to mincounts and minratio criteria
def get_good(UMIDATA, NTDATA):
	GlobalNTcounts = 0
	for Umerge in UMIDATA:
		for varnbr in range(0, len(UMIDATA[Umerge])):
			if UMIDATA[Umerge][varnbr][4] != 0:
				myfactor = UMIDATA[Umerge][varnbr][4]
			else:
				myfactor = 0.1
			if UMIDATA[Umerge][varnbr][5] == "Yes" and UMIDATA[Umerge][varnbr][2] >= mincounts and (UMIDATA[Umerge][varnbr][3] / float(myfactor)) >= minratio:
				if NTDATA[UMIDATA[Umerge][varnbr][0]] == None:
					NTDATA[UMIDATA[Umerge][varnbr][0]] = [UMIDATA[Umerge][varnbr][1], 1 , (UMIDATA[Umerge][varnbr][2] / float(myfactor))] # my_ntseq (Key), my_protseq (value 0), counts(value 1), ratio(value 2)
					GlobalNTcounts += 1
				else:
					NTDATA[UMIDATA[Umerge][varnbr][0]][1] += 1
					GlobalNTcounts += 1	
	return NTDATA, GlobalNTcounts
	del NTDATA, GlobalNTcounts

# Write every reliable NT sequence to a file
def write_NRNTs2file(in_file_name, NTDATA, FilteredNTDATA):
	NTf = in_file_name.rsplit(".")[ 0 ] + "_NT.csv"
	print("      Now will going to write the csv file : %s " % (NTf))
	logfile.write("      Now will going to write the csv file : %s \n" % (NTf))
	NTfile = open(NTf, "w")
	NTfile.write("nt_sequence,aa_sequence,counts,ratio,frequency\n")
	FilteredNTDATA[in_file_name] = []			
	for NT in NTDATA:
		if NTDATA[NT][1] >= minmols:
			freq = NTDATA[NT][1] / float(GlobalNTcounts)
			NTfile.write("%s,%s,%d,%f,%0.8f\n" % (NT, NTDATA[NT][0], NTDATA[NT][1], NTDATA[NT][2], freq ))
			FilteredNTDATA[in_file_name].append([NT, NTDATA[NT][0], NTDATA[NT][1], NTDATA[NT][2], freq]) # nt seq [0], prot_seq [1], count of molecules[2], ratio best/second[3], frequency[4]
	NTfile.write("Total molecule counts,%d,,\n" % (GlobalNTcounts))
	NTfile.close()

	print("      Found %d reliable starting molecules !" % (len(FilteredNTDATA[in_file_name])))
	logfile.write("      Found %d reliable starting variants !\n" % (len(FilteredNTDATA[in_file_name])))

	print("    Ended : %s" % (str(datetime.datetime.now())))
	logfile.write("    Ended : %s\n" % (str(datetime.datetime.now())))
	return FilteredNTDATA
	del FilteredNTDATA
	
# Calculate enrichment and write results to file
def enrich_calc(jobs,EnrichCompList):
	for entry in range(0, len(jobs)):
		Enrichf = jobs[entry][0].rsplit(".")[ 0 ] + "_vs_" + jobs[entry][1].rsplit(".")[ 0 ] + "_Enrichment.csv"
		MyEnrichfile = open(Enrichf, "w")
		MyEnrichfile.write("nt_sequence,aa_sequence,counts_before,counts_after,frequency_before,frequency_after,enrichment,nt_mutations,aa_mutations\n")
		print("Now will going to write the csv file : %s " % (Enrichf))
		logfile.write("Now will going to write the csv file : %s \n" % (Enrichf))
		File_0 = EnrichCompList[entry][0].replace("_NT.csv", ".fastq")
		File_1 = EnrichCompList[entry][1].replace("_NT.csv", ".fastq")
		for varnbr_0 in range(0, len(GatheredData[File_0])):
			for varnbr_1 in range(0, len(GatheredData[File_1])):
				if GatheredData[File_0][varnbr_0][ 0 ]  == GatheredData[File_1][varnbr_1][ 0 ]:
					enrichment = float(GatheredData[File_1][varnbr_1][ 4 ]) / float(GatheredData[File_0][varnbr_0][ 4 ])
					# FinalDATA[NT] = [aa_sequence[0],counts_before[1],counts_after[2],frequency_before[3],frequency_after[4],enrichment[5],nt_mutations[6],aa_mutations[7]] 
					FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]] = [ GatheredData[File_0][varnbr_0][1], GatheredData[File_0][varnbr_0][2], GatheredData[File_1][varnbr_1][2], GatheredData[File_0][varnbr_0][4], GatheredData[File_1][varnbr_1][4], enrichment, [], [] ]
					wt_aaseq = translate_dna(wt_ntseq)
					FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6] = find_mutations(wt_ntseq, GatheredData[File_0][varnbr_0][ 0 ]) # find the mutations (position=mutation[0], residue 1=mutation[1] and residue 2=mutation[2])
					FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7] = find_mutations(wt_aaseq, GatheredData[File_0][varnbr_0][ 1 ]) # find the mutations (position=mutation[0], residue 1=mutation[1] and residue 2=mutation[2])
					MyEnrichfile.write("%s,%s,%d,%d,%0.8f,%0.8f,%0.8f," % (GatheredData[File_0][varnbr_0][ 0 ], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][0], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][1], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][2], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][3], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][4], FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][5]))
					for mutanbr in range(0, len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6])):
						if mutanbr < (len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6])-1) and len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6]) != 1:
							if FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][0]  != "WT":
								muta = [FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][1], str(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][0]+1), FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][2]]
								MyEnrichfile.write(''.join(muta))
								MyEnrichfile.write("|")
							else:
								MyEnrichfile.write("WT")
						else:
							if FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][0] != "WT":
								muta = [FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][1], str(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][0]+1), FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6][mutanbr][2]]
								MyEnrichfile.write(''.join(muta))
							else:
								MyEnrichfile.write("WT")
					MyEnrichfile.write(",")
					for mutanbr in range(0, len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7])):
						if mutanbr < (len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7])-1) and len(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][6]) != 1:
							if FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][0] != "WT":
								muta = [FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][1], str(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][0]+1), FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][2]]
								MyEnrichfile.write(''.join(muta))
								MyEnrichfile.write("|")
							else:
								MyEnrichfile.write("WT")
						else:
							if FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][0] != "WT":
								muta = [FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][1], str(FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][0]+1), FinalDATA[GatheredData[File_0][varnbr_0][ 0 ]][7][mutanbr][2]]
								MyEnrichfile.write(''.join(muta))
							else:
								MyEnrichfile.write("WT")
					MyEnrichfile.write("\n")
		MyEnrichfile.close()

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

# Do the job. Get fastq file list and for each fastq file, get sequences, isolate the region of interest between 5' and 3' anchors, 
# translate, clean non-expected sequences and count their occurences. If a hydrophobicity scale is choosen then a score is calculated
# for each sequence using it.
# The output is a fasta file corresponding to the isolated region after translation, an weblogo file and the ranked sequences in csv format 

expected_aas = re.compile(r"%s" % expected);
expected_UMI01 = re.compile(r"%s" % UMI1_regex);
expected_UMI02 = re.compile(r"%s" % UMI2_regex);
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

FilteredNTDATA = {}
FinalDATA = {}
GatheredData = {}

## MPI related block starts here
# MPI=>Master define chuncks
if COMM.rank == 0:
	# Split into however many cores are available.
	jobs = split(onlyfiles, COMM.size)
else:
	jobs = None

# MPI=> Scatter jobs across cores.
jobs = COMM.scatter(jobs, root=0)

# Log file update
logfile = open(o, "a")

for in_file_name in jobs:
	if ImpUMIData != "Yes" and ImpUMIData != "yes":
#		i = 0
#		my_regions_list = []
#		myallowedseqs_list = []
		FileTime = datetime.datetime.now()		
		# Log file update
		logfile.write("Reading FASTQ file " + in_file_name + " from " + MPI.Get_processor_name() + ", task " + str(COMM.rank+1) + ". Current date and time is : " + str(datetime.datetime.now()) + "\n")
		print("Reading FASTQ file " + in_file_name + " from " + MPI.Get_processor_name() + ", task " + str(COMM.rank+1) + ". Current date and time is : " + str(datetime.datetime.now()))
		print("    Phase I: Identification of sequence regions, translation, filtering and creation of a non-reduntant (NR) set of merged UMIs.")
		logfile.write("    Phase I: Identification of sequence regions, translation, filtering and creation of a non-reduntant (NR) set of merged UMIs.\n")
		splited_filename = in_file_name.rsplit( ".")
		if splited_filename[len(splited_filename) - 1] != "gz":
			f1_lines = get_file_content(mypath+in_file_name)
			gz = 'FALSE' 
		else:
			f1_lines = get_gzfile_content(mypath+in_file_name)
			gz = 'TRUE'		

		# Find NR merged UMIs put them in a set
		NRUMIs = set()
		NRUMIs=get_NRUMIs_set(f1_lines)

		# Log file update
		logfile.write("      Finished reading FASTQ file. Found %d NR merged UMIs from %d reads : %s \n" % (len(NRUMIs), len(f1_lines)/4, str(datetime.datetime.now())))
		print("      Finished reading FASTQ file. Found %d NR merged UMIs from %d reads : %s " % (len(NRUMIs), len(f1_lines)/4, str(datetime.datetime.now())))

		#Find and count the number of variants for each merged UMI
		UMIDATA_0 = get_Var4NRUMIs(f1_lines, NRUMIs)
		
		# release memory
		del f1_lines 
		del NRUMIs

		# Log file update
		logfile.write("      Finished parsing FASTQ sequences : " + str(datetime.datetime.now()) + "\n")
		print("      Finished parsing FASTQ sequences : " + str(datetime.datetime.now()))
		print("    Phase II: Now will find the master variant for each UMI")
		logfile.write("    Phase II: Now will find the master variant for each UMI\n")

		# Find master and the second best
		UMIDATA_0 = find_master(in_file_name, UMIDATA_0)

	else:
		UMIsmergeFN = in_file_name.rsplit(".")[ 0 ] + "_UMIsmerge.csv"
		# Log file update
		logfile.write("Importing log file %s : %s\n" % (UMIsmergeFN, str(datetime.datetime.now())))
		print("Importing file %s : %s" % (UMIsmergeFN, str(datetime.datetime.now())))
		UMIDATA_0 = ImpUMIDATA(UMIsmergeFN)
		print("    Done importing file %s : %s." % (UMIsmergeFN, str(datetime.datetime.now())))
		logfile.write("    Done importing file %s : %s.\n" % (UMIsmergeFN, str(datetime.datetime.now())))
		


        #Get reliable molecules
	
	# Build a set of NRNTs and create a new dictionar from it
	NRNTs=get_NRNTs_set(UMIDATA_0)
	NTDATA = dict.fromkeys(NRNTs)
	del NRNTs
	
	# Get masters that correspond to mincounts and minratio criteria
	NTDATA,GlobalNTcounts = get_good(UMIDATA_0, NTDATA)
	del UMIDATA_0
	# Log file update
	print("      Found %d molecules !" % (GlobalNTcounts) )
	logfile.write("      Found %d molecules !\n" % (GlobalNTcounts) )
	
	# Write every reliable NT sequence to a file
	allNTDATA=write_NRNTs2file(in_file_name, NTDATA, FilteredNTDATA)

	del NTDATA


COMM.barrier()
logfile.close()
logfile = open(o, "a")


## MPI=>Gather results on rank 0.
FilteredNTDATA = MPI.COMM_WORLD.gather(FilteredNTDATA, root=0)

if COMM.rank == 0:
	for itemnbr in range(0,len(FilteredNTDATA)):
		for MyF in FilteredNTDATA[itemnbr]:
			GatheredData[MyF] = FilteredNTDATA[itemnbr][MyF]
	del FilteredNTDATA

COMM.barrier()

# MPI=>Master define chuncks
if COMM.rank == 0:
	# Split into however many cores are available.
	jobs = split(EnrichCompList, COMM.size)
else:
	jobs = None

# MPI=> Scatter enrichment jobs across cores.
jobs = COMM.scatter(jobs, root=0)
GatheredData = COMM.bcast(GatheredData, root=0)

# Do enrichment calculations
enrich_calc(jobs,EnrichCompList)

COMM.barrier()
logfile.close()
logfile = open(o, "a")

if COMM.rank == 0:
	print("Ended : %s" % (str(datetime.datetime.now())))
	logfile.write("Ended : %s\n" % (str(datetime.datetime.now())))	

	logfile.close()
	quit()
