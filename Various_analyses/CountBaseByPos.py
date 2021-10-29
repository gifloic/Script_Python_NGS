#!/usr/bin/env python

# This code is intended to point mutations in non-redundant sequence list compared to an WT sequence

# User can specify the parameters that will used for the run in an option file   
# Here an example that can be used as template after deleting # at the begining of each line
#
#@options section
#path to the fastq files:../
#
# Oscar Ramos (mer. 24 oct. 2018 07:47:50 CEST )

import argparse
import datetime
import locale
from itertools import combinations

#import pandas as pd
#import matplotlib.pyplot as plt
#import numpy as np

# American global settings
locale.setlocale(locale.LC_ALL, 'C')

# Global variables section


# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='This code is intended to evalue the fitness associated to a residue in a given position. ', prefix_chars='-+')
parser.add_argument("-I", help="the path to the non-redundant file to be analyzed")
parser.add_argument("-wt", help="the sequence corresponding to wt")
parser.add_argument("-o", help="the name of the output file to be generated")
parser.add_argument("-op", help="the path to the options file.")

# Actually parse command line options
args = parser.parse_args()

if args.op is None: # for command line options
	SRC = args.I
if args.op is None: # for command line options
	WT = args.wt
if args.op is None: # for command line options
	out_file_name = args.o

else: # for options file
	optionsfile=logfile=open(args.op, "r")
	optionlinenbr = 0
	for line in optionsfile :
		if line[:1] != "@" and line!= "\n":
			if optionlinenbr == 1:
				SRC = (line.rsplit( ":")[ 1 ]).replace("\n", "")
				print(line)
			if optionlinenbr == 2:
				WT = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 3:
				out_file_name = (line.rsplit( ":")[ 1 ]).replace("\n", "")
				print(out_file_name)
		optionlinenbr += 1


# Transfer the content of a file to RAM  
def get_file_content(filename):
    with open(filename) as f:
	data = f.read()
	my_splited_splitines = []
	my_splitlines = data.splitlines()
	f.close()
	for i in range(len(my_splitlines)):
		if len(my_splitlines[0].rsplit(",")) <= len(my_splitlines[i].rsplit(",")): 	
			my_splited_splitines.append(my_splitlines[i].rsplit(","))
	return my_splited_splitines

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
					mutation = (pos+1, res1, res2)
					mutations.append(mutation)
		del mutation
		return mutations
	except ValueError:
        	return "error while trying to find mutations"

# Begin to log job events
logfile=open("job.log", "w")
logfile.write("***\nUSED OPTIONS\n")
print("***\nUSED OPTIONS")
logfile.write("File containing mutated sequences : " + SRC + "\n")
print("File containing mutated sequences : " + SRC  )
StartTime = datetime.datetime.now()
logfile.write("Job started : " + str(StartTime) + "\n")
print("Job started : " + str(StartTime))

# Read and parse data file. Load it to RAM
data_lines = get_file_content(SRC)
PosBaseDict = {}
n = 0

# Analysis
for line in data_lines:  										# for each read
	mutations=find_mutations(WT, line[0])  								# find the mutations to the reference (WT)
	for mutation in mutations:
		res2_flag = 0  									# for each mutated position
		if mutation[0] in PosBaseDict.keys():							# if the mutated position is already record in PosBaseDict
			for recordedmutnbr in range(0, len(PosBaseDict[mutation[0]])):			# for each mutatation recorded for this position
				if mutation[2] == PosBaseDict[mutation[0]][recordedmutnbr][0]: 		# verify if the mutated residued is already recorded
					count = PosBaseDict[mutation[0]][recordedmutnbr][1] + 1 			# if it is the case sum up 1 to the mutated base for this position
					PosBaseDict[mutation[0]][recordedmutnbr] = [PosBaseDict[mutation[0]][recordedmutnbr][0], count]
					res2_flag = 1
#					print("Pos and res found: adding to count")
			if res2_flag == 0:
				PosBaseDict[mutation[0]].append([mutation[2], 1])
				res2_flag = 1
#				print("Pos found and including new res")	
		else:											# if it is not the case record the mutated base for this position and assign the value = 1 for its count 
			PosBaseDict[mutation[0]] = []
			PosBaseDict[mutation[0]].append([mutation[2], 1])
			res2_flag = 1	
#			print("Pos not found: Including pos and res")
			if len(mutations) == 1:
				break		
	n += 1

print ("Finished counting mutated bases per position for %d reads: %s"  % ( n , str(datetime.datetime.now())))
logfile.write("Finished counting mutated bases per position for %d reads: %s\n"  % ( n , str(datetime.datetime.now())))

myoutfile=open(out_file_name, "wt")
acceptedreslist = ["A","C","G","T","N","-"]
myoutfile.write('nt')
for i in range(1, len(data_lines[0][0])+1):					# print the position numbers
	myoutfile.write(',%d' % (i))
myoutfile.write('\n')
for n in range(0, len(acceptedreslist)):					# for each accepted residue
	myoutfile.write('%s' % (acceptedreslist[n]))				# write the name of the current accepted residue to the output
	for x in range(1, len(data_lines[0][0])+1):     			# for each position in the sequence alignment
		if x in PosBaseDict.keys():           				# if the position corresponds was recorded in possibleres dict
			for w in range(0, len(PosBaseDict[x])):			# for each mutation recorded for the position
				if acceptedreslist[n] == PosBaseDict[x][w][0]: 	# if current accepted residue = current mutation for the position
					myoutfile.write(',%d' % (PosBaseDict[x][w][1]))
					
					break
			else:
				myoutfile.write(',%d' % (0))
		else:
			myoutfile.write(',%d' % (0))
	myoutfile.write('\n')
myoutfile.close()

print ("\tDone wrinting results to file %s: %s"  % (out_file_name, str(datetime.datetime.now())))
logfile.write("\tDone wrinting results to file %s: %s\n"  % (out_file_name, str(datetime.datetime.now())))

EndTime = datetime.datetime.now()
print("Batch job successfully ended : " + str(EndTime))
logfile.write("Batch job successfully ended : " + str(EndTime) + "\n")
print("Duration : " + str(EndTime - StartTime))
logfile.write("Duration : " + str(EndTime - StartTime) + "\n")
logfile.close()


#N = len(data_lines[0][0])
#A = (20, 35, 30, 35, 27)
#womenMeans = (25, 32, 34, 20, 25)
#menStd = (2, 3, 4, 1, 2)
#womenStd = (3, 5, 2, 3, 3)
#ind = np.arange(N)    # the x locations for the groups
#width = 0.35       # the width of the bars: can also be len(x) sequence

#p1 = plt.bar(ind, Acounts, width, yerr=menStd)
#p2 = plt.bar(ind, womenMeans, width,
#             bottom=menMeans, yerr=womenStd)

#plt.ylabel('Scores')
#plt.title('Scores by group and gender')
#plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
#plt.yticks(np.arange(0, 81, 10))
#plt.legend((p1[0], p2[0]), ('Men', 'Women'))

#plt.show()
