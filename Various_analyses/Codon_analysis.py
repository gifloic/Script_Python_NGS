#!/usr/bin/env python

#Options file example:
#@options section
#path to the input file:../VN1053-avant-Kana.assembled_vs_VN1053-apres-Kana.assembled.csv
#path to the output file:Out_results.csv
#path to the output file without WT codons:Out_results_noWTcodons.csv
#WT sequence:GATGACCTGGAATGGAAAATTATCTATGTGGGCAGCGCGGAATCTGAAGAATTTGATCAGATCCTGGACAGCGTCCTGGTGGGCCCGGTCCCGGCAGGTCGTCATATGTTTGTTTTCCAAGCAGATGCTCCGAACCCGTCGCTGATTCCGGAAACCGACGCAGTTGGTGTCACGGTGGTTCTGATTACCTGC
#Factor for residue position correction:35
#
# Oscar Ramos (13/02/2019 16:32:36)

import argparse
import collections
import datetime
import locale
import sys
from os import listdir
from os.path import isfile, join

# American global settings
locale.setlocale(locale.LC_ALL, 'C')

# Global variables section


# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='This code is intended to count codons. ', prefix_chars='-+')
parser.add_argument("-op", help="the path to the options file.")


# Actually parse command line options
args = parser.parse_args()

# Request and parse options file in order define required variables 
if args.op is None:
	print("Please provide an options file by using the -op option")
else:
	optionsfile=logfile=open(args.op, "r")
	optionlinenbr = 0
	for line in optionsfile :
		if line[:1] != "@" and line!= "\n":
			if optionlinenbr == 1:
				InFile = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 2:
				OutFile = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 3:
				OutFilenoWT = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 4:
				WTseq = (line.rsplit( ":")[ 1 ]).replace("\n", "")
			if optionlinenbr == 5:
				background = float((line.rsplit( ":")[ 1 ]).replace("\n", ""))
		optionlinenbr += 1

# Transfer the content of a file to RAM  
def get_file_content(filename):
	try:
		with open(filename) as f:
			data = f.read()
			my_parsed_splitines = []
			my_splitlines = data.splitlines()
			f.close()
			for i in range(len(my_splitlines)):
				if len(my_splitlines[0].rsplit(",")) <= len(my_splitlines[i].rsplit(",")): 	
					my_parsed_splitines.append(my_splitlines[i].rsplit(","))
		return my_parsed_splitines

	except ValueError:
        	return "error while trying to treat file "  + filename
	

# Take two sequences, find mutations and return a list in which each element (mutation) is composed by position in sequence (pos), residue 1 (res1), residue 2 (res2) 
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
		return mutations

	except ValueError:
        	return "error while trying to find mutations"

# Evaluate if a chain of characters corresponds to a number or not
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
    return False

# translate a DNA sequence to protein sequence
def translate_dna(sequence):
	try:
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
		return proteinsequence

	except ValueError:
        	return "error while trying to translate DNA sequence: " + sequence

# Get codon from a DNA sequence (the size of the sequence should be a multiple of 3)
def get_codons(sequence):
	try:
		pos_vs_codon = []
		for x in range(0,len(sequence),3):
			end = x+3
			codonnbr = (end/3)
			codonseq = variant[x:end]
			pos_vs_codon.append([codonnbr,codonseq])
		return pos_vs_codon
	
	except ValueError:
        	return "error while trying to get codons from DNA sequence: " + sequence

				
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
logfile=open("job_sc.log", "w")
logfile.write("***\nUSED OPTIONS\n")
print("***\nUSED OPTIONS")
logfile.write("Directory containing the files to be analyzed : " + InFile + "\n")
print("Directory containing the files to be analyzed : " + InFile)
StartTime = datetime.datetime.now()
logfile.write("Job started : " + str(StartTime) + "\n")
print("Job started : " + str(StartTime))

# Create an empty dictionnary
VariantDict = {}
CodonPosList = {}


# Read, parse and analyze data files (in this case the outuput from an enrichment analysis) 
data_lines = get_file_content(InFile)

# Identify mutations and list them to a output file 
myoutfile=open("Mutations.csv", "wt")
for linenbr in range(len(data_lines)):
	if linenbr != 0:
		mutations=find_mutations(WTseq, data_lines[linenbr][6])
	#	print (line[0] + "," + str(len(mutations)))
	#	if len(mutations) == 3:
	#		if (mutations[2][0] - mutations[0][0]) < 3:
		myoutfile.write('%s,%d,' % (data_lines[linenbr][6], len(mutations)))
		for mutation in mutations:
			myoutfile.write('%s%d%s,' % (mutation[1],mutation[0],mutation[2]))
		myoutfile.write("\n")
myoutfile.close
	
# Assume DNA variants as keys to the VariantDict and their frequencies as values
for linenbr in range (1, len(data_lines)):
	VariantDict[data_lines[linenbr][6]] = (float(data_lines[linenbr][1]), float(data_lines[linenbr][2])) # get only each variant sequence, freq1 and freq2

# Calculate the cumulated frequencies for all available positions and codons
for variant in VariantDict.keys():
	MySeqCodons=get_codons(variant)
	for x in range(0,len(MySeqCodons)):
		codonflag = 0
		if MySeqCodons[x][0] in CodonPosList.keys():
			for y in range(0, len(CodonPosList[MySeqCodons[x][0]])):
				if CodonPosList[MySeqCodons[x][0]][y][0] == MySeqCodons[x][1]:
#					print("Codon_%i %s found" % (MySeqCodons[x][0], MySeqCodons[x][1]))
					updateddata = [MySeqCodons[x][1], (CodonPosList[MySeqCodons[x][0]][y][1] + VariantDict[variant][0]), (CodonPosList[MySeqCodons[x][0]][y][2]+ VariantDict[variant][1])]
					CodonPosList[MySeqCodons[x][0]][y] = updateddata
#					print("Variant=%s\ncodonnbr_%i = %s" % (variant, MySeqCodons[x][0],str(CodonPosList[MySeqCodons[x][0]])))
					codonflag = 1
			if codonflag == 0:
#				print("Codon_%i %s not found" % (MySeqCodons[x][0], MySeqCodons[x][1]))
				CodonPosList[MySeqCodons[x][0]].append([MySeqCodons[x][1], float(VariantDict[variant][0]), float(VariantDict[variant][1])])
		else:
			CodonPosList[MySeqCodons[x][0]] = []
			CodonPosList[MySeqCodons[x][0]].append([MySeqCodons[x][1], float(VariantDict[variant][0]),  float(VariantDict[variant][1])])   # Key = codon position; Value 0 = codon triplet; Value 1 = freq1; Value 2 = freq2  
#			print("Variant=%s\ncodonnbr_%i = %s" % (variant, MySeqCodons[x][0],str(CodonPosList[MySeqCodons[x][0]])))
	#		print(variant[x:end])
	#		print(MySeqCodons[x][0])


#Add enrichment values to the dictinary CodonPosList 
for positionnbr in CodonPosList.keys():
	for mycodondatanbr in range(0, len(CodonPosList[positionnbr])):
		enrichment = CodonPosList[positionnbr][mycodondatanbr][2] / CodonPosList[positionnbr][mycodondatanbr][1]
		CodonPosList[positionnbr][mycodondatanbr] = [CodonPosList[positionnbr][mycodondatanbr][0], CodonPosList[positionnbr][mycodondatanbr][1],CodonPosList[positionnbr][mycodondatanbr][2], enrichment]


# Get the codons of the WT sequence specified in the options file 
WTposvscodon=get_codons(WTseq)

# Calculate and retain the lower and the higher enrichment values. The lower enrichment value will be taken as reference for normalizations 
MinEnrich = 0
MaxEnrich = 0
for pos in WTposvscodon:
	if pos[0] in CodonPosList.keys():
		for mycodondatanbr in range(0, len(CodonPosList[pos[0]])):
			if MinEnrich == 0 or MinEnrich > CodonPosList[pos[0]][mycodondatanbr][3]:
				MinEnrich = CodonPosList[pos[0]][mycodondatanbr][3]
			if MaxEnrich == 0 or MaxEnrich < CodonPosList[pos[0]][mycodondatanbr][3]:
				MaxEnrich = CodonPosList[pos[0]][mycodondatanbr][3]
print("Lower codon enrichment value = %s\nHigher codon enrichment value = %s" % (MinEnrich, MaxEnrich))
logfile.write("Lower codon enrichment value = %s\nHigher codon enrichment value = %s\n"  % (MinEnrich, MaxEnrich))

# Generate an output file including WT codons
myoutfile=open(OutFile, "wt")
myoutfile.write('Position,Codon,Residue,frequency_1,frequency_2,enrichment\n')
for pos in WTposvscodon:
	if pos[0] in CodonPosList.keys():
		for mycodondatanbr in range(0, len(CodonPosList[pos[0]])):
			myoutfile.write('%i,%s, %s, %0.8f,%0.8f,%0.8f\n' % (pos[0], pos[1]  + "=>" + CodonPosList[pos[0]][mycodondatanbr][0], translate_dna(pos[1]) + str(pos[0]+35) + translate_dna(CodonPosList[pos[0]][mycodondatanbr][0]), CodonPosList[pos[0]][mycodondatanbr][1], CodonPosList[pos[0]][mycodondatanbr][2], CodonPosList[pos[0]][mycodondatanbr][3]))
myoutfile.close()
print("Finished writing the output file: %s" % (OutFile))

# Generate an output file excluding WT codons. #Also generate an input file for CodonLogo in the Clustal format.
myoutfile=open(OutFilenoWT, "wt")
#myWeblogooutfile=open("CodonsForWeblogo.ali", "wt")
myoutfile.write('Position,Codon,Residue,frequency_1,frequency_2,enrichment\n')
for pos in WTposvscodon:
	if pos[0] in CodonPosList.keys():
		for mycodondatanbr in range(0, len(CodonPosList[pos[0]])):
			if pos[1] != CodonPosList[pos[0]][mycodondatanbr][0]:
				myCurrentSeq = ""				
				for n in range(1,pos[0]):
					myCurrentSeq = myCurrentSeq + "---"
				myCurrentSeq = myCurrentSeq + CodonPosList[pos[0]][mycodondatanbr][0]
				for m in range(pos[0],len(CodonPosList)):
					myCurrentSeq = myCurrentSeq + "---"
#				print (myCurrentSeq) 
#				print(int(enrichment / MinEnrich))
#				for copynbr in range(1, int(enrichment / MinEnrich)):
#					myWeblogooutfile.write('CodonMutation=%s:CodonPosition=%i:ProtMutation=%s:Copy=%i\t%s\n'% (pos[1]  + "=>" + CodonPosList[pos[0]][mycodondatanbr][0], pos[0], translate_dna(pos[1]) + str(pos[0]+35) + translate_dna(CodonPosList[pos[0]][mycodondatanbr][0]),copynbr,myCurrentSeq))
				myoutfile.write('%i,%s, %s, %0.8f,%0.8f,%0.8f\n' % (pos[0], pos[1]  + "=>" + CodonPosList[pos[0]][mycodondatanbr][0], translate_dna(pos[1]) + str(pos[0]+35) + translate_dna(CodonPosList[pos[0]][mycodondatanbr][0]), CodonPosList[pos[0]][mycodondatanbr][1], CodonPosList[pos[0]][mycodondatanbr][2], CodonPosList[pos[0]][mycodondatanbr][3]))
			sys.stdout.write("Codon position %i of %i | Retained codon %i of %i done! \r" % (pos[0], len(CodonPosList), mycodondatanbr, len(CodonPosList[pos[0]])))
			sys.stdout.flush()
myoutfile.close()
#myWeblogooutfile.close()
print("Finished writing the output file: %s" % (OutFilenoWT))
#print("Finished writing the output file: %s" % ("CodonsForWeblogo.ali"))

# Generate an output file excluding WT codons and containing for each codon the enrichment value per position.
mycodonlist = ['ACC', 'ATG', 'ACA', 'ACG', 'ATC', 'AAC', 'ATA', 'AGG', 'CCT', 'CTC', 'AGC', 'AAG', 'AGA', 'CAT', 'AAT', 'ATT', 'CTG', 'CTA', 'ACT', 'CAC', 'AAA', 'CCG', 'AGT', 'CCA', 'CAA', 'CCC', 'TAT', 'GGT', 'TGT', 'CGA', 'CAG', 'CGC', 'GAT', 'CGG', 'CTT', 'TGC', 'GGG', 'TAG', 'GGA', 'TAA', 'GGC', 'TAC', 'GAG', 'TCG', 'TTT', 'GAC', 'CGT', 'GAA', 'TCA', 'GCA', 'GTA', 'GCC', 'GTC', 'TGA', 'GCG', 'GTG', 'TTC', 'GTT', 'GCT', 'TTA', 'TTG', 'TCC', 'TGG', 'TCT']
myoutfile=open("enrichment.txt", "wt")
myenrichDict = {}

# create dictionary using codons as keys
for codon in mycodonlist:
	initialization_list = []
	for pos in WTposvscodon:
		initialization_list = initialization_list + [0]
	myenrichDict[codon] = initialization_list
for codon in mycodonlist:
	for pos in WTposvscodon:
		if pos[0] in CodonPosList.keys():
			for mycodondatanbr in range(1, len(CodonPosList[pos[0]])):
				if codon != pos[1] and codon == CodonPosList[pos[0]][mycodondatanbr][0]:
					myenrichDict[codon][pos[0]-1] = CodonPosList[pos[0]][mycodondatanbr][3]
#				sys.stdout.write("Codon position %i of %i | Retained codon %i of %i done! \r" % (pos[0], len(CodonPosList), mycodondatanbr, len(CodonPosList[pos[0]])))
#				sys.stdout.flush()
#sys.stdout.write("\n")

# generate the heading line for the output	
myoutfile.write("Codon")
for posnbr in range(0, len(myenrichDict["ACC"])):
	myoutfile.write(' Pos_%i' % (posnbr+1))
myoutfile.write("\n")

# output results to file line-by-line
for codon in mycodonlist:
	if codon in myenrichDict:
		myoutfile.write('%s ' % (codon))
		for positionvalue in myenrichDict[codon]:
			myoutfile.write('%0.8f ' % (positionvalue))
		myoutfile.write('\n')
myoutfile.close()

print("Finished writing the output file: %s" % ("enrichment.txt"))

	
logfile.write("Done writing all results: %s\n"  % (str(datetime.datetime.now())))
EndTime = datetime.datetime.now()
print("Batch job successfully ended : " + str(EndTime))
logfile.write("Batch job successfully ended : " + str(EndTime) + "\n")
print("Duration : " + str(EndTime - StartTime))
logfile.write("Duration : " + str(EndTime - StartTime) + "\n")
logfile.close()
