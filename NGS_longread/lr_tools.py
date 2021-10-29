#!/usr/bin/env python3

import os,re
import subprocess

def encode_barcode(s):
    """
    :param s: string of DNA
    :return: barcode ancoded as a regular expression pattern
    """
    """ 
    """
    nucD = {}
    nucD["("] = "("
    nucD[")"] = ")"
    nucD["A"] = "A"
    nucD["C"] = "C"
    nucD["G"] = "G"
    nucD["T"] = "T"
    nucD["U"] = "U"
    nucD["R"] = "AG"
    nucD["Y"] = "CT"
    nucD["S"] = "GC"
    nucD["W"] = "AT"
    nucD["K"] = "GT"
    nucD["M"] = "AC"
    nucD["B"] = "CGT"
    nucD["D"] = "AGT"
    nucD["H"] = "ACT"
    nucD["V"] = "ACG"
    nucD["N"] = "ACGT"

    pattern = ""
    for n in s:
        if len(nucD[n]) == 1:
            pattern += nucD[n]
        else:
            pattern += "[{}]".format(nucD[n])
    return pattern

def get_dicaa():
    d3to1 = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',\
             'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',\
             'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',\
             'TRP':'W','TYR':'Y','TER':'Z'}
    d1to3 = {}
    for key,val in d3to1.items():
        d1to3[val] = key
    return d3to1, d1to3

def get_reverse(s):
    """

    :param s: anysequence
    :return: the reversed order sequence
    """
    reversed_sequence = ""
    x = s[::-1]
    for c in x:
        if c=="(":
            reversed_sequence += ")"
        elif c==")":
            reversed_sequence += "("
        else:
            reversed_sequence += c
    return reversed_sequence


def get_complement(s):
    """
    :param s: string of DNA
    :return: complement
    """
    """ Returns complementary DNA sequence even with ambiguous nucleotides.
    """
    nucD = {}
    nucD["\n"] = ""
    nucD["("] = "(" # useful for regular expression tracking
    nucD[")"] = ")" # useful for regular expression tracking
    nucD["A"] = "T"
    nucD["C"] = "G"
    nucD["G"] = "C"
    nucD["T"] = "A"
    nucD["U"] = "A"
    nucD["R"] = "Y"
    nucD["Y"] = "R"
    nucD["S"] = "S"
    nucD["W"] = "W"
    nucD["K"] = "M"
    nucD["M"] = "K"
    nucD["B"] = "V"
    nucD["D"] = "H"
    nucD["H"] = "D"
    nucD["V"] = "B"
    nucD["N"] = "N"
    complement = ""
    for n in s:
        complement += nucD[n]
    return complement

def count_file_lines(input_file):
    """

    :param input_file:
    :return: integer with number of lines
    """
    return int(subprocess.check_output(["wc","-l","{}".format(input_file)]).split()[0])

def convert_list_bytes2string(l):
    x = []
    for e in l:
        if isinstance(e,bytes):
            try:
                e = e.decode()
            except:
                e = ' '
        elif isinstance(e, list):
            e = convert_list_bytes2string(e)
        x.append(e)
    return x

def convert_dic_bytes2string(d):
    """
    Python Code to convert ByteString key:value
        pair of dictionary to String.
    :param d: dictionnary of bytesString (often got from unpickling python2 dict)
    :return: dictionnary with String format in key:values
    """
    x = {}
    for key, value in d.items():
        if isinstance(key,bytes):
            key = key.decode()
        if isinstance(value,bytes):
            value = value.decode()
        elif isinstance(value, dict):  # recursivity if dict of dict
            value = convert_dic_bytes2string(value)
        elif isinstance(value,list): # recursivity if list
            value = convert_list_bytes2string(value)
        x[key] = value
    return x

def get_seq_from_fasta(fasta_file):
    """

    :param file:
    :return:
    """
    fin = open(fasta_file).readlines()
    seq = ""
    for l in fin:
        if l[0] == ">":
            header = l[1:].strip()
        else:
            if len(l.split()) > 0:
                seq += l.strip("\n")
    return seq

def get_mutation_in_read(read, coding_seq, index,  thresh_multi_mut=10):
    """

    :param read:
    :param coding_seq:
    :param index:
    :return:
    """
    len_coding_seq = len(coding_seq)
    mutation_code = ""
    if index + len(read) <= len(coding_seq):
        len_read_limit = len(read)
    else:
        len_read_limit = len(coding_seq) - index
    if read[:len_read_limit] == coding_seq[index:index+len(read[:len_read_limit])]:
        return "W"

    separator = ''
    FLAG_MUTANT = 0
    for ii,n in enumerate(read):
        index_in_coding_seq = index + ii
        if index_in_coding_seq >= len_coding_seq:
            break
        elif n == coding_seq[index_in_coding_seq]:
            pass
        else:
            FLAG_MUTANT += 1
            mutation_code += separator
            mutation_code += str(index_in_coding_seq) + n
            if separator == '_':
                FLAG_MUTANT += 1
            if index_in_coding_seq%3 != 2:
                separator = '.'

        if FLAG_MUTANT and index_in_coding_seq%3 == 2:
            separator = '_'
        if FLAG_MUTANT > thresh_multi_mut:
            return "T"

    return mutation_code

def get_aamut(dn2a, coding_seq, mutation_code, aa1letter = True):
    """
    Convert into aminoacids code the nucleotide code of the mutations
    :param dn2a:
    :param coding_seq:
    :param mutation_code:
    :return: list_aamut, list_nucmut
    """
    codons = mutation_code.split('_')
    list_aamut = []
    list_nucmut = []
    if aa1letter:
        d3to1, d1to3 = get_dicaa()
    for codon in codons:
        nmuts = codon.split('.')
        try:
            mut_index = [int(x[:-1]) for x in nmuts] # as [393,394]
        except ValueError:
            return None,None
        mut_nuc   = [x[-1] for x in nmuts] # as ['G','T']
        dmut = dict(zip(mut_index,mut_nuc))
        start_codon = 3*mut_index[0]//3
        wt_codon = coding_seq[start_codon:start_codon+3]
        codon_index = range(start_codon,start_codon+3)
        mut_codon = ""
        for ii,n in enumerate(codon_index):
            if n in dmut:
                mut_codon += dmut[n]
            else:
                mut_codon += wt_codon[ii]
        aa_wt = dn2a[wt_codon]
        if aa1letter:
            aa_wt = d3to1[aa_wt]
        aa_mut = dn2a[mut_codon]
        if aa1letter:
            aa_mut = d3to1[aa_mut]
        aa_mutname = aa_wt + str(mut_index[0]//3) + aa_mut
        nuc_mutname = wt_codon + str(start_codon) + mut_codon
        list_aamut.append(aa_mutname)
        list_nucmut.append(nuc_mutname)
    return list_aamut, list_nucmut


def parse_twist_devis_table(twist_file, dn2a, d3to1, start=None, stop=None, info_codon_index=None):
    """
    :param twist_file: csv table generated for the order at Twist
    :param dn2a: dic extracted from the nuc2aa.pickle
    :param d3to1: dic to translate 3-letter code aa into 1
    :param start: residue index of the beginning of the sequence in the csv file
    :param stop: residue index of the end of the sequence in the csv file
    :param info_codon_index:
    :return: dictionnary integrating all the information about the mutation landscape
            dic_mutnuc[mut_nuc_name] mut_nuc_name as "125G.126T.127T" or "125C.127A"
            dic_mutnuc[mut_nuc_name]['nuc_wt'] -> wt_nucleotide
            dic_mutnuc[mut_nuc_name]['AA_index'] -> resi
            dic_mutnuc[mut_nuc_name]['AA_wt'] -> wt_aminoacid
            dic_mutnuc[mut_nuc_name]['nuc_mut'] -> mutated_nucleotide
            dic_mutnuc[mut_nuc_name]['AA_mut'] -> mutated amino acid
             dic_mutnuc[mut_nuc_name]['code_type'] -> Either '30'->WT AA with 3 nuc mutations,
                                                or      '20'->WT AA with 2 nuc mutations,
                                                or      '10'->WT AA with 1 nuc mutations,
                                                or mutated AA with '3', '2' or '1' nucleotide
           * Historically full information is given by:
            dic_mutnuc[mut_nuc_name]['full_mut_AA_name'] -> 'nuc_mutation@@aa_mutation'
            * Shorter information (no wt nuc) is given by:
            dic_mutnuc[mut_nuc_name]['short_mut_AA_name'] -> 'nuc_mutation@@aa_mutation'


    Variables defined in the script but not returned
        list_aai -> Asf1 absolute index in structure
        dic_mut[<Asf1 absolute index in structure>]['i2d'] -> <Asf1 index from 0 to n>
        dic_ftwist[<first column of the twist file=legend or codon>] -> line following the 1st col
    """
    ftwist = open(twist_file).readlines()
    print("Reading the mutation profile table : ",twist_file)
    dic_mut = {}
    dic_ftwist = {}

    dic_mutnuc = {}
    set_mutnuc = set()
    list_mut_nuc_name = []

    if start:
        START = start
    if stop:
        STOP = stop
    if info_codon_index:
        CODON_INDEX_INFO_LINE = info_codon_index

    if not start:
        START = 1
    print("The first index of the AA considered in mutation table is :", START)

    set_possible_codons = []  # set of list[codon,AA_3letter]
    FOUND_AAindex = False
    for l in ftwist:
        s = [x.strip() for x in l.split(';')]
        if s[0] == '':
            continue
        dic_ftwist[s[0]] = s[1:]
        if len(s[0]) == 3:
            set_possible_codons.append([s[0], s[1]])
        if s[0] == 'AA_index':
            FOUND_AAindex = True
            list_aai = s[2:] # List of indexes matching coding sequence
            for ii in range(1, len(list_aai) + 1):
                dic_mut[list_aai[ii - 1]] = {}
                # matches REAL resi index to index d which is used to retrieve the information in the other columns
                dic_mut[list_aai[ii - 1]]['i2d'] = ii
    if not FOUND_AAindex:
        print("The mutation profile table file is lacking the AA_index line. Check the correct format !")
    if not stop:
        STOP = int(list_aai[ii - 1])

    print("The last index of the AA considered in mutation table is :",STOP)

    # Definition for IP scanned region
    for resi in range(START, STOP + 1):
        ini_nuc = dic_ftwist['NUC_wt'][dic_mut[str(resi)]['i2d']]
        # ini_nuc is wt codon of every mutated position resi
        try:
            codon_index = int(dic_ftwist[CODON_INDEX_INFO_LINE][dic_mut[str(resi)]['i2d']])
        except:
            codon_index = resi*3

        n = 0
        for key in dic_ftwist:
            # key is <first column of the twist file>
            if len(key) != 3:
                continue
            """
            We read only the codons
            """
            # for every resi we get the equivalent column index and get whether codon is mutated or not
            mut_nuc = dic_ftwist[key][dic_mut[str(resi)]['i2d']]
            if mut_nuc == '':
                # codon not mutated
                continue
            n += 1
            list_mutnuc = []
            list_mutnuc_withwt = []
            #ini_nuc is wt codon
            for jj, nn in enumerate(ini_nuc):
                if nn != mut_nuc[jj]:
                    nuc_index = codon_index + jj
                    list_mutnuc.append('{}{}'.format(nuc_index, mut_nuc[jj]))
                    list_mutnuc_withwt.append('{}{}{}'.format(nn, nuc_index, mut_nuc[jj]))
            if len(list_mutnuc) > 0:
                mut_nuc_name = '.'.join(list_mutnuc)
                set_mutnuc.add(mut_nuc_name)
                dic_mutnuc[mut_nuc_name] = {}
                list_mut_nuc_name.append(mut_nuc_name)
                dic_mutnuc[mut_nuc_name]['nuc_wt'] = ini_nuc
                dic_mutnuc[mut_nuc_name]['AA_index'] = resi
                dic_mutnuc[mut_nuc_name]['AA_wt'] = dic_ftwist['AA_wt'][dic_mut[str(resi)]['i2d']]
                dic_mutnuc[mut_nuc_name]['nuc_mut'] = mut_nuc
                dic_mutnuc[mut_nuc_name]['AA_mut'] = dn2a[mut_nuc]
                dic_mutnuc[mut_nuc_name]['type'] = dn2a[mut_nuc]
                label1 = '.'.join(list_mutnuc_withwt)
                label2 = '{}{}{}'.format(d3to1[dic_mutnuc[mut_nuc_name]['AA_wt']], resi,
                                         d3to1[dic_mutnuc[mut_nuc_name]['AA_mut']])
                if dic_mutnuc[mut_nuc_name]['AA_wt'] == dic_mutnuc[mut_nuc_name]['AA_mut']:
                    dic_mutnuc[mut_nuc_name]['code_type'] = '{}0'.format(len(list_mutnuc))
                else:
                    dic_mutnuc[mut_nuc_name]['code_type'] = '{}'.format(len(list_mutnuc))

                dic_mutnuc[mut_nuc_name]['full_mut_AA_name'] = '@@'.join([label1, label2])
                dic_mutnuc[mut_nuc_name]['short_mut_AA_name'] = '@@'.join([mut_nuc_name,label2])

            else:
                print('problem with {},{}'.format(resi, key))
    print('The number of mutants in the sequenced region should be : {}'.format(len(set_mutnuc)))
    dic_mutnuc['order'] = list_mut_nuc_name[:]
    return dic_mutnuc

