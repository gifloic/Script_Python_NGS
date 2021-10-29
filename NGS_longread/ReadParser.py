#!/usr/bin/env python3

import lr_tools as lt

import sys
import os
import re

from datetime import datetime
import time
import tempfile
import argparse
import pickle
import fileinput
from pprint import pprint
from copy import deepcopy
from mpi4py import MPI
from itertools import islice

import fileinput as FI
import h5py

try:
    import configparser
except:
    import ConfigParser as configparser
from configparser import ExtendedInterpolation



class Parser:
    def __init__(self, fastq,
                 exp_index,
                 exp_orientation,
                 config_file,
                 search_barcode,
                 use_mpi,
                 outputfile,
                 verbose = False, cleanup = True):
        self.fastq = fastq
        self.use_mpi = use_mpi
        self.exp_index = exp_index
        self.exp_orientation = exp_orientation
        self.config_file = config_file
        self.search_barcode = search_barcode # can be 0, -1, 1 depending if nobarcode, reverse barcode or formard barcode has to be found
        self.verbose = verbose
        self.cleanup = cleanup


        self.d3to1, self.d1to3 = lt.get_dicaa()

        self.start_mpi()
        if self.rank == 0:
            print(">>> Counting number of lines")
        self.number_of_lines = self.get_nlines()
        self.read_config_file()

        if self.rank == 0:
            print(">>> Extracting the expected mutation patterns")
        self.get_dic_expect_mutations()

        if self.rank == 0:
            print("CONSTANT_PRIMER:",self.constant_primer)
            print("PATTERN_BARCODE:",self.pattern_barcode)

        if self.rank == 0:
            print(">>> Parsing the fastq")
        self.split_input(outputfile, rank=self.rank, size=self.size)


    def start_mpi(self):
        if self.use_mpi:
            # MPI communication
            self.comm = MPI.COMM_WORLD
            self.rank = self.comm.Get_rank()
            self.size = self.comm.Get_size()
        else:
            self.comm = None
            self.rank = 0
            self.size = 1

    def get_nlines(self):
        if self.rank == 0:
            number_of_lines = lt.count_file_lines(self.fastq)
            #print("I finished")
        else:
            number_of_lines = None
        if self.use_mpi:
            self.comm.Barrier()
            number_of_lines = self.comm.bcast(number_of_lines, root=0)
        return number_of_lines

    def get_dic_expect_mutations(self):
        """

        :return:
        """
        self.dmut = lt.parse_twist_devis_table(self.mutfile, self.dn2a, self.d3to1, start=None, stop=None, info_codon_index=None)

    def split_input(self, outputfile, rank=0, size=1):
        """Split an input with Nlines so that every rank creates its fraction and works on it
        """
        if rank == 0:
            print("Splitting the input file into {} parts".format(size))

        self.outputfile = open(outputfile + ".rank{}".format(rank), 'w')
         # Ditributing the work for every rank
        full_input = self.fastq
        size_chunks = (int(self.number_of_lines) / 4) / self.size + 1 # we want to recover blocks of 4 lines
        index_start = int(rank * (size_chunks * 4))
        if rank == size - 1:
            index_stop = None
        else:
            index_stop = int((rank + 1) * (size_chunks * 4))
        if self.verbose:
            print("!! rank {}, will parse lines from {} to {}.".format(rank, index_start, index_stop))

        hit = 0
        read_index = 0
        rank_index = 10 + rank
        with open(full_input) as file_obj:
            n = 1
            object_to_iterate = islice(file_obj, index_start, index_stop)
            # for line in islice(file_obj, index_start, index_stop):
            for line in object_to_iterate:
                # Parse the fastq file
                l = line[: -1]
                if l[0] == "@":
                    n=1
                    FLAG_REACH_END_READ = False
                    read_index += 1
                    rdint = int('%d%09d' % (rank_index, read_index))
                    rdhex = hex(rdint)
                    header = l
                    n += 1
                    KEEP_READ = False
                    PRIMER_DEFECT = False
                    BARCODE_DEFECT = False
                elif n == 4:
                    n += 1
                elif n == 3:
                    n += 1
                elif n == 2:
                    n += 1
                    # extract umi
                    seq = line.strip()
                    if self.exp_orientation == 'reverse':
                        seq = lt.get_complement(lt.get_reverse(seq))
                    if self.exp_orientation == 'forward':
                        pattern1 = r"[NACGT]{{1,{}}}".format(self.length_tail) + self.constant_primer + "([ATGC]+)"
                        match = re.search(pattern1, seq)
                    elif self.exp_orientation == 'reverse':
                        pattern1 = "([ATGC]+)" + self.constant_primer + r"[NACGT]{{1,{}}}".format(self.length_tail)
                        match = re.search(pattern1, seq)
                    if match:  # primer is detected
                        KEEP_READ = True
                        start_seq, stop_seq = match.span(1)
                        if self.exp_orientation == 'forward':
                            start_read_coding_seq = self.stop_primer_in_ampliseq - self.start_codseq_in_ampliseq
                        if self.exp_orientation == 'reverse':
                            start_read_coding_seq = len(self.coding_seq) - len(seq[:(stop_seq -
                                                                            (self.start_primer_in_ampliseq
                                                                             - self.stop_codseq_in_ampliseq))])

                        if self.verbose:
                            print("First filter for constant primer detection passed")
                            print(start_seq, stop_seq,
                                                  '\n',seq[start_seq:start_seq+10],'...',seq[stop_seq-10:stop_seq])
                            print(start_read_coding_seq,start_read_coding_seq%3,start_read_coding_seq//3)
                            print("len_coding_seq:",len(self.coding_seq))
                            print("####\n")
                            print(self.coding_seq[start_read_coding_seq:start_read_coding_seq+10])
                            print(seq[:10])

                        if self.search_barcode:
                            if self.exp_orientation == 'forward':
                                # The barcode is expressed 5'3' downstream the constant part
                                match_barcode = re.search(self.pattern_barcode, seq[start_seq:])
                            elif self.exp_orientation == 'reverse':
                                # The barcode is expressed 5'3' upstream the constant part
                                match_barcode = re.search(self.pattern_barcode, seq[:stop_seq])
                            if match_barcode:
                                umi_all_positions = match_barcode.group(1)
                                if self.barcode_pos2remove:
                                    umi = ""
                                    for ii,pos in enumerate(umi_all_positions):
                                        if ii in self.barcode_pos2remove:
                                            continue
                                        umi += pos
                                else:
                                    umi = umi_all_positions
                                KEEP_READ = True
                            else:
                                KEEP_READ = False
                                BARCODE_DEFECT = True
                        else:
                            umi = None
                    else:
                        KEEP_READ = False
                        PRIMER_DEFECT = True
                    #
                    # Get the mutation code in nucleotides
                    #
                    print(self.coding_seq[174:180])

                    if KEEP_READ:
                        if self.exp_orientation == 'forward':
                            mcode = lt.get_mutation_in_read(seq[start_seq:],
                                                        self.coding_seq,
                                                        start_read_coding_seq,
                                                        thresh_multi_mut = self.thresh_multi_mut)
                        if self.exp_orientation == 'reverse':
                            mcode = lt.get_mutation_in_read(seq[:stop_seq],
                                                        self.coding_seq,
                                                        start_read_coding_seq,
                                                        thresh_multi_mut = self.thresh_multi_mut)
                        if self.verbose:
                            print("## Keep:{}".format(header))
                            print("UMI:",umi)
                            print("SEQ:",seq[start_seq:start_seq+10])
                            print("mcode:",mcode)
                            print(lt.get_aamut(self.dn2a, self.coding_seq.strip(), mcode))
                            print("\n")
                    else:
                        if self.verbose:
                            print("## Discard:{}\n".format(header))
                    #

                    if self.search_barcode:
                        if KEEP_READ:
                            if mcode in self.dmut: # for instance compatible with twist or trinex mutational pattern
                                mcode = "O\t" + mcode
                            else:
                                mcode = "N\t" + mcode
                            self.outputfile.write("{}\t{}\t{}\n".format(rdhex, mcode, umi))
                        else:
                            if PRIMER_DEFECT:
                                self.outputfile.write("{}\tD\tP\n".format(rdhex))
                            if BARCODE_DEFECT:
                                self.outputfile.write("{}\tD\tB\n".format(rdhex))
                    else:
                        if KEEP_READ:
                            if mcode in self.dmut:
                                mcode = "O\t" + mcode
                            else:
                                mcode = "N\t" + mcode
                            self.outputfile.write("{}\t{}\t-\n".format(rdhex, mcode))
                        else:
                            self.outputfile.write("{}\tD\t-\n".format(rdhex))
        if rank == 0:
            print(
                "Start of the mutant analysis at position in the coding sequence:{}".format(start_read_coding_seq))
            print("First codon analyzed indexed as aminoacid: {}".format(start_read_coding_seq // 3 + 1))
        self.outputfile.close()
        if rank == 0:
            listfiles = []
            for r in range(size):
                listfiles.append(outputfile + ".rank{}".format(r))
            with open(outputfile,"w") as fout, fileinput.input(listfiles) as fin:
                for line in fin:
                    fout.write(line)
            for fname in listfiles:
                os.remove(fname)


    def read_config_file(self):
        """

        :return:
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)
        print(self.config_file)
        print(open(self.config_file).read())


        reverse_primers = config.get("Parameters", "REVERSE").split()
        forward_primers = config.get("Parameters", "FORWARD").split()
        if self.exp_orientation   == 'forward':
            self.constant_primer  = forward_primers[int(self.exp_index)-1]
        elif self.exp_orientation == 'reverse':
            self.constant_primer  = lt.get_complement(lt.get_reverse(reverse_primers[int(self.exp_index)-1]))
        else:
            print("Error in defining the constant primer. Please check")
            sys.exit(1)

        #
        # Parameters for barcode retrieving
        #
        self.barcode             = config.get("Parameters", "BARCODE")
        if self.search_barcode  == 1:
            self.pattern_barcode = lt.encode_barcode(self.barcode)
        elif self.search_barcode  == -1:
            # Since barcode is defined in the 3'5', we reverse it because all reads are converted 5'3'
            self.pattern_barcode = lt.encode_barcode(lt.get_complement(lt.get_reverse(self.barcode)))
        else:
            self.pattern_barcode = None
        try:
            self.barcode_pos2remove  = [int(x) for x in config.get("Parameters", "POS2REMOVE").split(',')]
        except:
            self.barcode_pos2remove  = None

        ndephaser          = int(config.get("Parameters", "DEPHASER_MAX_LENGTH"))
        ntolerance         = int(config.get("Parameters", "DEPHASER_TOLERANCE"))
        self.length_tail   = ndephaser + ntolerance

        self.coding_seq    = lt.get_seq_from_fasta(config.get("Files", "CODING_SEQ"))
        self.amplic_seq    = lt.get_seq_from_fasta(config.get("Files", "AMPLICON_SEQ"))

        self.mutfile       = config.get("Files", "TWIST_MUT_TABLE")
        faa2nuc            = config.get("Files", "AA2NUC")
        self.da2n = pickle.load(open(faa2nuc,'rb'),encoding='bytes') # due to py2 to 3
        self.da2n = lt.convert_dic_bytes2string(self.da2n) # due to py2 to 3
        fnuc2aa   = config.get("Files", "NUC2AA")
        self.dn2a = pickle.load(open(fnuc2aa,'rb'))#,encoding='bytes')

        mcodseq2ampli = re.search(self.coding_seq, self.amplic_seq)
        self.start_codseq_in_ampliseq, self.stop_codseq_in_ampliseq = mcodseq2ampli.span()

        if self.exp_orientation == 'forward':
            mprimer2ampli = re.search(self.constant_primer, self.amplic_seq)
        if self.exp_orientation == 'reverse':
            mprimer2ampli = re.search(self.constant_primer,self.amplic_seq)
        self.start_primer_in_ampliseq, self.stop_primer_in_ampliseq = mprimer2ampli.span()
        print("###", self.start_primer_in_ampliseq, self.stop_primer_in_ampliseq)
        #nuc2abundance = config.get("Files", "NUC2ABUND")
        self.thresh_multi_mut = int(config.get("Parameters", "THRESH_MULTI_MUT"))



def Main():

    try:
        default_configfile = os.path.join(os.path.dirname(os.path.abspath(os.readlink(os.path.abspath(__file__)))),
                                  "NGS_lread.ini")
    except:
        default_configfile = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(os.path.abspath(__file__)))),
                                  "NGS_lread.ini")


    usage = "Script to parse the reads from a long read experiment.\n" + \
            "output is a dictionary encoded as a hdf5 file"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('-f', '--fastq', action="store", help='input fastq file')
    parser.add_argument('-o', '--outputfile', action="store", help='output mutant analysis file')
    parser.add_argument('-i', '--exp_index', action="store", help='index (1, 2, 3, ...) of the experiment defines the constant primer used in the PCR')
    parser.add_argument('-e', '--exp_orientation', choices = ['forward', 'reverse'], action="store", help='experiment is reverse or forward')
    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('-b', '--barcode', default=False, action="store_true", help='search barcode as in configuration file')
    parser.add_argument('-r', '--reverse_barcode', default=False, action="store_true", help='search reverse complementary of configuration file')
    parser.add_argument('--mpi', action="store_true", default=False, help='triggers mpi running')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.barcode:
        search_barcode = 1
    elif args.reverse_barcode:
        search_barcode = -1
    else:
        search_barcode = 0

    if args.fastq:
        """
        Here generate the code for generating the hhr file first
        """
        P = Parser(args.fastq,
                   args.exp_index,
                   args.exp_orientation,
                   args.config_file,
                   search_barcode,
                   args.mpi,
                   args.outputfile,
                   verbose = args.verbose,
                   cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
