#!/usr/bin/env python3

import lr_tools as lt

import sys
import os
import re

from datetime import datetime
import time
import tempfile
import subprocess
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



class Process:
    def __init__(self, input_name, primer_range, conditions, working_dir, nameout, mpi, config_file, verbose = False, cleanup = False):
        self.input_name = input_name
        self.primer_range = [int(x) for x in primer_range.split('-')]
        self.conditions = [int(x) for x in conditions.split(',')]
        self.working_dir = working_dir
        self.nameout = nameout
        self.mpi = mpi
        self.config_file = config_file
        self.verbose  = verbose
        self.cleanup  = cleanup

        self.read_config_file()
        self.run_process()

    def run_process(self):

        os.chdir(self.working_dir)
        for read_orientation in ['forward', 'reverse']:
            if len(self.conditions) > 1:

                self.comparecond_out = "{}_{}_compare_conditions_{}.out".format(self.nameout,read_orientation, 'vs'.join([str(x) for x in self.conditions[::-1]]))
            self.dmut_cond = {}
            for condition in self.conditions:
                self.fileout = "{}_{}_{}_barcodes.out".format(self.nameout,read_orientation,condition)
                self.mutant_out = "{}_{}_{}_mutations.out".format(self.nameout,read_orientation,condition)
                self.dbcode = {}
                self.dmutcode = {}
                self.dmutcode['order'] = []
                for primer in range(self.primer_range[0],self.primer_range[1]+1):

                    print("Analyzing primer ",primer)
                    process_path = '{}_{}_{}'.format(self.input_name,primer,condition)
                    fbcode = '{}_{}_assembled_{}.out'.format(self.input_name,primer,read_orientation)
                    path_fcode = os.path.join(process_path,fbcode)
                    mutcode = '{}_{}_assembled_mutants_{}.out'.format(self.input_name,primer,read_orientation)
                    path_mutcode = os.path.join(process_path,mutcode)


                    # Analyze the common barcodes between all the experiements
                    # store in self.dbcode dictionnary for every primer the name of the corresponding mutant
                    with (open(path_fcode)) as fin:
                        for line in fin:
                            s = line.split()
                            bc = s[0]
                            mut1 = ' '.join(s[1:4])
                            if bc not in self.dbcode:
                                self.dbcode[bc] = {}
                            self.dbcode[bc][primer] = mut1

                    # Analyze the mutant counts for every primer output
                    with (open(path_mutcode)) as fin:
                        for line in fin:
                            s = line.split()
                            nucmut = s[0]
                            aamut = s[1]
                            aalist = [int(s[1][1:-1]),s[1][-1]]
                            count = s[2]
                            distrib = s[3]
                            if nucmut not in self.dmutcode:
                                self.dmutcode[nucmut] = {}
                                self.dmutcode['order'].append([aalist[0],aalist[1],nucmut])
                            self.dmutcode[nucmut][primer] = [aamut,count,distrib]
                self.dmutcode['order'].sort()
                #pprint(self.dmutcode['435A.436A.437C'])
                # Output file gathering tog ether all the barcode and their mutants assigned
                with (open(self.fileout,"w")) as fout:
                    fout.write("# COLUMNS = 1) Barcode\t2) Primer 1 <Counts ExpectCode VariantCode> | 3) Primer 2 <...> | ...\n#\n")
                    fout.write("# ExpectCode : W = wild-type ; N = Not-expected-mutation ; O = Okay, expected\n")
                    fout.write("# VariantCode : W = wild-type ; T = Number of mutations above Threshold ;\n"
                               "#               B = barcode not found ; P = Constant primer not found\n#\n")
                    for bc in self.dbcode:
                        fout.write("{}\t".format(bc))
                        for primer in range(self.primer_range[0],self.primer_range[1]+1):
                            if primer in self.dbcode[bc]:
                                mut = self.dbcode[bc][primer]
                                fout.write("{} | ".format(mut))
                            else:
                                fout.write("- - - | ")
                        fout.write("\n")

                primer_start = self.primer_range[0]
                primer_end = self.primer_range[1]
                with (open(self.mutant_out,"w")) as fout:
                    fout.write("# COLUMNS = 1) Nucleotide Mutant\n"
                               "#           2) AminoAcid Mutant\n"
                               "#           3) Mutant counts for every primer (from {} to {})\n"
                               "#           4) Distrib counts per bcodes (after pipe from {} to {})\n"
                               .format(primer_start,primer_end, primer_start, primer_end))
                    fout.write("#\n")
                    for mutlist in self.dmutcode['order']:
                        mut = mutlist[-1]
                        if mut not in self.dmut_cond:
                            self.dmut_cond[mut] = {}
                        self.dmut_cond[mut][condition] = []
                        fout.write("{}\t".format(mut))
                        first_primer = True
                        for primer in range(self.primer_range[0],self.primer_range[1]+1):
                            [aamut, count, distrib] = self.dmutcode[mut][primer]
                            self.dmut_cond[mut][condition].append(count)
                            if first_primer:
                                fout.write("{}\t".format(aamut))
                                first_primer = False
                            fout.write("{} ".format(count))
                        fout.write("|\t")
                        for primer in range(self.primer_range[0],self.primer_range[1]+1):
                            [aamut, count, distrib] = self.dmutcode[mut][primer]
                            fout.write("{} ".format(distrib))
                        fout.write("\n")

            with(open(self.comparecond_out,"w")) as fout:
                fout.write("# COLUMNS = 1) Nucleotide Mutant\n"
                           "#           2) AminoAcid Mutant\n"
                           "#           3) Mutant counts before selection (maximum value from all primers)\n"
                           "#           4) Mutant counts after selection (maximum value from all primers)\n")
                fout.write("#\n")

                for mutlist in self.dmutcode['order']:
                    mut = mutlist[-1]
                    list_primer = list(self.dmutcode[nucmut].keys())
                    aamut = self.dmutcode[mut][list_primer[0]][0]
                    fout.write("{}\t{}\t".format(mut, aamut))
                    for condition in self.conditions:
                        count = max([int(x) for x in self.dmut_cond[mut][condition]])
                        fout.write("{}\t".format(count))
                    fout.write("\n")




    def read_config_file(self):
        """

        :return:
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)

        self.exe_dir    = config.get("Paths", "EXE_DIR")
        self.qqsub_exe  = config.get("Paths", "QSUB_EXE")

        self.use_qqsub   = eval(config.get("Parameters", "USE_QSUB"))

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
    parser.add_argument('-i', '--input_name', action="store", help='Root name of the input files')
    parser.add_argument('-p', '--primer_range', action="store", default='1-5', help='How many forward primers were used ex:1-5')
    parser.add_argument('-s', '--selection_conditions', action="store", help='Integer defining before (1) or after (2) selection or the two conditions. ex: 1,2')
    parser.add_argument('-d', '--working_directory', action="store", default='./', help='Working directory')
    parser.add_argument('-o', '--nameout', action="store", default='output', help='output_file')
    parser.add_argument('--mpi', action="store_true", default=False, help='triggers mpi running')
    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.config_file == default_configfile:
        os.system("cp {} {}".format(default_configfile,os.getcwd()))

    if args.input_name and args.primer_range and args.selection_conditions:
        """
        Here generate the code for generating the hhr file first
        """
        P = Process(args.input_name,
                      args.primer_range,
                      args.selection_conditions,
                      args.working_directory,
                      args.nameout,
                      args.mpi,
                      args.config_file,
                      verbose = args.verbose,
                      cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
