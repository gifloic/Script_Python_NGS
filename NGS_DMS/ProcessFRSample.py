#!/usr/bin/env python3
"""
This script combines the Parser and the Assembler to process forward and backward reads and returns a mutation analysis file
"""


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
    def __init__(self, input_name, primer_index, condition, working_dir, ncpu, mpi, config_file, verbose = False, cleanup = False):
        self.input_name = input_name
        self.condition = condition
        self.working_dir = working_dir
        self.ncpu = ncpu
        self.mpi = mpi
        self.config_file = config_file
        self.verbose  = verbose
        self.cleanup  = cleanup

        self.read_config_file()
        self.run_process()

    def run_process(self):

        os.chdir(self.working_dir)
        # Standard case
        process_path = '{}_{}'.format(self.input_name,self.condition)
        if not os.path.isdir(process_path):
            os.mkdir(process_path)

        os.chdir('ini_files/')

        f1 = '{}_1'.format(process_path) # forward file
        f2 = '{}_2'.format(process_path) # reverse file
        for fq in [f1,f2]:
            print('>>> Uncompresssing {} ...'.format(fq))
            if not os.path.isfile(fq+'.fastq'):
                subprocess.call(['gunzip','-k',fq+'.fastq.gz'])

        os.chdir('../{}'.format(process_path))
        for fq in [f1,f2]:
            if not os.path.isfile('{}.fastq'.format(fq)):
                os.symlink('../ini_files/{}.fastq'.format(fq),'{}.fastq'.format(fq))

        if os.path.basename(self.config_file) == "NGS_lread.ini":
            flag_config = ""
        else:
            flag_config = " -c {} ".format(self.config_file)

        cmd1 = '{}/ReadParser.py ' \
               '-f {}.fastq -o {}.out -e forward {} --mpi'.format( self.exe_dir,
                                                                        f1, f1, flag_config)
        cmd2 = '{}/ReadParser.py ' \
               '-f {}.fastq -o {}.out -e reverse {} --mpi'.format(self.exe_dir,
                                                                          f2, f2, flag_config)

        if self.use_qqsub:
            cmd1 = """{} -s "{}" -n 1 -p {} -t rp_cmd1 -w """.format(self.qqsub_exe, cmd1, self.ncpu)
            cmd2 = """{} -s "{}" -n 1 -p {} -t rp_cmd2 -w """.format(self.qqsub_exe, cmd2, self.ncpu)
        elif self.mpi:
            cmd1 = """mpirun -np {} python3 {}""".format(cmd1, self.ncpu)
            cmd2 = """mpirun -np {} python3 {}""".format(cmd2, self.ncpu)
        else:
            cmd1 = """python3 {}""".format(cmd1)
            cmd2 = """python3 {}""".format(cmd2)

        print(cmd1)
        print('>>> Parsing {} ...'.format(f1))
        os.system(cmd1)
        #subprocess.call(cmd1.split())

        print(cmd2)
        print('>>> Parsing {} ...'.format(f2))
        os.system(cmd2)
        #subprocess.call(cmd2.split())

        cmd3 = '{}/ReadAssembler.py -f {}.out -r {}.out -o {}_{}_assembled {}'.format(self.exe_dir, f1, f2,
                                                                                              self.input_name,
                                                                                              self.primer_index,
                                                                                              flag_config)
        if self.use_qqsub:
            cmd3 = """{} -s "{}" -n 1 -p 1 -t rp_cmd3 -w""".format(self.qqsub_exe, cmd3)
        else:
            cmd3 = 'python3 {}'.format(cmd3)

        print('>>> Assembling forward and reverse into {}_{}_assembled ...'.format(self.input_name,self.primer_index))
        #subprocess.call(cmd3.split())
        os.system(cmd3)
        if self.cleanup:
            for fq in [f1,f2]:
                print('>>> Removing the uncompressed fastq files')
                os.remove('../ini_files/{}.fastq'.format(fq))


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
                                  "NGS_DMS.ini")
    except:
        default_configfile = os.path.join(os.path.dirname(os.path.abspath(os.path.realpath(os.path.abspath(__file__)))),
                                  "NGS_DMS.ini")


    usage = "Script to parse the reads from a NGS DMS experiment.\n" + \
            "output is a dictionary encoded as a hdf5 file"
    parser = argparse.ArgumentParser(usage)
    parser.add_argument('-i', '--input_name', action="store", help='Root name of the input files')
    parser.add_argument('-s', '--selection_condition', action="store", help='Integer defining before (1) or after (2) selection')
    parser.add_argument('-d', '--working_directory', action="store", default='./', help='Working directory')
    parser.add_argument('-n', '--ncpu', action="store", default=40, help='Number of cpus')
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

    if args.input_name and args.primer_index and args.selection_condition:
        """
        Here generate the code for generating the hhr file first
        """
        P = Process(args.input_name,
                      args.primer_index,
                      args.selection_condition,
                      args.working_directory,
                      args.ncpu,
                      args.mpi,
                      args.config_file,
                      verbose = args.verbose,
                      cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
