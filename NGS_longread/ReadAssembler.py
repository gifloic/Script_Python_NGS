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



class Assembler:
    def __init__(self, fforward , freverse, fileout, config_file, tandem_barcode = False, verbose = False, cleanup = True):
        self.fforward = fforward
        self.freverse = freverse
        self.root_fname = self.get_root_filename(fileout)
        self.config_file = config_file
        self.verbose  = verbose
        self.cleanup  = cleanup

        self.tandem_barcode = tandem_barcode

        self.d3to1, self.d1to3 = lt.get_dicaa()
        self.read_config_file()
        self.dmut = self.get_dic_expect_mutations()

        self.dreverse = self.parse_file(freverse, get_barcode = True)
        print("dreverse:",len(self.dreverse))
        """
        if self.tandem_barcode:
            self.dforward = self.parse_file(fforward, get_barcode = self.tandem_barcode)
            #self.dbarcode = self.compare_reverse_vs_forward()
            self.fdbarcode, self.rdbarcode = self.match_mutant_to_barcodes(tandem_barcode = self.tandem_barcode)
            self.d_obs_mutants = self.dump_barcode_dictionary()
        else:
        
        """
        self.dforward = self.parse_file(fforward, get_barcode = self.tandem_barcode)
        print("dforward:",len(self.dforward))
        self.fdbarcode, self.rdbarcode = self.match_mutant_to_barcodes()
        print("fdbarcode:",len(self.fdbarcode),"rdbarcode:",len(self.rdbarcode))
        self.fdobsmutants = self.count_mutants(self.root_fname,self.fdbarcode,'forward')
        self.rdobsmutants = self.count_mutants(self.root_fname,self.rdbarcode,'reverse')
        print("fdobsmutants:",len(self.fdobsmutants),"rdobsmutants:",len(self.rdobsmutants))
        self.analyze_master_mutant()

    def get_root_filename(self,fileout):
        """

        :return:
        """
        s = fileout.split('.')
        if len(s) == 1:
            return fileout
        else:
            return '.'.join(fileout.split('.')[:-1])

    def get_dic_expect_mutations(self):
        """

        :return:
        """
        return lt.parse_twist_devis_table(self.mutfile, self.dn2a, self.d3to1, start=None, stop=None, info_codon_index=None)

    def parse_file(self,file, get_barcode = False):
        d_read2barcodemut = {}
        with open(file) as fin:
            for line in fin:
                if line[0] == "#":
                    continue
                s = line.split()
                read_hexidx  = s[0]
                variant_code = s[2]

                in_mut_table = s[1]
                if get_barcode:
                    try:
                        barcode = s[3]
                    except IndexError:
                        barcode = None
                    d_read2barcodemut[read_hexidx] = [variant_code, in_mut_table, barcode]
                else:
                    d_read2barcodemut[read_hexidx] = [variant_code, in_mut_table]
        return d_read2barcodemut

    def compare_reverse_vs_forward(self):
        dbarcode = {}
        for ii,k in enumerate(self.dreverse):
            if k in self.dforward:
                forward = self.dforward[k]
                reverse = self.dreverse[k]
                if forward[-1] == reverse[-1]: # compare the two barcodes
                    barcode = forward[-1]
                    if not barcode: # Both forward and reverse are None
                        continue
                    if not barcode in dbarcode:
                        dbarcode[barcode] = {}
                    mutation = forward[0]
                    if not mutation in dbarcode[barcode]:
                        dbarcode[barcode][mutation] = 0
                    dbarcode[barcode][mutation] += 1
            else:
                if self.verbose:
                    print("READ {} NOT FOUND in forward file".format(k))
        return dbarcode

    def dump_barcode_dictionary(self):
        dmut = {}
        with open(self.fileout,"w") as fout:
            for bc in self.fdbarcode:
                list_mutant = []
                for mut in self.fdbarcode[bc]:
                    list_mutant.append([self.fdbarcode[bc][mut],mut])
                list_mutant.sort(reverse=True)
                if list_mutant[0][0] <2:
                    continue
                fout.write("{}\t".format(bc))
                for ii,m in enumerate(list_mutant):
                    if ii >= 4:
                        break
                    if m[1] == "W":
                        canonical = "W"
                    elif m[1] in self.dmut:
                        canonical = "O"
                    else:
                        canonical = "N"
                    fout.write("{} {} {} | ".format(m[0],canonical,m[1]))
                    if ii == 0:
                        if canonical in set(["W","O"]):
                            if not m[1] in dmut:
                                dmut[m[1]] = set()
                            dmut[m[1]].add(m[0])
                fout.write("\n")
        return dmut

    def match_mutant_to_barcodes(self,tandem_barcodes = False):
        fdbarcode = {}
        rdbarcode = {}
        for k in self.dreverse:
            if k in self.dforward:
                forward = self.dforward[k]
                reverse = self.dreverse[k]
                barcode = self.dreverse[k][-1]
                if not barcode: # Both forward and reverse are None
                    continue
                if tandem_barcodes:
                    if forward[-1] != barcode:
                        continue
                if not barcode in fdbarcode:
                    fdbarcode[barcode] = {}
                if not barcode in rdbarcode:
                    rdbarcode[barcode] = {}
                fmutation = forward[0]
                rmutation = reverse[0]
                if not fmutation in fdbarcode[barcode]:
                    fdbarcode[barcode][fmutation] = 0
                if not rmutation in rdbarcode[barcode]:
                    rdbarcode[barcode][rmutation] = 0
                rdbarcode[barcode][rmutation] += 1
                fdbarcode[barcode][fmutation] += 1
            else:
                if self.verbose:
                    print("READ {} NOT FOUND in forward file".format(k))
        return fdbarcode, rdbarcode

    def count_mutants(self,fileoutroot, dbarcode, code_forw_rev):
        """

        :return:
        """
        dmut = {}
        fout_name = '{}_{}.out'.format(fileoutroot,code_forw_rev)
        with open(fout_name,"w") as fout:
            for bc in dbarcode:
                list_mutant = []
                for mut in dbarcode[bc]:
                    list_mutant.append([dbarcode[bc][mut],mut])
                list_mutant.sort(reverse=True)
                if list_mutant[0][0] < 2:
                    continue
                fout.write("{}\t".format(bc))
                for ii,m in enumerate(list_mutant):
                    if ii >= 4:
                        break
                    if m[1] == "W":
                        canonical = "W"
                    elif m[1] in self.dmut:
                        canonical = "O"
                    else:
                        canonical = "N"
                    fout.write("{} {} {} | ".format(m[0],canonical,m[1]))
                    if ii == 0:
                        if canonical in set(["W","O"]):
                            if not m[1] in dmut:
                                dmut[m[1]] = set()
                            dmut[m[1]].add(m[0])
                fout.write("\n")
        return dmut

    def analyze_master_mutant(self):
        """

        :return:
        """
        ffout_name = self.root_fname + '_mutants_forward.out'
        rfout_name = self.root_fname + '_mutants_reverse.out'
        d_obs_mutants = [self.fdobsmutants, self.rdobsmutants]
        files = [ffout_name, rfout_name]
        for ii in range(2):
            with open(files[ii],"w") as fout:
                for mut_nuc_name in self.dmut['order']:
                    aa_mutant_name = self.dmut[mut_nuc_name]['full_mut_AA_name'].split("@@")
                    if mut_nuc_name in d_obs_mutants[ii]:
                        distrib_largest = list(d_obs_mutants[ii][mut_nuc_name])
                        distrib_largest.sort(reverse=True)
                        str_distrib_largest = '_'.join([str(x) for x in distrib_largest[:]])
                        mutant_number = len(d_obs_mutants[ii][mut_nuc_name])
                        fout.write("{}\t{}\t{}\t{}\n".format(mut_nuc_name, aa_mutant_name[1],mutant_number, str_distrib_largest))
                    else:
                        fout.write("{}\t{}\t0\t0\n".format(mut_nuc_name, aa_mutant_name[1]))

    def read_config_file(self):
        """

        :return:
        """
        config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
        config.read(self.config_file)

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
    parser.add_argument('-f', '--forward_parsed', action="store", help='input forward file')
    parser.add_argument('-r', '--reverse_parsed', action="store", help='input reverse file')
    parser.add_argument('-b', '--tandem_barcode', action="store_true", help='if barcode is present in forward and reverse')
    parser.add_argument('-o', '--root_outfilename', action="store", help='output barcode analysis file')

    parser.add_argument('-c', '--config_file', default=default_configfile, action="store", help='name of the configuration file')
    parser.add_argument('--verbose', action="store_true", default=False, help="print extra info")
    parser.add_argument('--cleanup', action="store_true", default=False, help="cleanup temporary files after execution")

    if len(sys.argv) == 1:
        print("type -h or --help for more help information")
        sys.exit(1)

    args = parser.parse_args()

    if args.forward_parsed and args.reverse_parsed:
        """
        Here generate the code for generating the hhr file first
        """
        A = Assembler(args.forward_parsed,
                      args.reverse_parsed,
                      args.root_outfilename,
                      args.config_file,
                      tandem_barcode = args.tandem_barcode,
                      verbose = args.verbose,
                      cleanup = args.cleanup)

    else:
        print("type -h or --help for more help information")
        sys.exit(1)


if __name__ == "__main__":


    Main()
