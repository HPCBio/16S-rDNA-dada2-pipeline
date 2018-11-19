#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:58:57 2013
"""
__author__ = "Patricio Jeraldo"
__copyright__ = "Copyright 2014-2016, Mayo Foundation for Medical Education and Research "
__credits__ = ["Patricio Jeraldo", "Nicholas Chia"]
__license__ = "GPL"
__version__ = "2.0.3.3"
__maintainer__ = "Patricio Jeraldo"
__email__ = "jeraldo.patricio@mayo.edu"
__status__ = "Production"

import argparse
from Bio import SeqIO

parser= argparse.ArgumentParser(description="Convert Stockholm-formatted alignment files to fasta files. Only for small (< 100000 reads) alignments.")

parser.add_argument('in_stk', help="Input Stockholm file.")
parser.add_argument('out_fasta', help="Output fasta file.")

args= parser.parse_args()

SeqIO.convert(args.in_stk, "stockholm", args.out_fasta, "fasta")
