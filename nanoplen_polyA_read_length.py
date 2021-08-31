#!/usr/bin/env python

"""
Reads the nanopolish output file and prints the reference sequence-id and the
poly-A length for each "PASS" mapped read and polyA length greaterthan equalto
 1. An additional column will be printed at the beginning via the  agrgument
 '--opt-col-name', with values in rows from argument '--opt-col-val'.
"""

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input file from Nanopolish with polyA information")
parser.add_argument("-r", "--ref-col-name", default="reference",
                    help="Name for column for reference ids (default: %(default)s)")
parser.add_argument("-c", "--polya-len-col-name", default="polya_length",
                    help="Name for read length values column (default: %(default)s)")
parser.add_argument("-n", "--opt-col-name", required=False,
                    help="Name of an optional column (default: %(default)s)")
parser.add_argument("-v", "--opt-col-val", required=False,
                    help="Value for respective optional column (default: %(default)s)")
parser.add_argument("-d", "--delimiter", default="\t",
                    help="Column delimiter of output-file, default : TAB")
args = parser.parse_args()

if args.ifile == "-":
    ifile = sys.stdin
else:
    ifile = open(args.ifile, "r")

#Load data to dataframe
df = pd.read_csv(ifile, sep=args.delimiter, usecols=['contig', 'polya_length', 'qc_tag'])

#select only PASS values (valid values as per literature) and length >= 1 bp.
df = df[((df["qc_tag"] == "PASS") & (df["polya_length"] >= 1))]

# Option to rename the column names
if args.ref_col_name:
    df = df.rename(columns={"contig": args.ref_col_name})
if args.polya_len_col_name:
    df = df.rename(columns={"polya_length": args.polya_len_col_name})

# Add an extra column that is provided by user, column-name and its row-values
if args.opt_col_name and args.opt_col_val:
    df.insert(0, args.opt_col_name, args.opt_col_val)

#Drop extra column with QC tag
df = df.drop('qc_tag', 1)

#Print the output to std-out
df.to_csv(sys.stdout, sep=args.delimiter, index=False)
