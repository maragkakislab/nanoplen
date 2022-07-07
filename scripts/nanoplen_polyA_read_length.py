#!/usr/bin/env python

"""
Reads the nanopolish output file, prints the reference sequence-id and the
poly-A length for each "PASS" mapped read and polyA length greaterthan equalto
1. An additional column will be printed at the beginning via the  argument
'--opt-col-name', and values in its rows from argument '--opt-col-val'.
"""

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input file from Nanopolish with polyA information")
parser.add_argument("-n", "--opt-col-name", required=False,
                    help="Name of optional column)")
parser.add_argument("-v", "--opt-col-val", required=False,
                    help="Values for the optional column")
parser.add_argument("-d", "--delimiter", default="\t",
                    help="Column delimiter of output-file, default : TAB")
args = parser.parse_args()

if args.ifile == "-":
    ifile = sys.stdin
else:
    ifile = open(args.ifile, "r")

#Load data to dataframe
df = pd.read_csv(ifile,
                 sep=args.delimiter,
                 usecols=['contig', 'polya_length', 'qc_tag'])

#Select only 'PASS' values (valid values as per literature) and length >= 1 bp.
df = df[((df["qc_tag"] == "PASS") & (df["polya_length"] >= 1))]

#Drop extra column with QC tag
df = df.drop('qc_tag', 1)

# Add extra column (optional) provided by user
if args.opt_col_name and args.opt_col_val:
    df.insert(0, args.opt_col_name, args.opt_col_val)
elif not args.opt_col_name and not args.opt_col_val:
    pass
else:
    print("Please provide input for for both or none optional agruments")
    print(sys.exit(1))

#Print the output to std-out
df.to_csv(sys.stdout, sep=args.delimiter, index=False)
