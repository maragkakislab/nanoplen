#!/usr/bin/env python


"""
Reads a nanopolish file and prints the reference sequence id and the poly-A
length for each "PASS" mapped read. An additional column will be
printed at the beginning if requested with '--opt-col-name'.
"""

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--ifile",
                    help="Input file from Nanopolish with polyA information")
parser.add_argument("-r", "--ref_col_name", default="reference",
                    help="Name for column for reference ids, default:reference")
parser.add_argument("-c", "--polyAlen_col_name", default="length",
                    help="Name for read length values column, e.g length")
parser.add_argument("-n", "--opt_col_name", nargs=1,
                    help="Name of an optional column e.g. sample")
parser.add_argument("-v", "--opt_col_val", nargs=1,
                    help="Value for respective optional column e.g. lib. name")
parser.add_argument("-t", "--index-type", default='number',
                    help="Option to drop columns either by name or number, default: name")
parser.add_argument("-l", "--drop_cols", nargs='+', default=[0, 2, 3, 4, 5, 6, 7, 9],
                    help="Name or number of column(s) to drop from input file")
parser.add_argument("-d", "--delimiter", default="\t",
                    help="Column delimiter of output-file, default : tab(\t)")
args = parser.parse_args()

if args.ifile == "-":
    ifile = sys.stdin
else:
    ifile = open(args.ifile, "r")

#Load data to dataframe
df = pd.read_csv(ifile, sep=args.delimiter)

#select only PASS values (valid values as per literature)
df = df.drop(df[df['qc_tag'] != "PASS"].index)

# Dropping uncessary columns
if args.index_type == "name":
    df.drop(args.drop_cols, axis=1, inplace=True)
elif args.index_type == "number":
    index_list = list(map(int, args.drop_cols))
    df.drop(df.columns[index_list], axis=1, inplace=True)

# Check for polyA length greater than value '1'
df.drop(df[df['polya_length'] < 1].index, inplace=True)

# Option to rename the column names
if args.ref_col_name:
    df = df.rename(columns={"contig": args.ref_col_name})
if args.polyAlen_col_name:
    df = df.rename(columns={"polya_length": args.polyAlen_col_name})


if args.opt_col_name and args.opt_col_val:
    df[args.opt_col_name[0]] = args.opt_col_val[0]
    df = df[[df.columns[2], df.columns[0], df.columns[1]]]
else:
    print("Please provide input for both optional agruments")
    print(sys.exit(0))


#Print the output to std-out
df.to_csv(sys.stdout, sep=args.delimiter, index=False)
