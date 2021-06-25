#!/usr/bin/env python

import sys
import argparse
import pysam

"""
Reads a SAM/BAM file and prints the reference sequence name and the read length
for each mapped read. If a header exists in the SAM/BAM then the references
without reads will have NaN for the length. An additional column will be
printed at the beginning if requested with --opt-col-name.
"""

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input SAM/BAM file")
parser.add_argument("-f", "--sam", action='store_true',
                    help="Use this Option if input file is in SAM format")
parser.add_argument("-r", "--ref-col-name", default="reference",
                    help="Name for column for reference ids, default:reference")
parser.add_argument("-c", "--len-col-name", default="length",
                    help="Name for column with the read length values, e.g length")
parser.add_argument("-n", "--opt-col-name", nargs="+",
                    help="Name of an optional column e.g. sample_name")
parser.add_argument("-v", "--opt-col-val", nargs="+",
                    help="Value for the optional column; same for all rows")
parser.add_argument("-s", "--col_delimiter", default="\t",
                    help="Column delimiter of output-file, default : tab(\t)")
parser.add_argument("-x", "--no-zeros", action='store_true',
                    help="Skip references with 0 reads")
args = parser.parse_args()

#Arguments about input file type, BAM format is default
ifiletype = "rb"
if args.sam:
    ifiletype = "r"

delim = args.col_delimiter

#Creating the dictionary for reads present
reference_present = {}

#Opening and processing the input file
bamfile = pysam.AlignmentFile(args.ifile, ifiletype)

#Print the header of file
header = [str(args.ref_col_name),
          str(args.len_col_name)
          ]

#option to add another column e.g. library name
if args.opt_col_name and args.opt_col_val:
    header = args.opt_col_name + header
elif args.opt_col_name  and not args.opt_col_val:
    print("Please provide input for both optional agruments", sys.exit(0))
elif args.opt_col_val and not args.opt_col_name:
    print("Please provide input for both optional agruments", sys.exit(0))
print(delim.join(header))

#Processing data for the reads (mapped)
for seq in bamfile:
    if seq.is_unmapped:
        continue
    reference = seq.reference_name
    query_length = seq.query_length
    if reference not in reference_present:
        reference_present[reference] = []
    # Print output of reference ids, read length (and optional column)
    id_info = [reference, str(query_length)]
    if args.opt_col_name and args.opt_col_val:
        id_info = args.opt_col_val + id_info
    print(delim.join(id_info))

# If transcripts with 0 counts are not explicitly excluded, then add them to
# the dictionary. Assumes that the SAM header exists. (At the end of outfile)
#Dictionary of ids with no reads
ids_not_present = {}

if not args.no_zeros:
    header = bamfile.header
    if 'SQ' in header:
        for elem in header['SQ']:
            ref = elem['SN']
            if ref not in reference_present:
                ids_not_present[ref] = []

for key in ids_not_present.keys():
    id_no_info = [key,"NaN"]
    if args.opt_col_name and args.opt_col_val:
        id_no_info = args.opt_col_val + id_no_info
    print(delim.join(id_no_info))
