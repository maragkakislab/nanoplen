#!/usr/bin/env python

import argparse
# import statistics as st
import pysam

"""
Script to calculate the read length statistics of BAM/SAM file and all
reference Ids from inputfile (assuming the header is present).  The lengths
 for the reference will be NaN if no read is present in the inputfile. Users
can omit the lengths with NaN reads using '--no-zeros' option. Pysam package
 (tested on pysam 0.15.4) is required for the proper working of this script.
"""

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input filename with read information in SAM/BAM format")
parser.add_argument("-f", "--sam", action='store_true',
                    help="Input file format only if SAM; default BAM")
parser.add_argument("-r", "--ref-col-name", default="reference",
                    help="Reference output column-name, default:reference")
parser.add_argument("-c", "--len-col-name", default="count",
                    help="Read count output column name, e.g Length")
parser.add_argument("-n", "--opt-col-name", nargs="+",
                    help="Name of an optional column e.g. sample_name")
parser.add_argument("-v", "--opt-col-val", nargs="+",
                    help="Value for the optional column; same for all rows")
parser.add_argument("-s", "--col_delimiter", default="\t",
                    help="Delimiter to seperate the columns of the output-file, default : tab(\t)")
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
    header += args.opt_col_name
print(delim.join(header))

#Processing data for the reads (mapped)
for seq in bamfile:
    if seq.is_unmapped:
        continue
    reference = seq.reference_name
    # ref_len = seq.reference_length    # For checking reference length
    query_length = seq.query_length
    if reference not in reference_present:
        reference_present[reference] = []
    # Print output of reference ids, read length (and optional column)
    id_info = [reference, str(query_length)]
    if args.opt_col_name and args.opt_col_val:
        id_info += args.opt_col_val
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
                reference_present[ref] = []
                ids_not_present[ref] = []

for ref, count in ids_not_present.items():
    id_no_info = [ref,"NaN"]
    if args.opt_col_name and args.opt_col_val:
        id_no_info += args.opt_col_val
    print(delim.join(id_no_info))
