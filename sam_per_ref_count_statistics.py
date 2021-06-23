#!/usr/bin/env python

import argparse
import statistics as st
import pysam

"""
Script to calculate the read counts, length statistics (like mean length,
median length, st.dev, variance, min length, max length,total reads and
mapped reads) of BAM/SAM file and all reference Ids from inputfile (assuming
the header is present). The counts for the reference will be zero if no read
is present in the inputfile. User can omit the counts with zero reads using
'--no-zeros' option in argparse. Pysam package (tested on pysam 0.15.4) is
required for the proper working of this script.
"""

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input filename with read information in SAM/BAM format")
parser.add_argument("-f", "--sam", action='store_true',
                    help="Input file format only if SAM; default BAM")
parser.add_argument("-r", "--ref-col-name", default="reference",
                    help="Reference output column-name, default:reference")
parser.add_argument("-c", "--cnt-col-name", default="count",
                    help="Read count output column name, default: count")
parser.add_argument("-n", "--opt-col-name", default="sample",
                    help="Name of an optional column e.g. sample_name")
parser.add_argument("-v", "--opt-col-val", default="sample_library-I",
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

#Creating the dictionary for reads
reference_counts = {}

#Initiating -- total and mapped reads
total_reads = 0
mapped_reads = 0


#Opening and processing the input file
bamfile = pysam.AlignmentFile(args.ifile, ifiletype)

#Processing data for the reads (mapped)
for seq in bamfile:
    total_reads += 1
    if seq.is_unmapped:
        continue
    mapped_reads += 1
    reference = seq.reference_name
    ref_len = seq.reference_length
    query_length = seq.query_length
    if reference not in reference_counts:
        reference_counts[reference] = []
    reference_counts[reference] += [query_length]

# If transcripts with 0 counts are not explicitly excluded, then add them to
# the dictionary. Assumes that the SAM header exists.
if not args.no_zeros:
    header = bamfile.header
    if 'SQ' in header:
        for elem in header['SQ']:
            ref = elem['SN']
            if ref not in reference_counts:
                reference_counts[ref] = []

delim = args.col_delimiter

#Print the header of file
header = [args.ref_col_name,
          args.cnt_col_name,
          "mean_length",
          "median_length",
          "stdev_length",
          "var_length",
          "max_length",
          "min_length",
          "total_reads",
          "mapped_reads"]

if args.opt_col_name and args.opt_col_val:
    header += [args.opt_col_name]
print(delim.join(header))

#Print the value per read
for ref, count in reference_counts.items():
    if len(count) == 1:
        row = [ref,
               str(len(count)),
               str(st.mean(count)),
               str(st.median(count)),
               str('inf'),
               str('inf'),
               str(max(count)),
               str(min(count)),
               str(total_reads),
               str(mapped_reads)]
    elif len(count) > 1:
        row = [ref,
               str(len(count)),
               str(st.mean(count)),
               str(st.median(count)),
               str(st.stdev(count)),
               str(st.variance(count)),
               str(max(count)),
               str(min(count)),
               str(total_reads),
               str(mapped_reads)]
    else:
        row = [ref,
               str(len(count)),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str('NaN'),
               str(total_reads),
               str(mapped_reads)]

    if args.opt_col_name and args.opt_col_val:
        row += [args.opt_col_val]
    print(delim.join(row))
