#!/usr/bin/env python

"""
Pivot (convert table form long-to-wide format), log2-transform and create plots
as per given parameters/columns.
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
matplotlib.use('pdf')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input file from nanoplen with ids, mean-lengths etc")
parser.add_argument("-t", "--treatment", default="mean_length.alt",
                    help="Treatment library column name (default: %(default)s)")
parser.add_argument("-c", "--control", default="mean_length.control",
                    help="Control library column name (default: %(default)s)")
parser.add_argument("-l", "--lowerlimit", default=500, type=float,
                    help="Lower limit of X-axis (default: %(default)f)")
parser.add_argument("-u", "--upperlimit", default=5000, type=float,
                    help="Upper limit of X-axis (default: %(default)f)")
parser.add_argument("-f", "--ofigure", help="output figure name format in pdf")


args = parser.parse_args()

#Load longformat data to dataframe
if args.ifile == "-":
    ifile = sys.stdin
else:
    ifile = open(args.ifile, "r")

#table to dataframe
df = pd.read_csv(ifile, sep=r'\,|\t', engine='python') #regex patterns as the delimiter comma or tab

df = df[df != 0].dropna()

df_melted = pd.melt(df, id_vars =['name'], value_vars =[args.control,args.treatment])
df_melted = df_melted.rename(columns={'variable': 'Sample', 'value': 'Length'})

mean_length = df_melted["Length"].mean()

#Plotting Begins
with PdfPages(args.ofigure) as pdf:
    sns.displot(data=df_melted,
                x="Length",
                hue="Sample",
                kind="ecdf")
    plt.suptitle('ECDF', y=0.95, x=0.80)
    plt.xlabel('Length', size = 20)
    plt.ylabel('Density', size = 20)
    plt.xlim(args.lowerlimit, args.upperlimit)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    g = sns.JointGrid(data=df,
                      x=args.control,
                      y=args.treatment,
                      xlim=(args.lowerlimit, args.upperlimit),
                      ylim=(args.lowerlimit, args.upperlimit)
                      )
    g.plot_joint(sns.histplot)
    g.plot_marginals(sns.boxplot)
    g.ax_joint.axline((mean_length, mean_length),
                        (mean_length + 0.5, mean_length + 0.5),
                        color='red',
                        linestyle=':')
    g.set_axis_labels(args.control, args.treatment, size = 20)
    plt.xlabel('Length', size = 20)
    plt.ylabel('Proportion', size = 20)
    plt.tight_layout()
    pdf.savefig()
    plt.close()
