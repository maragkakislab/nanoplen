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
# import warnings
# from pandas.core.common import SettingWithCopyWarning
# warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
matplotlib.use('pdf')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-i", "--ifile",
                    help="Input file from nanoplen with ids, mean-lengths etc")
parser.add_argument("-t", "--treatment", default="mean_length.alt",
                    help="Treatment library column name (default: %(default)s)")
parser.add_argument("-c", "--control", default="mean_length.control",
                    help="Control library column name (default: %(default)s)")
parser.add_argument("-l", "--lowerlimit", default=8, type=float,
                    help="Lower limit of X-axis (default: %(default)f)")
parser.add_argument("-u", "--upperlimit", default=13, type=float,
                    help="Upper limit of X-axis (default: %(default)f)")
# parser.add_argument("-o", "--ofile", help="output folder redirection")
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

control_log2 = str(args.control) + "_log2" #define the name of new column
df[control_log2] = np.log2(df[args.control])
treat_log2 = str(args.treatment) + "_log2" #define the name of new column
df[treat_log2] = np.log2(df[args.treatment])


df_melted_ = pd.melt(df, id_vars =['name'], value_vars =['mean_length.control','mean_length.alt'])
df_melted = df_melted.rename(columns={'variable': 'Sample', 'value': 'Length'})
df_melted["Length_log2"] = np.log2(df_melted["Length"])

# max_length = df_melted["Length"].max()
# min_length = df_melted["Length"].min()
mean_length = df_melted["Length"].mean()
# max_length = df_melted["Length"].max()
# min_length = df_melted["Length"].min()
mean_length = df_melted["Length_log2"].mean()


#
# if args.lowerlimit:
#     min_length=args.lowerlimit
# if args.upperlimit:
#     max_length=args.upperlimit

#Plotting Begins
with PdfPages(args.ofigure) as pdf:
    sns.displot(data=df_melted,
                x="Length_log2",
                hue="Sample",
                kind="kde")
    plt.suptitle('kernel density estimate', y=0.95, x=0.80)
    plt.xlabel('Length (log scale)', size = 20)
    plt.ylabel('Density', size = 20)
    plt.xlim(args.lowerlimit, args.upperlimit)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    pdf.savefig()
    plt.close()

    sns.displot(data=df_melted,
                x="Length",
                hue="Sample",
                kind="kde")
    plt.suptitle('kernel density estimate', y=0.95, x=0.80)
    plt.xlabel('Length', size = 20)
    plt.ylabel('Density', size = 20)
    # plt.xlim(args.lowerlimit, args.upperlimit)
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    pdf.savefig()
    plt.close()




    # sns.displot(data=df_melted,
    #             x="Length_log2",
    #             hue="Sample",
    #             kind="ecdf")
    # plt.suptitle('Cumulative Distribution Function', y=0.95, x=0.80)
    # plt.xlabel('Length (log scale)', size = 20)
    # plt.ylabel('Proportion', size = 20)
    # plt.xlim(args.lowerlimit, args.upperlimit)
    # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    # pdf.savefig()
    # plt.close()


    g = sns.JointGrid(data=df,
                      x=control_log2,
                      y=treat_log2,
                      )
    g.plot_joint(sns.histplot)
    g.plot_marginals(sns.boxplot)
    g.ax_joint.axline((5, 5),
                        (6, 6),
                        color='grey',
                        linestyle=':')
    g.set_axis_labels(args.control, args.treatment, size = 20)
    g.plot.xlabel('Length (log scale)', size = 20)
    plt.ylabel('Proportion', size = 20)
    pdf.savefig()
    plt.close()




    # g2 = sns.JointGrid(data=df_non_zero,
    #                    x=df_non_zero.columns[1],
    #                    y=df_non_zero.columns[0],
    #                    xlim=(0, max_avglen+100),
    #                    ylim=(0, max_avglen+100))
    # g2.plot_joint(sns.histplot)
    # g2.plot_marginals(sns.histplot)
    # g2.ax_joint.axline((0, 0), (1, 1), color='grey', linestyle=':')
    # g2.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
    # pdf.savefig()
    # plt.close()
    #
    # #plot the figure-log2 values -- scatter_density_boxplot
    # g3 = sns.JointGrid(data=df_non_zero_log2,
    #                    x=df_non_zero_log2.columns[1],
    #                    y=df_non_zero_log2.columns[0],
    #                    xlim=(min_avglenlog2-1, max_avglenlog2+1),
    #                    ylim=(min_avglenlog2-1, max_avglenlog2+1))
    # g3.plot_joint(sns.histplot)
    # g3.plot_marginals(sns.boxplot)
    # g3.ax_joint.axline((mean_avglenlog2, mean_avglenlog2), (mean_avglenlog2+1,
    # mean_avglenlog2+1), color='grey', linestyle=':')
    # g3.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
    # # plt.savefig(args.ofile + "scatter_density_boxplot_log2.pdf", dpi=300)
    # pdf.savefig()
    # plt.close()
    #
    # #plot the figure-log2 values -- scatter_density_boxplot
    # # Limited axis --> 8, 13
    # g5 = sns.JointGrid(data=df_non_zero_log2,
    #                    x=df_non_zero_log2.columns[1],
    #                    y=df_non_zero_log2.columns[0],
    #                    xlim=(args.lowerlimit, args.upperlimit),
    #                    ylim=(args.lowerlimit, args.upperlimit))
    # g5.plot_joint(sns.histplot)
    # g5.plot_marginals(sns.boxplot)
    # g5.ax_joint.axline((args.lowerlimit+1, args.lowerlimit+1),
    # (args.lowerlimit+2, args.lowerlimit+2), color='grey', linestyle=':')
    # g5.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
    # # plt.savefig(args.ofile + "scatter_density_boxplot_log2.pdf", dpi=300)
    # pdf.savefig()
    # plt.close()
    #
    # #plot the figure-log2 values -- scatter_density_boxplot --> NO Limits
    # g6 = sns.JointGrid(data=df_non_zero_log2,
    #                    x=df_non_zero_log2.columns[1],
    #                    y=df_non_zero_log2.columns[0],
    #                    # xlim=(8, 13),
    #                    # ylim=(8, 13)
    #                    )
    # g6.plot_joint(sns.histplot)
    # g6.plot_marginals(sns.boxplot)
    # g6.ax_joint.axline((mean_avglenlog2, mean_avglenlog2), (mean_avglenlog2+1,
    # mean_avglenlog2+1), color='grey', linestyle=':')
    # g6.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
    # pdf.savefig()
    # plt.close()
    #
    # #scatter_histogram_plot
    # g4 = sns.JointGrid(data=df_non_zero_log2,
    #                    x=df_non_zero_log2.columns[1],
    #                    y=df_non_zero_log2.columns[0],
    #                    xlim=(min_avglenlog2 -1, max_avglenlog2+1),
    #                    ylim=(min_avglenlog2-1, max_avglenlog2+1))
    # g4.plot_joint(sns.histplot)
    # g4.plot_marginals(sns.histplot)
    # g4.ax_joint.axline((mean_avglenlog2, mean_avglenlog2),
    # (mean_avglenlog2+1, mean_avglenlog2+1),
    # color='grey', linestyle=':')
    # g4.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
    # pdf.savefig()
    # plt.close()

################################################################################


# #Convert the datafrme to log2 scale (techniclaly the count/values column)
# new_log2_col = str(args.df_values) + "_log2" #define the name of new column
# df[new_log2_col] = np.log2(df[args.df_values])
#
# #statistics about our dataframe
# max_avglen = df[args.df_values].max()
# max_avglenlog2 = df[new_log2_col].max()
#
# min_avglen = df[args.df_values].min()
# min_avglenlog2 = df[new_log2_col].min()
#
# mean_avglen = df[args.df_values].mean()
# mean_avglenlog2 = df[new_log2_col].mean()

#convert dataframe to wideformat
# df_wide_format = df.pivot(index=args.df_index,
#                           columns=args.df_column,
#                           values=args.df_values)
# # Reorder the columns to keep condition before control (avoids unpredictable
# # alphanumeric ordering).
# df_wide_format = df_wide_format[[args.condition, args.control]]
#
# df_wide_format_log2 = df.pivot(index=args.df_index,
#                                columns=args.df_column,
#                                values=new_log2_col)
# # Reorder the columns to keep condition before control.
# df_wide_format_log2 = df_wide_format_log2[[args.condition, args.control]]
#
# #Drop all values with "zero" (NaN) values
# df_non_zero = df_wide_format.dropna()
# df_non_zero_log2 = df_wide_format_log2.dropna()

# #Rename columns
# df_non_zero_log2.rename(columns={df_non_zero_log2.columns[0]:df_non_zero_log2.columns[0]+'_log2',
#     df_non_zero_log2.columns[1]:df_non_zero_log2.columns[1]+'_log2'}, inplace=True)
#
# df_outfile = pd.concat([df_non_zero, df_non_zero_log2.reindex(df_non_zero.index)], axis=1)
# # df_outfile.to_csv(args.ofile + ".tab", sep='\t')
# df_outfile.to_csv(args.ofile, sep='\t')

# Drawing plots
# with PdfPages(args.ofigure) as pdf:
#     sns.displot(data=df_melted,
#                 x=args.df_values,
#                 hue="sample",
#                 kind="kde")
#     plt.suptitle('Normal Scale (kde)', y=0.95, x=0.80)
#     pdf.savefig()
#     plt.close()
#
#     sns.displot(data=df,
#                 x=args.df_values,
#                 hue="sample",
#                 kind="ecdf")
#     plt.suptitle('Normal Scale (ecdf)', y=0.95, x=0.80)
#     pdf.savefig()
#     plt.close()
#
#     sns.displot(data=df,
#                 x=new_log2_col,
#                 hue="sample",
#                 kind="kde")
#     plt.suptitle('Log2 Scale (kde)', y=0.95, x=0.80)
#     pdf.savefig()
#     plt.close()
#
#     sns.displot(data=df,
#                 x=new_log2_col,
#                 hue="sample",
#                 kind="ecdf")
#     plt.suptitle('Log2 Scale (ecdf)', y=0.95, x=0.80)
#     pdf.savefig()
#     plt.close()
#
#     # x-axis limited 8, 13
#     sns.displot(data=df,
#                 x=new_log2_col,
#                 hue="sample",
#                 kind="ecdf")
#     plt.suptitle('Log2 Scale (ecdf)', y=0.95, x=0.80)
#     plt.xlim(args.lowerlimit, args.upperlimit)
#     pdf.savefig()
#     plt.close()
#
#     # x-axis limited -- Custom
#     sns.displot(data=df,
#                 x=new_log2_col,
#                 hue="sample",
#                 kind="ecdf")
#     plt.suptitle('Log2 Scale (ecdf)', y=0.95, x=0.80)
#     plt.xlim(args.lowerlimit-0.5, args.upperlimit+0.5)
#     pdf.savefig()
#     plt.close()
#
#     # x-axis limited -- Custom
#     sns.displot(data=df,
#                 x=new_log2_col,
#                 hue="sample",
#                 kind="ecdf")
#     plt.suptitle('Log2 Scale (ecdf)', y=0.95, x=0.80)
#     plt.xlim(args.lowerlimit-0.75, args.upperlimit+0.75)
#     pdf.savefig()
#     plt.close()
#
#
#     g = sns.JointGrid(data=df_non_zero,
#                       x=df_non_zero.columns[1],
#                       y=df_non_zero.columns[0],
#                       xlim=(0, max_avglen+100),
#                       ylim=(0, max_avglen+100))
#     g.plot_joint(sns.histplot)
#     g.plot_marginals(sns.boxplot)
#     g.ax_joint.axline((0, 0), (1, 1), color='grey', linestyle=':')
#     g.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     pdf.savefig()
#     plt.close()
#
#     g2 = sns.JointGrid(data=df_non_zero,
#                        x=df_non_zero.columns[1],
#                        y=df_non_zero.columns[0],
#                        xlim=(0, max_avglen+100),
#                        ylim=(0, max_avglen+100))
#     g2.plot_joint(sns.histplot)
#     g2.plot_marginals(sns.histplot)
#     g2.ax_joint.axline((0, 0), (1, 1), color='grey', linestyle=':')
#     g2.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     pdf.savefig()
#     plt.close()
#
#     #plot the figure-log2 values -- scatter_density_boxplot
#     g3 = sns.JointGrid(data=df_non_zero_log2,
#                        x=df_non_zero_log2.columns[1],
#                        y=df_non_zero_log2.columns[0],
#                        xlim=(min_avglenlog2-1, max_avglenlog2+1),
#                        ylim=(min_avglenlog2-1, max_avglenlog2+1))
#     g3.plot_joint(sns.histplot)
#     g3.plot_marginals(sns.boxplot)
#     g3.ax_joint.axline((mean_avglenlog2, mean_avglenlog2), (mean_avglenlog2+1,
#     mean_avglenlog2+1), color='grey', linestyle=':')
#     g3.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     # plt.savefig(args.ofile + "scatter_density_boxplot_log2.pdf", dpi=300)
#     pdf.savefig()
#     plt.close()
#
#     #plot the figure-log2 values -- scatter_density_boxplot
#     # Limited axis --> 8, 13
#     g5 = sns.JointGrid(data=df_non_zero_log2,
#                        x=df_non_zero_log2.columns[1],
#                        y=df_non_zero_log2.columns[0],
#                        xlim=(args.lowerlimit, args.upperlimit),
#                        ylim=(args.lowerlimit, args.upperlimit))
#     g5.plot_joint(sns.histplot)
#     g5.plot_marginals(sns.boxplot)
#     g5.ax_joint.axline((args.lowerlimit+1, args.lowerlimit+1),
#     (args.lowerlimit+2, args.lowerlimit+2), color='grey', linestyle=':')
#     g5.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     # plt.savefig(args.ofile + "scatter_density_boxplot_log2.pdf", dpi=300)
#     pdf.savefig()
#     plt.close()
#
#     #plot the figure-log2 values -- scatter_density_boxplot --> NO Limits
#     g6 = sns.JointGrid(data=df_non_zero_log2,
#                        x=df_non_zero_log2.columns[1],
#                        y=df_non_zero_log2.columns[0],
#                        # xlim=(8, 13),
#                        # ylim=(8, 13)
#                        )
#     g6.plot_joint(sns.histplot)
#     g6.plot_marginals(sns.boxplot)
#     g6.ax_joint.axline((mean_avglenlog2, mean_avglenlog2), (mean_avglenlog2+1,
#     mean_avglenlog2+1), color='grey', linestyle=':')
#     g6.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     pdf.savefig()
#     plt.close()
#
#     #scatter_histogram_plot
#     g4 = sns.JointGrid(data=df_non_zero_log2,
#                        x=df_non_zero_log2.columns[1],
#                        y=df_non_zero_log2.columns[0],
#                        xlim=(min_avglenlog2 -1, max_avglenlog2+1),
#                        ylim=(min_avglenlog2-1, max_avglenlog2+1))
#     g4.plot_joint(sns.histplot)
#     g4.plot_marginals(sns.histplot)
#     g4.ax_joint.axline((mean_avglenlog2, mean_avglenlog2),
#     (mean_avglenlog2+1, mean_avglenlog2+1),
#     color='grey', linestyle=':')
#     g4.set_axis_labels(df_non_zero.columns[1], df_non_zero.columns[0])
#     pdf.savefig()
#     plt.close()
