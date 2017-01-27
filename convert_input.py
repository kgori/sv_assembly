#!/usr/bin/env python

"""
Convert brass BEDPE format to TIGRA-SV (breakdancer) format

Input should be in BEDPE format, e.g.:

Chr1	12934	13024	Chr5	50399	51211	.	25	-	+

Output will be in Tigra/Breakdancer format, e.g.:

Chr1 12979	-	Chr5	50805	+	INV	.
"""

import pandas as pd
import os, sys

# Basic argument handling - argparse is overkill
if not len(sys.argv) == 3:
    sys.stderr.write("Usage: {} <infile.bedpe> <outfile.bedpe>\n".format(os.path.basename(sys.argv[0])))
    sys.exit()

infile = sys.argv[1]
if not os.path.exists(infile):
    sys.stderr.write("Input file not found {}\n".format(infile))
    sys.exit()

outfile = sys.argv[2]
dirname = os.path.dirname(outfile) or '.'
if not os.path.isdir(dirname):
  sys.stderr.write("Output directory not found {}\n".format(outfile))


data = pd.read_csv(infile, sep='\t', header=None)

data.columns = ["LeftName", "LeftStart", "LeftEnd",
                "RightName", "RightStart", "RightEnd",
                "Label", "Coverage",
                "LeftStrand",
                "RightStrand"]


tigra = pd.DataFrame(columns=["LeftName", "LeftPosition", "LeftStrand",
                              "RightName", "RightPosition", "RightStrand",
                              "Type", "Size"],
                     index=data.index).fillna(0)

tigra["LeftName"] = data["LeftName"]
tigra["LeftPosition"] = pd.Series(data[["LeftStart", "LeftEnd"]].mean(1), dtype=int)  # Guess the breakpoint is in the middle of the range
tigra["LeftStrand"] = data["LeftStrand"]
tigra["RightName"] = data["RightName"]
tigra["RightPosition"] = pd.Series(data[["RightStart", "RightEnd"]].mean(1), dtype=int)
tigra["RightStrand"] = data["RightStrand"]

# Specify all interchromosomal breaks
interchromosomal = data["LeftName"] != data["RightName"]
tigra.loc[interchromosomal, "Type"] = "CTX"
tigra.loc[interchromosomal, "Size"] = 1000  # Size is meaningless for these breaks, but will be skipped if < '-M' flag value (default 3)
tigra.loc[~interchromosomal, "Size"] = tigra[~interchromosomal]["RightPosition"] - tigra[~interchromosomal]["LeftPosition"]

duplications = (data["LeftStrand"] == "-") & (data["RightStrand"] == "+") & (~interchromosomal)
tigra.loc[duplications, "Type"] = "ITX"

inversions = ((data["LeftStrand"] == "+") & (data["RightStrand"] == "+") | (data["LeftStrand"] == "-") & (data["RightStrand"] == "-")) & (~interchromosomal)
tigra.loc[inversions, "Type"] = "INV"

tigra.loc[tigra["Type"]==0, "Type"] = "DEL"

tigra.to_csv(outfile, header=None, index=None, sep='\t')

