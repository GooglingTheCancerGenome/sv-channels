#!/usr/bin/env python

"""
Description:
  MakeChannels.py is a command-line tool to generate 'channels' from a BAM file.

Usage:
  MakeChannels.py -h|--help
  MakeChannels.py -v|--version
  MakeChannels.py [-t TYPE] [-c CHR] [-o OUT_FILE] BAM

Arguments:
  BAM_FILE         Input file in BAM format.

Options:
  -h, --help
  -v, --version
  -t TYPE          Select a channel type: [default: all]
                     all = includes all the types below
                     cr = clipped reads
                     crp = clipped read positions
                     crd = clipped read distance
                     srd = split read distance
                     cov = per-base coverage
  -c CHROM         Select a chromosome. [default: all]
  -o OUTFILE       Specify the output file. [default: all.npz]
"""

from __future__ import print_function
from docopt import docopt
from pprint import pprint
from extras import Alignment

import numpy as np


__authors__ = ["Arnold Kuzniar", "Luca Santuari"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.1.0"
__status__ = "alpha"


def main():
    args = docopt(__doc__, version=__version__)
    bamfile = args["BAM"]
    outfile = args["-o"]
    chrom = args["-c"]
    pprint(args)
    aln = Alignment(bamfile)
    for chr, cov in aln.get_coverage():
        outfile = str(chr) + '.npy'
        np.save(outfile, cov)
    #clr = aln.get_clipped_reads(chrom)
    #pprint(clr)
    #print(time.asctime(time.localtime(time.time())))

if __name__ == "__main__":
    main()
