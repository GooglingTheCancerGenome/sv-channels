#!/usr/bin/env python

"""
Description:
  MakeChannels.py is a command-line tool to generate different 'channels' from
  whole-genome sequencing (WGS) data for Deep Learning.

Usage:
  MakeChannels.py -h|--help
  MakeChannels.py -v|--version
  MakeChannels.py -t TYPE [-c CHR] [-o OUT_FILE] [-m MODE] BAM_FILE

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
  -c CHROM         Select a chromosome. [default: 17]
  -o OUT_FILE      Specify the output file. [default: all.npy.bz2]
  -m MODE          Run the tool using (s)imulated or (r)eal data. [default: r]
"""

from __future__ import print_function
from docopt import docopt
from pprint import pprint
from extras import Alignment


__authors__ = ["Arnold Kuzniar", "Luca Santuari"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.1.0"
__status__ = "alpha"


def main():
    args = docopt(__doc__, version=__version__)
    bam_file = args["BAM_FILE"]
    chrom = args["-c"]
    pprint(args)
    aln = Alignment(bam_file)
    print(aln.get_length(chrom))
    cov = aln.get_coverage(chrom)
    pprint(cov)
    print(time.asctime(time.localtime(time.time())))
    #clr = aln.get_clipped_reads(chrom)
    #pprint(clr)
    #print(time.asctime(time.localtime(time.time())))

if __name__ == "__main__":
    main()
