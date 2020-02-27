#!/usr/bin/env python

"""
Description:
  sv-channels.py is a command-line tool to generate 'channels' from a BAM file.

Usage:
  sv-channels.py -h|--help
  sv-channels.py -v|--version
  sv-channels.py [-t TYPE] [-o PATH] BAM_FILE

Arguments:
  BAM_FILE             Input file in BAM format.

Options:
  -h, --help
  -v, --version
  -t, --type TYPE      Select a channel type: [default: cov]
                         cov = sequence coverage per position
  -o, --outpath PATH   Specify the output dir [default: .]
"""

from __future__ import print_function
from docopt import docopt
from pprint import pprint
from extras import Alignment

import os
import numpy as np


__authors__ = ['Arnold Kuzniar', 'Luca Santuari']
__license__ = 'Apache License, Version 2.0'
__version__ = '0.1.0'
__status__ = 'alpha'


def main():
    args = docopt(__doc__, version=__version__)
    # pprint(args)
    bamfile = args['BAM_FILE']
    outpath = args['--outpath']
    file_ext = '.npy'

    with Alignment(bamfile) as bam:
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        for seqid, cov in bam.get_seqcov():
            outfile = os.path.join(outpath, str(seqid) + file_ext)
            np.save(outfile, cov)


if __name__ == '__main__':
    main()
