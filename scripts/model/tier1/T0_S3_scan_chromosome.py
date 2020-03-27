import argparse
import errno
import glob
import gzip
import json
import logging
import os
from time import time

import bcolz
import numpy as np
from keras.models import load_model


def load_bcolz_array(channels_dir, chrom):

    carray_fname = glob.glob(
        os.path.join(channels_dir, 'chr_array', '*_' + chrom + '_carray'))
    assert len(
        carray_fname) == 1, 'Not a single carray folder found:' + carray_fname
    carray_fname = carray_fname[0]
    assert os.path.exists(carray_fname), carray_fname + ' not found'
    chr_array = bcolz.open(rootdir=carray_fname)
    # logging.info('Array shape: {}'.format(chr_array[c].shape))

    return chr_array


def scan_chromosome(args):

    win = args.window
    half_win = int(args.window / 2)
    chrom = args.chr

    model = load_model(args.model)

    bc_array = load_bcolz_array(args.inputdir, chrom)

    step = 10**6
    batch_size = 10**5

    chr_len = bc_array.shape[0] // win * win

    dim1 = (step // win) * (chr_len // step) + (chr_len % step) // win
    dim2 = int(model.outputs[0].shape.dims[1])

    res_array = np.empty(shape=(dim1, dim2))

    for i in np.arange(0, chr_len // step + 1):

        vstart = i * step + args.shift
        vend = min((i + 1) * step + args.shift, chr_len)
        d0 = step // win

        print('Scanning chr {} from {} to {}'.format(chrom, vstart, vend))

        split = (chr_len % step) // win if i == (chr_len //
                                                 step) else step // win

        B = np.array(np.split(bc_array[vstart:vend, :], split))

        # print(B.shape)
        probs = model.predict_proba(B, batch_size=batch_size, verbose=False)
        # print(probs)
        res_array[i * d0:(i + 1) * d0] = probs

    predicted = res_array.argmax(axis=1)

    center_pos = np.arange(half_win + args.shift, len(predicted) * win, win)

    mapclasses = {'DEL_start': 0, 'DEL_end': 1, 'noSV': 2}
    bp_start = center_pos[np.where(
        np.array(predicted) == mapclasses['DEL_start'])]
    bp_end = center_pos[np.where(np.array(predicted) == mapclasses['DEL_end'])]

    # Write
    np.savez_compressed(file=args.output,
                        start=bp_start,
                        end=bp_end,
                        probs=res_array)


def main():

    parser = argparse.ArgumentParser(description='Train model',
                                     usage='''T0_S2_training.py [<args>]
        ''')
    parser.add_argument(
        '-inputdir',
        type=str,
        default='/Users/lsantuari/Documents/Processed/channel_maker_output/T1',
        help="Specify channel directory")
    parser.add_argument('-window',
                        type=int,
                        default=200,
                        help="Specify window size")
    parser.add_argument('-shift', type=int, default=0, help="Specify shift")
    parser.add_argument('-chr',
                        type=str,
                        default='17',
                        help="Specify chromosome")
    parser.add_argument('-model',
                        type=str,
                        default='model.hdf5',
                        help="Specify model")
    parser.add_argument('-output',
                        type=str,
                        default='17_predictions.npz',
                        help="Specify output")

    args = parser.parse_args()

    scan_chromosome(args)


if __name__ == '__main__':
    main()
