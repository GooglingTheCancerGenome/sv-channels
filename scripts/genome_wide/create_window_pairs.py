import argparse
import gzip
import itertools
import json
import logging
import os
from collections import Counter
from time import time

import bcolz
import numpy as np


def get_range(dictionary, begin, end):
    return dict(itertools.islice(dictionary.items(), begin, end))


def load_chr_array(channel_data_dir, chrlist):
    chr_array = dict()
    for c in chrlist:
        carray_file = os.path.join(
            channel_data_dir, 'chr_array', c + '_carray')
        logging.info("Loading file %s" % carray_file)
        chr_array[c] = bcolz.open(rootdir=carray_file)
        logging.info("Array shape: %s" % str(chr_array[c].shape))
        logging.info("Array data type: %s" % str(chr_array[c].dtype))
    return chr_array


def get_labels(label_file):
    with gzip.GzipFile(label_file, 'r') as fin:
        return json.loads(fin.read().decode('utf-8'))


def split_labels(labels):
    p = {}
    n = {}
    for k, v in labels.items():
        if v == 'DEL':
            p[k] = v
        elif v == 'noDEL':
            n[k] = v
    return p, n


def unfold_win_id(win_id):
    chr1, pos1, chr2, pos2, strand_info = win_id.split('_')
    pos1 = int(pos1)
    pos2 = int(pos2)
    return chr1, pos1, chr2, pos2, strand_info


def get_windows(carrays_dir, outDir, chrom_list, win, label_file_path, mode, npz_mode, padding_len):
    if win % 2 != 0:
        win += 1
    chr_array = load_chr_array(carrays_dir, chrom_list)
    n_channels = chr_array[chrom_list[0]].shape[1]
    logging.info("%d channels" % n_channels)
    labels = get_labels(label_file_path)
    logging.info("%d labels found: %s" %
                 (len(labels), str(Counter(labels.values()))))

    if mode == 'training':
        labels_positive = {k: v for k, v in labels.items() if v == 'DEL'}
        labels_negative = {k: v for k, v in labels.items() if v == 'noDEL'}
        labels_negative = get_range(
            labels_negative, 0, len(labels_positive.keys()))
        labels_set = {'positive': labels_positive, 'negative': labels_negative}
    elif mode == 'test':
        labels_set = {'test': labels}
    win_hlen = int(int(win) / 2)

    for labs_name, labs in labels_set.items():
        logging.info("Creating %s..." % str(labs_name))
        n_r = 10 ** 5
        last_t = time()
        i = 1
        padding = np.zeros(shape=(padding_len, n_channels), dtype=np.half)
        if npz_mode:
            numpy_array = []
        logging.info('Creating np.arrays win1 and win2...')
        for chr1, pos1, chr2, pos2, strand_info in map(unfold_win_id, labs.keys()):
            if not i % n_r:
                logging.info("%d window pairs processed (%f window pairs / s)" %
                             (i, n_r / (time() - last_t)))
                last_t = time()
            partial_array = list()
            d = chr_array[chr1][pos1 - win_hlen:pos1 + win_hlen, :]
            partial_array.append(d)
            partial_array.append(padding)
            d = chr_array[chr2][pos2 - win_hlen:pos2 + win_hlen, :]
            partial_array.append(d)

            try:
                full_array = np.concatenate(partial_array, axis=0)
                if npz_mode:
                    numpy_array.append(full_array)
            except ValueError:
                print('{}:{}-{}:{}'.format(chr1, pos1, chr2, pos2))
                for d in numpy_array:
                    print(d.shape)

        if npz_mode:
            numpy_array = np.stack(numpy_array, axis=0)
            logging.info("Numpy array shape: %s" % str(numpy_array.shape))
            for i in np.arange(numpy_array.shape[2]):
                logging.info("windows array: non-zero elements at index %d:%d" %
                             (i, np.argwhere(numpy_array[i, :] != 0).shape[0]))
            numpy_array= numpy_array.astype(np.half)
            np.savez(file=os.path.join(outDir, 'windows'),
                     data=numpy_array,
                     labels=labs)


def main():
    parser = argparse.ArgumentParser(
        description='Create windows from chromosome arrays')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chrlist',
                        type=str,
                        default='12,22',
                        help="List of chromosomes to consider")
    parser.add_argument('-ca',
                        '--carraydir',
                        type=str,
                        default='.',
                        help="chr_array directory")
    parser.add_argument('-p',
                        '--outputpath',
                        type=str,
                        default='./cnn/win200/split_reads/windows/DEL',
                        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='windows.log',
                        help='File in which to write logs.')
    parser.add_argument('-w',
                        '--window',
                        type=int,
                        default=200,
                        help="Specify window size")
    parser.add_argument('-lb',
                        '--labels',
                        type=str,
                        default='./cnn/win200/split_reads/windows/DEL/labels.json.gz',
                        help="Specify label file")
    parser.add_argument('-m',
                        '--mode',
                        type=str,
                        default='test',
                        help="training/test mode")
    parser.add_argument('-npz',
                        '--save_npz',
                        type=bool,
                        default=True,
                        help="save in npz format?")
    parser.add_argument('-pd',
                        '--padding',
                        type=int,
                        default=10,
                        help="Length of the padding in between windows")
    args = parser.parse_args()
    output_dir = args.outputpath
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    get_windows(carrays_dir=args.carraydir,
                outDir=output_dir,
                chrom_list=args.chrlist.split(','),
                win=args.window,
                label_file_path=args.labels,
                mode=args.mode,
                npz_mode=args.save_npz,
                padding_len=args.padding)
    logging.info('Elapsed time create_windows = %f seconds' % (time() - t0))


if __name__ == '__main__':
    main()
