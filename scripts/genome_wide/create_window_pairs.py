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


def get_chr_list(chrom):
    # if sampleName in ['T1', 'G1', 'ART_INDEL_HET', 'ART_INDEL_HOM']:
    #     chrlist = ['17']
    # else:
    #     chrlist = list(map(str, range(1, 23)))
    #     chrlist.extend(['X'])
    #
    return [chrom]


def load_chr_array(channel_data_dir, chrlist):
    chr_array = dict()

    for c in chrlist:
        carray_file = os.path.join(channel_data_dir, 'chr_array',
                                   c + '_carray')
        logging.info('Loading file: {}'.format(carray_file))
        assert os.path.exists(carray_file), carray_file + ' not found'
        chr_array[c] = bcolz.open(rootdir=carray_file)
        logging.info('Array shape: {}'.format(chr_array[c].shape))

    return chr_array


def get_labels(label_file):
    # 'label_window_pairs_on_svcallset.json.gz')
    # 'label_window_pairs_on_split_read_positions.json.gz')
    with gzip.GzipFile(label_file, 'r') as fin:
        labels = json.loads(fin.read().decode('utf-8'))

    return labels


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


def get_window_by_id(win_id, chr_array, padding, win_hlen):
    chr1, pos1, chr2, pos2 = win_id.split('_')
    pos1 = int(pos1)
    pos2 = int(pos2)

    dask_arrays = list()
    dask_arrays.append(chr_array[chr1][pos1 - win_hlen:pos1 + win_hlen, :])
    dask_arrays.append(padding)
    dask_arrays.append(chr_array[chr2][pos2 - win_hlen:pos2 + win_hlen, :])
    return da.concatenate(dask_arrays, axis=0)


def get_windows(carrays_dir, outDir, chrom_list, win, label_file_path, mode, npz_mode, padding_len):

    chr_array = load_chr_array(carrays_dir, chrom_list)
    n_channels = chr_array[chrom_list[0]].shape[1]
    logging.info('{} channels'.format(n_channels))

    labels = get_labels(label_file_path)
    # labels = get_range(labels, 0, 10000)

    # if sampleName == 'T1':
    #     labels = {k: v for k, v in labels.items() if same_chr_in_winid(k)}

    logging.info('{} labels found: {}'.format(len(labels),
                                              Counter(labels.values())))

    if mode == 'training':

        labels_positive = {k: v for k, v in labels.items() if v == 'DEL'}
        labels_negative = {k: v for k, v in labels.items() if v == 'noDEL'}
        labels_negative = get_range(labels_negative, 0,
                                    len(labels_positive.keys()))
        labels_set = {'positive': labels_positive, 'negative': labels_negative}
        # labels_set = {'negative': labels_negative}

    elif mode == 'test':

        labels_set = {'test': labels}

    win_hlen = int(int(win) / 2)

    for labs_name, labs in labels_set.items():

        logging.info('Creating {}...'.format(labs_name))

        n_r = 10 ** 5
        # print(n_r)
        last_t = time()
        i = 1

        outfile = os.path.join(outDir, 'windows')

        bcolz_array = bcolz.carray(bcolz.zeros(shape=(0, int(win) * 2 +
                                                      padding_len, n_channels),
                                               dtype=np.float32),
                                   mode='w',
                                   rootdir=outfile + '_carray')

        padding = bcolz.zeros(shape=(padding_len, n_channels),
                              dtype=np.float32)

        if npz_mode:
            numpy_array = []

        logging.info('Creating dask_arrays_win1 and dask_arrays_win2...')
        for chr1, pos1, chr2, pos2, strand_info in map(unfold_win_id, labs.keys()):
            # logging.info("chr1={} pos1={} chr2={} pos2={}".format(
            #    chr1, pos1, chr2, pos2))
            if not i % n_r:
                # Record the current time
                now_t = time()
                # print(type(now_t))
                logging.info(
                    "%d window pairs processed (%f window pairs / s)" %
                    (i, n_r / (now_t - last_t)))
                last_t = time()

            dask_array = list()
            d = chr_array[chr1][pos1 - win_hlen:pos1 + win_hlen, :]
            dask_array.append(d)
            dask_array.append(padding)
            d = chr_array[chr2][pos2 - win_hlen:pos2 + win_hlen, :]
            dask_array.append(d)

            try:

                dask_array = np.concatenate(dask_array, axis=0)

                if npz_mode:
                    numpy_array.append(dask_array)

            except ValueError:

                print('{}:{}-{}:{}'.format(chr1, pos1, chr2, pos2))

                for d in dask_array:
                    print(d.shape)

            # print(type(dask_array))
            bcolz_array.append(dask_array)
            i += 1

        # bcolz_array.append(dask_array)
        bcolz_array.attrs['labels'] = labs
        bcolz_array.flush()
        logging.info(bcolz_array.shape)

        if npz_mode:
            numpy_array = np.stack(numpy_array, axis=0)
            logging.info('Numpy array shape: {}'.format(numpy_array.shape))
            np.savez(file=os.path.join(outDir, 'windows'),
                     data=numpy_array,
                     labels=labs)


def main():
    '''
    Main function for parsing the input arguments and calling the function to create windows
    :return: None
    '''


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
                        default='./labels/win200/DEL/split_reads',
                        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='windows.log',
                        help='File in which to write logs.')
    parser.add_argument('-w',
                        '--window',
                        type=str,
                        default=200,
                        help="Specify window size")
    parser.add_argument('-lb',
                        '--labels',
                        type=str,
                        default='./labels/win200/DEL/split_reads/labels.json.gz',
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

    cmd_name = 'windows'
    output_dir = os.path.join(args.outputpath, cmd_name)
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    # output_file = os.path.join(output_dir, args.out)

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

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    logging.info('Elapsed time create_windows = %f seconds' % (time() - t0))


if __name__ == '__main__':
    main()
