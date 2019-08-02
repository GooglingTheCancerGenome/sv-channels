import dask.array as da
import h5py
import os, errno
import logging
import json
import argparse
from time import time
import gzip
import numpy as np
from collections import Counter
import itertools


def get_range(dictionary, begin, end):
    return dict(itertools.islice(dictionary.items(), begin, end + 1))


def create_dir(directory):
    '''
    Create a directory if it does not exist. Raises an exception if the directory exists.
    :param directory: directory to create
    :return: None
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def get_chr_list():
    # chrlist = list(map(str, range(1, 23)))
    # chrlist.extend(['X'])
    chrlist = ['17']

    return chrlist


def load_chr_array(channel_data_dir, sampleName):
    chrlist = get_chr_list()
    chr_array = dict()

    for c in chrlist:
        chrname = '/chr' + c
        hdf5_file = os.path.join(channel_data_dir, sampleName, 'chr_array', sampleName + '_' + c + '.hdf5')
        logging.info('Loading file: {}'.format(hdf5_file))
        assert os.path.exists(hdf5_file), hdf5_file + ' not found'
        f = h5py.File(hdf5_file)
        d = f[chrname]
        chr_array[c] = da.from_array(d, chunks=("auto", -1))

    return chr_array


def get_labels(channel_data_dir, sampleName, win):
    label_file = os.path.join(channel_data_dir, sampleName, 'labels_win' + str(win), 'labels.json.gz')

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
    chr1, pos1, chr2, pos2 = win_id.split('_')
    pos1 = int(pos1)
    pos2 = int(pos2)
    return chr1, pos1, chr2, pos2


def get_window_by_id(win_id, chr_array, padding, win_hlen):
    chr1, pos1, chr2, pos2 = win_id.split('_')
    pos1 = int(pos1)
    pos2 = int(pos2)

    dask_arrays = list()
    dask_arrays.append(chr_array[chr1][pos1 - win_hlen:pos1 + win_hlen, :])
    dask_arrays.append(padding)
    dask_arrays.append(chr_array[chr2][pos2 - win_hlen:pos2 + win_hlen, :])
    return da.concatenate(dask_arrays, axis=0)


def get_windows(sampleName, outDir, win, cmd_name, mode):

    def same_chr_in_winid(win_id):
        chr1, pos1, chr2, pos2 = win_id.split('_')
        return chr1 == chr2

    outfile_dir = os.path.join(outDir, sampleName, cmd_name)

    chr_array = load_chr_array(outDir, sampleName)
    n_channels = chr_array['17'].shape[1]
    labels = get_labels(outDir, sampleName, win)

    if sampleName == 'T1':
        labels = {k: v for k, v in labels.items() if same_chr_in_winid(k)}

    logging.info('{} labels found: {}'.format(len(labels), Counter(labels.values())))

    padding_len = 10
    win_hlen = int(int(win) / 2)

    dask_arrays_win1 = list()
    dask_arrays_win2 = list()

    for chr1, pos1, chr2, pos2 in map(unfold_win_id, labels.keys()):
        dask_arrays_win1.append(chr_array[chr1][pos1 - win_hlen:pos1 + win_hlen, :])
        dask_arrays_win2.append(chr_array[chr2][pos2 - win_hlen:pos2 + win_hlen, :])

    padding = da.zeros(shape=(len(labels.keys()), padding_len, n_channels), dtype=np.float32)
    dask_array = list()
    dask_array.append(da.stack(dask_arrays_win1, axis=0))
    dask_array.append(padding)
    dask_array.append(da.stack(dask_arrays_win2, axis=0))
    dask_array = da.concatenate(dask_array, axis=1)
    dask_array = dask_array.rechunk({0: 'auto', 1: -1, 2: -1})
    print(dask_array)

    outfile = os.path.join(outfile_dir, 'windows')

    da.to_hdf5(outfile + '.hdf5',
               {'/data': dask_array},
               compression='lzf')

    with gzip.GzipFile(outfile + '_labels.json.gz', 'wb') as fout:
        fout.write(json.dumps(labels).encode('utf-8'))


def main():
    '''
    Main function for parsing the input arguments and calling the function to create windows
    :return: None
    '''

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/'+
    #   'artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create windows from chromosome arrays')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-s', '--sample', type=str, default='T1',
                        help="Specify sample")
    parser.add_argument('-l', '--logfile', default='windows.log',
                        help='File in which to write logs.')
    parser.add_argument('-w', '--window', type=str, default=200,
                        help="Specify window size")
    parser.add_argument('-m', '--mode', type=str, default='training',
                        help="training/test mode")

    args = parser.parse_args()

    cmd_name = 'windows'
    output_dir = os.path.join(args.outputpath, args.sample, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, args.logfile)
    # output_file = os.path.join(output_dir, args.out)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    get_windows(sampleName=args.sample,
                outDir=args.outputpath,
                win=args.window,
                cmd_name=cmd_name,
                mode=args.mode
                )

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    logging.info('Elapsed time create_windows = %f seconds' % (time() - t0))


if __name__ == '__main__':
    main()
