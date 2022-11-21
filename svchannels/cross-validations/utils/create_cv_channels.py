import argparse
import numpy as np
import zarr
import os
import json
import gzip
from time import time


def load_channels(channels_dir):
    samples_list = os.listdir(channels_dir)

    X = []
    y = []
    print('List of samples:')
    for s in samples_list:
        print(s)

        channels = os.path.join(channels_dir, s, 'channels', 'channels.zarr.zip')
        labels = os.path.join(channels_dir, s, 'labels', 'labels.json.gz')

        X_partial = zarr.load(channels)

        with gzip.GzipFile(labels, 'r') as fin:
            y_partial = json.loads(fin.read().decode('utf-8'))

        k_s = [s.split('/') for s in y_partial.keys()]
        y_array = [[s, k[0], k[1], k[2], k[3], '/'.join(k), v]
                   for k, v in zip(k_s, y_partial.values())]
        # print(y_array)
        X.extend(X_partial)
        y.extend(y_array)

    X = np.stack(X, axis=0)
    y = np.array(y)

    # select only chromosomes from chr1 to chr22
    chr_list = ['chr' + str(i) for i in np.arange(1, 23)]
    chr_list_y1 = y[:, 1]
    chr_list_y2 = y[:, 3]
    idx_chrlist = [i for i in np.arange(len(chr_list_y1)) if chr_list_y1[i] in chr_list and
                   chr_list_y2[i] in chr_list]
    X = X[idx_chrlist, ]
    y = y[idx_chrlist, ]

    print('----------------')
    print('shape of X is {}'.format(X.shape))
    print('shape of y is {}'.format(y.shape))

    return X, y


def main():

    random_seed = 42
    np.random.seed(random_seed)

    parser = argparse.ArgumentParser(description='Load channels from multiple samples from'
                                                 'a directory and save them in npz format')

    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default='svchan',
                        help="File with the channels (X) and labels (y) in npz format")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='multisample_svchan',
                        help="Output filename")
    args = parser.parse_args()

    t0 = time()

    X, y = load_channels(args.input)
    np.savez(args.output, X=X, y=y)

    print('Elapsed time = %f seconds' %
          (time() - t0))


if __name__ == '__main__':
    main()
