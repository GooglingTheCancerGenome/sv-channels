"""
Leave one chromosome and one sample out procedure
"""

import os
import argparse
import logging
import zarr
import gzip
import json
import tensorflow as tf
import numpy as np
from time import time
from sklearn.utils.class_weight import compute_class_weight
from tensorflow.keras.utils import to_categorical

from model_functions import create_model


def load_windows(win_file, lab_file):

    X = zarr.load(win_file)
    with gzip.GzipFile(lab_file, 'r') as fin:
        y = json.loads(fin.read().decode('utf-8'))
    return X, y


def train(args):

    randomState = 42
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}

    batch_size = args.batch_size
    max_epoch = args.epochs

    # chr_list = [str(i) for i in np.arange(1, 23)]

    windows = args.windows.split(',')
    labels = args.labels.split(',')
    samples = args.samples.split(',')

    X = []
    y = []
    win_pos = []
    samples_list = []

    for w, l, s in zip(windows, labels, samples):
        if s != args.sample_out:
            partial_X, partial_y = load_windows(w, l)
            X.extend(partial_X)
            y.extend(partial_y.values())
            win_pos.extend(partial_y.keys())
            # add sample name
            samples_list.extend([s] * len(partial_y))

    X = np.stack(X, axis=0)

    first_chrom = [w.split('/')[0] for w in win_pos]

    print('Running training leaving chromosome {} out'.format(args.test_chr))

    model_dir = os.path.dirname(args.model)
    model_base = os.path.basename(args.model)
    path_best_model = model_dir + '/' + args.test_chr + '.' + model_base

    chrom_idx = [i for i, k in enumerate(first_chrom) if k != args.test_chr]
    chrom_idx = np.asarray(chrom_idx)

    train_X = X[chrom_idx]
    y_nochrom = [y[i] for i in chrom_idx]

    y_nochrom = np.array([mapclasses[i] for i in y_nochrom])
    classes = np.array(np.unique(y_nochrom))
    y_lab = np.asarray(y_nochrom)

    class_weights = compute_class_weight('balanced', classes, y_lab)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    train_y = to_categorical(y_lab, num_classes=2)

    hparams = np.load(args.hparams)

    # cnn_init_learning_rate = 1e-4
    # cnn_regularization_rate = 1e-1
    # cnn_filters = 4
    # cnn_layers = 1
    # cnn_kernel_size = 7
    # cnn_fc_nodes = 4
    # dropout_rate = 0.2

    cnn_init_learning_rate = hparams[5]
    cnn_regularization_rate = hparams[6]
    cnn_filters = np.int32(hparams[0])
    cnn_layers = np.int32(hparams[1])
    cnn_kernel_size = np.int32(hparams[2])
    cnn_fc_nodes = np.int32(hparams[3])
    dropout_rate = hparams[4]

    outputdim = 2

    print()
    print('cnn_filters: ', cnn_filters)
    print('cnn_layers: ', cnn_layers)
    print('cnn_kernel_size: ', cnn_kernel_size)
    print('cnn_fc_nodes: ', cnn_fc_nodes)
    print('cnn_init_learning_rate: ', cnn_init_learning_rate)
    print('cnn_regularization_rate: ', cnn_regularization_rate)
    print('dropout_rate: ', dropout_rate)
    print()

    model = create_model(train_X, outputdim,
                         learning_rate=cnn_init_learning_rate,
                         regularization_rate=cnn_regularization_rate,
                         filters=cnn_filters,
                         layers=cnn_layers,
                         kernel_size=cnn_kernel_size,
                         fc_nodes=cnn_fc_nodes,
                         dropout_rate=dropout_rate)
    print(model.summary())

    history = model.fit(x=train_X, y=train_y,
                        epochs=max_epoch,
                        batch_size=batch_size,
                        shuffle=True,
                        validation_split=args.validation_split,
                        class_weight=class_weights,
                        verbose=0)

    accuracy = history.history['val_accuracy'][-1]
    print()
    print('Accuracy: {0:.2%}'.format(accuracy))

    model.save(path_best_model)

    tf.keras.backend.clear_session()


def main():
    parser = argparse.ArgumentParser(description='Optimize model')

    parser.add_argument('-w',
                        '--windows',
                        type=str,
                        default='sv_chan.zarr,sv_chan.zarr',
                        help="Comma separated list of training data")
    parser.add_argument('-lab',
                        '--labels',
                        type=str,
                        default='labels/labels.json.gz,labels/labels.json.gz',
                        help="Comma separated list of JSON.GZ file for labels")
    parser.add_argument('-sm',
                        '--samples',
                        type=str,
                        default='SAMPLE1,SAMPLE2',
                        help="Comma separated list of sample names")
    parser.add_argument('-so',
                        '--sample_out',
                        type=str,
                        default='SAMPLE',
                        help="Comma separated list of sample names")
    parser.add_argument('-l',
                        '--logfile',
                        default='optimize.log',
                        help='File in which to write logs.')
    parser.add_argument('-e',
                        '--epochs',
                        type=int,
                        default=50,
                        help="Number of epochs")
    parser.add_argument('-b',
                        '--batch_size',
                        type=int,
                        default=32,
                        help="Batch size")
    parser.add_argument('-test',
                        '--test_chr',
                        type=str,
                        default='chr1',
                        help="Chromosome used for testing")
    parser.add_argument('-val',
                        '--validation_split',
                        type=float,
                        default=0.3,
                        help="Proportion of training set used for validation")
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Type of SV")
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='best_model.keras',
                        help="Best model")
    parser.add_argument('-p',
                        '--hparams',
                        type=str,
                        default='hyperparams.npy',
                        help="File with hyperparameters")
    args = parser.parse_args()

    log_dir = os.path.dirname(args.logfile)
    log_base = os.path.basename(args.logfile)
    path_log = log_dir + '/' + args.test_chr + '.' + log_base

    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(format=log_format,
                        filename=path_log,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()

    train(args)

    logging.info('Elapsed time = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
