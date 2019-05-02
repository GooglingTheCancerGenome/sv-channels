# Imports
import gzip
import os
import argparse
import logging
import pickle
import errno

import numpy as np

from matplotlib import pyplot as plt

from collections import Counter, defaultdict

# Keras imports
from keras.utils.np_utils import to_categorical
from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution1D, Flatten, MaxPooling1D
from keras.optimizers import Adam
from keras.models import load_model
from keras import backend as K

from mcfly import modelgen, find_architecture

from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import label_binarize
from sklearn.utils import class_weight

# TensorFlow import
import tensorflow as tf

# Pandas import
import pandas as pd

from NA12878_mixed_data_CV import get_channel_labels

from model_functions import create_model, train_model, evaluate_model, unfold_win_id

HPC_MODE = True

sample_name = 'NA12878'
label_type = 'Mills2011'

# date = '260319'
date = '050419'

if HPC_MODE:

    datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari'
    channel_dir = datapath_prefix + '/Processed/Test/' + \
                  date + '/TestData_' + date + '/' + sample_name + '/TrainingData/'
else:
    channel_dir = '/Users/lsantuari/Documents/Processed/Test/test_060219'

chr_list = list(map(str, np.arange(1, 23)))
chr_list.append('X')

mapclasses = {'DEL': 0, 'noSV': 1}


def data(sample_name, label_type, suffix):

    data_output_file = os.path.join(channel_dir, '_'.join([sample_name, label_type, suffix]))

    with gzip.GzipFile(data_output_file + '.npz.gz', 'rb') as f:

        npzfiles = np.load(f)
        X = npzfiles['X']
        y = npzfiles['y']
        y_binary = npzfiles['y_binary']
        z = npzfiles['z']

    # print(X.shape)
    # print(y.shape)

    return X, y, y_binary, z


def save_data(sample_name, label_type, suffix, X, y, y_binary, z):
    data_output_file = os.path.join(channel_dir, '_'.join([sample_name, label_type, suffix]))
    np.savez(data_output_file, X=X, y=y, y_binary=y_binary, z=z)
    os.system('gzip -f ' + data_output_file + '.npz')


def create_training_and_test():

    def subset_data(X, y, y_binary, z, idx):
        return X[idx], y[idx], y_binary[idx], z[idx]

    training_data_output_file = os.path.join(channel_dir, '_'.join([sample_name, label_type, 'pairs_train'+'.npz.gz']))
    test_data_output_file = os.path.join(channel_dir, '_'.join([sample_name, label_type, 'pairs_test'+'.npz.gz']))

    if os.path.exists(training_data_output_file) and \
            os.path.exists(test_data_output_file):

        logging.info('Loading existing data for training and test set...')
        X_train, y_train, y_train_binary, z_train = data(sample_name, label_type, 'pairs_train')
        X_test, y_test, y_test_binary, z_test = data(sample_name, label_type, 'pairs_test')

    else:

        logging.info('Loading data...')
        X, y, y_binary, z = data(sample_name, label_type, suffix='pairs')
        logging.info('Data loaded')

        logging.info('X:%s, y:%s, y_binary:%s, z:%s' % (
            X.shape, y.shape, y_binary.shape, z.shape
        ))

        # select only split read positions on the same chromosome
        idx = np.where(np.array([coord_filter_by_chr(win_id) for win_id in z]))[0]
        X, y, y_binary, z = subset_data(X, y, y_binary, z, idx)

        test_chroms = ['1', '2', '3']
        train_chroms = [c for c in chr_list if c not in test_chroms]

        # select only split read positions on the test chromosomes
        idx = np.where(np.array([select_coord_by_chr(win_id, test_chroms) for win_id in z]))[0]
        X_test, y_test, y_test_binary, z_test = subset_data(X, y, y_binary, z, idx)
        save_data(sample_name, label_type, 'pairs_test', X_test, y_test, y_test_binary, z_test)

        # select only split read positions on the train chromosomes
        idx = np.where(np.array([select_coord_by_chr(win_id, train_chroms) for win_id in z]))[0]
        X_train, y_train, y_train_binary, z_train = subset_data(X, y, y_binary, z, idx)
        save_data(sample_name, label_type, 'pairs_train', X_train, y_train, y_train_binary, z_train)

    logging.info('Train => X:%s, y:%s, y_binary:%s, z:%s' % (
        X_train.shape, y_train.shape, y_train_binary.shape, z_train.shape
    ))
    logging.info('Test => X:%s, y:%s, y_binary:%s, z:%s' % (
        X_test.shape, y_test.shape, y_test_binary.shape, z_test.shape
    ))

    return (X_train, y_train, y_train_binary, z_train), (X_test, y_test, y_test_binary, z_test)


def coord_filter_by_chr(win_id):
    chr1, pos2, chr2, pos2 = unfold_win_id(win_id)

    if chr1 == chr2:
        return True
    else:
        return False


def select_coord_by_chr(win_id, chroms):
    chr1, pos2, chr2, pos2 = unfold_win_id(win_id)

    if chr1 == chr2 and chr1 in chroms:
        return True
    else:
        return False


def run_cv(index, output):

    filename, file_extension = os.path.splitext(output)

    def oversample(X, y):

        logging.info('Oversampling...')
        cnt_lab = Counter(y)

        max_v = max([v for k, v in cnt_lab.items()])

        data_balanced = []
        labels_balanced = []

        for l in cnt_lab.keys():
            # print(l)
            iw = np.where(y == l)
            # ii = iw[0][:min_v]
            ii = np.random.choice(a=iw[0], size=max_v, replace=True)
            data_balanced.extend(X[ii])
            labels_balanced.extend(y[ii])

        logging.info(Counter(labels_balanced))

        X = np.array(data_balanced)
        y = np.array(labels_balanced)

        return X, y

    def make_label_binary(y, mapclasses):

        y_num = np.array([mapclasses[c] for c in y], dtype='int')
        y_binary = to_categorical(y_num)
        return y_binary

    results = pd.DataFrame()
    metrics = dict()

    (X_train, y_train, y_train_binary, z_train), \
    (X_test, y_test, y_test_binary, z_test) = create_training_and_test()

    classlabels = mapclasses.keys()

    X_train, y_train = oversample(X_train, y_train)

    # split into train/validation sets
    xtrain, xval, ytrain, yval = train_test_split(X_train, y_train,
                                                  test_size=0.2, random_state=2,
                                                  stratify=y_train,
                                                  shuffle=True)

    # make binary labels
    ytrain_binary = make_label_binary(ytrain, mapclasses)
    yval_binary = make_label_binary(yval, mapclasses)

    logging.info('Test run:{0}'.format(str(int(index)+1)))

    # Create a new model
    model = create_model(xtrain, ytrain_binary)

    history, model = train_model(model, xtrain, ytrain_binary,
                                 xval, yval_binary)

    accuracy_history = history.history['acc']
    val_accuracy_history = history.history['val_acc']
    logging.info("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
        val_accuracy_history[-1]))

    # score_test = model.evaluate(X_hold_out_test, y_hold_out_test_binary, verbose=False)
    # print('Test loss and accuracy of best model: ' + str(score_test))

    intermediate_results, metrics[str(index + 1)] = evaluate_model(model, X_test, y_test,
                                                                   y_test_binary, z_test,
                                                                   results, index, filename,
                                                                   train_set_size=xtrain.shape[0],
                                                                   validation_set_size=xval.shape[0],
                                                                   mapclasses=mapclasses)
    results = results.append(intermediate_results)

    results.to_csv(os.path.join(date, filename + '_' + index + file_extension), sep='\t')


def main():

    parser = argparse.ArgumentParser(description='Test models on window pairs based on SR positions')
    parser.add_argument('-r', '--run', default=1,
                        help='Number of run.')
    parser.add_argument('-o', '--output', default='CV_results.csv',
                        help='File in which to write output.')
    parser.add_argument('-l', '--logfile', default='CV_results.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    filename, file_extension = os.path.splitext(args.logfile)
    logfilename = filename + file_extension

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    run_cv(index=args.run, output=args.output)


if __name__ == '__main__':
    main()
