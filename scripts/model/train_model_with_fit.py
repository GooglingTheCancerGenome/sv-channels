import sys
import dask.array as da
import h5py
import os, errno
import numpy as np
import matplotlib.pyplot as plt
import json
import gzip
import pandas as pd
import logging
from time import time
import argparse

from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution1D, Lambda, \
    Convolution2D, Flatten, \
    Reshape, LSTM, Dropout, TimeDistributed, BatchNormalization, MaxPooling1D
from keras.regularizers import l2
from keras.optimizers import Adam
from keras.utils import to_categorical

from keras.callbacks import TensorBoard
from keras.models import load_model

from sklearn.model_selection import train_test_split
from collections import Counter

from model_functions import create_model, train_model, evaluate_model, create_dir


def get_channel_labels():
    # Fill labels for legend

    labels = list()
    labels.append("coverage")
    labels.append("discordant_reads_F")
    labels.append("discordant_reads_R")
    labels.append("mean_read_quality")

    labels.append("median_base_quality")
    labels.append("SNV_frequency")

    labels.append("#left_clipped_reads")
    labels.append("#right_clipped_reads")

    labels.append("#left split reads")
    labels.append("#right split reads")

    labels.append("#CIGAR_D_left_reads")
    labels.append("#CIGAR_D_right_reads")
    labels.append("#CIGAR_I_right_reads")

    labels.append("INV_before")
    labels.append("INV_after")
    labels.append("DUP_before")
    labels.append("DUP_after")
    labels.append("TRA_opposite")
    labels.append("TRA_same")

    for direction in ['Forward', 'Reverse']:
        for clipped in ['Left', 'Right', 'All']:
            for value in ['median']:
                # for value in ['median']:
                labels.append(direction + '_' + clipped + '_ClippedRead_distance_' + value)

    for clipped in ['L', 'R']:
        for value in ['median']:
            labels.append(clipped + '_SplitRead_distance_' + value)

    labels.append("Mappability")

    for nuc in ['A', 'T', 'C', 'G', 'N']:
        labels.append("One_hot_encoding_" + nuc)

    # for k, l in enumerate(labels):
    #     print(str(k) + ': ' + l)

    return labels


def plot_channels(outDir, X, z, l):
    title_plot = str(z) + '_' + str(l)
    print('Plotting %s' % title_plot)

    number_channels = X.shape[1]
    # print(number_channels)
    label = get_channel_labels()
    # print(len(label))

    fig = plt.figure(figsize=(6, 4))
    fig.suptitle(str(z) + ' ' + l, fontsize=20)

    for j in range(number_channels - 1, -1, -1):

        if sum(X[:, j]) != 0:
            X_win = (X[:, j] - min(X[:, j])) / max(X[:, j])
        else:
            X_win = X[:, j]

        Z = [x + j + 1 for x in X_win]
        plt.plot(Z, label=label[j], linewidth=0.9)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size': 5})
        plt.yticks(range(0, len(label) + 1, 1))
        plt.tick_params(axis='both', which='major', labelsize=5)
        plt.axvline(x=200, color='r', linewidth=0.05, alpha=0.5)
        plt.axvline(x=210, color='r', linewidth=0.05, alpha=0.5)

    plt.savefig(os.path.join(outDir, title_plot + '.png'),
                format='png', dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()


def get_hpc_flag():
    fileDir = os.path.dirname(os.path.abspath(__file__))
    parentDir = os.path.dirname(fileDir)
    newPath = os.path.join(parentDir, 'genome_wide')
    # print('appending path: ' + newPath)
    sys.path.append(newPath)
    # print(os.listdir(newPath))

    from functions import get_config_file

    config = get_config_file()
    HPC_MODE = config["DEFAULT"]["HPC_MODE"]

    return HPC_MODE


def data(channel_dir, sampleName):
    y = []
    numpy_array = []
    win_ids = []

    # class_dict = {'positive': 'DEL', 'negative': 'noDEL'}

    #for label_type in ['positive', 'negative']:

        # fn = os.path.join(channel_dir, 'windows', label_type + '.hdf5')
        # d = h5py.File(fn)
        #
        # fn = os.path.join(channel_dir, 'windows', label_type + '_labels.json.gz')
        #
        # with gzip.GzipFile(fn, 'r') as fin:
        #     labels = json.loads(fin.read().decode('utf-8'))

    for label_type in ['test']:

        carray_file = os.path.join(channel_dir,
                                   'windows', label_type + '_win200_carray')
        logging.info('Loading file: {}'.format(carray_file))
        assert os.path.exists(carray_file), carray_file + ' not found'
        X = bcolz.open(rootdir=carray_file)

        labels = X.attrs['labels']
        y.extend(labels.values())
        win_ids.extend(labels.keys())

        # if label_type == 'positive':
        #
        #     numpy_array.append(d['data'])
        #     y.extend(labels.values())
        #     win_ids.extend(labels.keys())
        #
        # elif label_type == 'negative':
        #
        #     rnd_idx = np.random.choice(d['data'].shape[0], size=numpy_array[0].shape[0], replace=True)
        #
        #     labs = list(d['labels'].item().values())
        #     labs_keys = list(d['labels'].item().keys())
        #
        #     # print(len(labs))
        #     numpy_array.append(d['data'][rnd_idx])
        #     y.extend(list(map(lambda i: labs[i], rnd_idx)))
        #     win_ids.extend(list(map(lambda i: labs_keys[i], rnd_idx)))

    #X = np.concatenate(numpy_array, axis=0)
    #X = X[:, :, np.array([0,1,2,7,8,26,27])]
    # X = np.delete(X,33,2)
    print(X.shape)
    print(Counter(y))

    my_classes = {'DEL': 0, 'noDEL': 1}
    y = np.array([my_classes[i] for i in y])
    win_ids = np.array(win_ids)

    # Shuffle arrays
    # new_indices = np.arange(X.shape[0])
    # np.random.shuffle(new_indices)
    # # print(new_indices)
    # X = X[new_indices]
    # y = y[new_indices]
    # win_ids = win_ids[new_indices]

    return X, y, win_ids


# def create_model(dim_length, dim_channels, class_number):
#
#     layers = 2
#     filters = [4] * layers
#     fc_hidden_nodes = 6
#     learning_rate = 4
#     regularization_rate = 1
#     kernel_size = 7
#     drp_out1 = 0
#     drp_out2 = 0
#
#     outputdim = class_number  # number of classes
#
#     weightinit = 'lecun_uniform'  # weight initialization
#
#     model = Sequential()
#     model.add(
#         BatchNormalization(
#             input_shape=(
#                 dim_length,
#                 dim_channels)))
#
#     for filter_number in filters:
#         # model.add(MaxPooling1D(pool_size=5, strides=None, padding='same'))
#
#         model.add(Convolution1D(filter_number, kernel_size=kernel_size, padding='same',
#                                 kernel_regularizer=l2(regularization_rate),
#                                 kernel_initializer=weightinit))
#         model.add(BatchNormalization())
#         model.add(Activation('relu'))
#
#     model.add(Flatten())
#     model.add(Dropout(drp_out1))
#     model.add(Dense(units=fc_hidden_nodes,
#                     kernel_regularizer=l2(regularization_rate),
#                     kernel_initializer=weightinit))  # Fully connected layer
#     model.add(Activation('relu'))  # Relu activation
#     model.add(Dropout(drp_out2))
#     model.add(Dense(units=outputdim, kernel_initializer=weightinit))
#     model.add(BatchNormalization())
#     model.add(Activation("softmax"))  # Final classification layer
#
#     model.compile(loss='categorical_crossentropy',
#                   optimizer=Adam(lr=learning_rate),
#                   metrics=['accuracy'])
#
#     return model


def train(channel_data_dir, sampleName):
    win_len = 200
    padding_len = 10

    dim = win_len * 2 + padding_len
    batch_size = 1
    epochs = 20

    # Parameters
    params = {'dim': dim,
              'batch_size': batch_size,
              'epochs': epochs,
              'n_classes': 2,
              'n_channels': 33,
              'shuffle': True}

    # Datasets
    X, y, win_ids = data(channel_data_dir, sampleName)

    plots_dir = os.path.join(channel_data_dir, 'plots_'+sampleName)
    create_dir(plots_dir)

    # for i, window in enumerate(X):
    #     plot_channels(plots_dir, window, y[i], win_ids[i])

    # split into train/validation sets
    X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                        test_size=0.3, random_state=2,
                                                        stratify=y,
                                                        shuffle=True)

    y_train_binary = to_categorical(y_train, num_classes=2)
    y_test_binary = to_categorical(y_test, num_classes=2)

    c = dict(Counter(y_train))
    total_labels = sum(c.values())
    class_weights = {k: v / total_labels for k, v in c.items()}

    model = create_model(X_train, y_train_binary)

    history, model = train_model(model, X_train, y_train_binary,
                                 X_test, y_test_binary)

    # # Design model
    # model = create_model(params['dim'], params['n_channels'], params['n_classes'])
    #
    # tbCallBack = TensorBoard(log_dir=os.path.join(channel_data_dir, 'Graph'),
    #                          histogram_freq=0,
    #                          write_graph=True,
    #                          write_images=True)
    #
    # # Train model on dataset
    # history = model.fit(X_train, y_train_binary,
    #                     validation_split=0.2,
    #                     batch_size=params['batch_size'],
    #                     epochs=params['epochs'],
    #                     shuffle=True,
    #                     # class_weight=class_weights,
    #                     verbose=1,
    #                     callbacks=[tbCallBack]
    #                     )
    #
    return model, history, X_train.shape[0], X_test.shape[0]


def train_and_test_model(sampleName_training, sampleName_test, outDir):

    HPC_MODE = get_hpc_flag()

    # channel_data_dir =
    channel_data_dir_training = \
        os.path.join('/Users/lsantuari/Documents/Processed/channel_maker_output',
                     sampleName_training) if not HPC_MODE else \
            os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data',
                         sampleName_training)

    model_fn = 'model_'+sampleName_training+'.hdf5'
    if os.path.exists(model_fn):
        model = load_model(model_fn)
    else:
        model, history, train_set_size, validation_set_size = train(channel_data_dir_training, sampleName_training)
        model.save(model_fn)

    results = pd.DataFrame()

    channel_data_dir_test = \
        os.path.join('/Users/lsantuari/Documents/Processed/channel_maker_output',
                     sampleName_test) if not HPC_MODE else \
            os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data',
                         sampleName_test)

    X_test, y_test, win_ids_test = data(channel_data_dir_test, sampleName_test)
    ytest_binary = to_categorical(y_test, num_classes=2)
    # print(win_ids_test[0])

    mapclasses = {'DEL': 0, 'noDEL': 1}

    intermediate_results, metrics = evaluate_model(model, X_test, y_test, ytest_binary, win_ids_test,
                                                   results, 1, 'results', mapclasses, outDir)

    results = results.append(intermediate_results)

    results.to_csv(os.path.join(outDir, 'results.csv'), sep='\t')

    # get_channel_labels()


def main():

    parser = argparse.ArgumentParser(description='Train and test model')
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-t', '--training_sample', type=str, default='NA12878',
                        help="Specify training sample")
    parser.add_argument('-x', '--test_sample', type=str, default='NA24385',
                        help="Specify training sample")
    parser.add_argument('-l', '--logfile', default='windows.log',
                        help='File in which to write logs.')
    parser.add_argument('-m', '--mode', type=str, default='training',
                        help="training/test mode")

    args = parser.parse_args()

    cmd_name = 'cnn'
    output_dir = os.path.join(args.outputpath, args.training_sample, cmd_name)
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

    train_and_test_model(sampleName_training=args.training_sample,
                         sampleName_test=args.test_sample,
                         outDir=output_dir
                         )

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    logging.info('Elapsed time training and testing = %f seconds' % (time() - t0))


if __name__ == '__main__':
    main()
