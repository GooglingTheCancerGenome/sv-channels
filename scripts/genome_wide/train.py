from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import gzip
import json
import logging
import os
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.utils.class_weight import compute_class_weight

from tensorflow.keras.callbacks import (EarlyStopping, ModelCheckpoint,
                                        TensorBoard)
from tensorflow.keras.layers import (Activation, BatchNormalization,
                                     Convolution1D, Dense, Dropout, Flatten,
                                     Lambda, Reshape, TimeDistributed)
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras.utils import to_categorical

from model_functions import (  # create_model_with_mcfly, train_model_with_mcfly
    evaluate_model, get_data)


def get_labels(channel_data_dir, win):
    label_file = os.path.join(channel_data_dir, 'labels_win' + str(win),
                              'labels.json.gz')

    with gzip.GzipFile(label_file, 'r') as fin:
        labels = json.loads(fin.read().decode('utf-8'))

    return labels


def train_and_test_data(sampleName, npz_mode, svtype):
    # Datasets
    X, y, win_ids = get_data(sampleName, npz_mode, svtype)

    X = np.array(X)
    y = np.array(y)

    # split into train/validation sets
    X_train, X_test, y_train, y_test, win_ids_train, win_ids_test = train_test_split(
        X, y, win_ids, test_size=0.3, random_state=2, stratify=y, shuffle=True)

    return X_train, X_test, y_train, y_test, win_ids_train, win_ids_test


def create_model(dim_length, dim_channels, outputdim):

    weightinit = 'lecun_uniform'  # weight initialization

    learning_rate = 10 ** (-model_params['learning_rate_exp'])
    regularization_rate = 10 ** (-model_params['regularization_rate_exp'])

    model = Sequential()

    model.add(BatchNormalization(input_shape=(dim_length, dim_channels)))

    filters = [model_params['cnn_filters']] * model_params['cnn_layers']

    for filter_number in filters:
        # model.add(MaxPooling1D(pool_size=5, strides=None, padding='same'))

        model.add(
            Convolution1D(filter_number,
                          kernel_size=model_params['kernel_size'],
                          padding='same',
                          kernel_regularizer=l2(regularization_rate),
                          kernel_initializer=weightinit))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

    model.add(Flatten())

    model.add(
        Dense(units=model_params['fc_nodes'],
              kernel_regularizer=l2(regularization_rate),
              kernel_initializer=weightinit))  # Fully connected layer
    model.add(Activation('relu'))  # Relu activation

    # Adding one more FC layer
    model.add(
        Dense(units=model_params['fc_nodes'],
              kernel_regularizer=l2(regularization_rate),
              kernel_initializer=weightinit))  # Fully connected layer
    model.add(Activation('relu'))  # Relu activation

    model.add(Dense(units=outputdim, kernel_initializer=weightinit))
    model.add(BatchNormalization())
    model.add(Activation("softmax"))  # Final classification layer

    model.compile(loss='categorical_crossentropy',
                  optimizer=Adam(lr=learning_rate),
                  metrics=['accuracy'])

    # i = 0
    # for model, params, model_types in [model]:
    #     logging.info('model ' + str(i))
    #     i = i + 1
    #     logging.info(params)
    #     logging.info(model.summary())

    return model


def train(model_fn, params, X_train, y_train, y_train_binary):
    # Design model
    logging.info('Creating model...')
    model = create_model(params['dim'], params['n_channels'],
                         params['n_classes'])
    logging.info(model.summary())
    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0,
                              patience=3,
                              verbose=1,
                              restore_best_weights=True)

    checkpoint = ModelCheckpoint(model_fn,
                                 monitor='val_loss',
                                 mode='min',
                                 save_best_only=True,
                                 verbose=1)

    # csv_logger = CSVLogger(os.path.join(channel_data_dir, 'training.log'))
    #
    # tbCallBack = TensorBoard(log_dir=os.path.join(channel_data_dir, 'Graph'),
    #                          histogram_freq=0,
    #                          write_graph=True,
    #                          write_images=True)

    callbacks = [earlystop, checkpoint]

    class_weights = compute_class_weight('balanced', np.unique(y_train), y_train)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    logging.info('Fitting model...')

    # Train model on dataset
    history = model.fit(
        X_train,
        y_train_binary,
        validation_split=model_params['validation_split'],
        batch_size=model_params['batch_size'],
        epochs=model_params['epochs'],
        shuffle=True,
        class_weight=class_weights,
        verbose=1,
        callbacks=callbacks)

    return model, history, X_train.shape[0], int(X_train.shape[0] *
                                                 model_params['validation_split'])


def cv_train_and_evaluate(X, y, y_binary, win_ids, train_indices, test_indices, model_dir, svtype):
    # Generate batches from indices
    X_train, X_test = X[train_indices], X[test_indices]
    y_train, y_test = y[train_indices], y[test_indices]
    y_train_binary, y_test_binary = y_binary[train_indices], y_binary[
        test_indices]
    win_ids_train, win_ids_test = win_ids[train_indices], win_ids[
        test_indices]

    # Parameters
    params = {
        'dim': X_train.shape[1],
        'n_classes': len(mapclasses.keys()),
        'n_channels': X_train.shape[2],
        'shuffle': True
    }

    os.makedirs(model_dir, exist_ok=True)
    model_fn = os.path.join(model_dir, 'model.hdf5')

    model, history, train_set_size, validation_set_size = train(
        model_fn, params, X_train, y_train, y_train_binary)

    model.save(model_fn)

    results = pd.DataFrame()

    intermediate_results, metrics = evaluate_model(model, X_test, y_test_binary, win_ids_test,
                                                   results, mapclasses, model_dir, svtype)

    results = results.append(intermediate_results)
    results.to_csv(os.path.join(model_dir, 'metrics.csv'), sep='\t')


def cross_validation(training_windows, outDir, npz_mode, svtype, kfold):
    X, y, win_ids = get_data(training_windows, npz_mode, svtype)
    y_binary = to_categorical(y, num_classes=len(mapclasses.keys()))

    # Instantiate the cross validator
    skf = StratifiedKFold(n_splits=kfold, shuffle=True, random_state=1)

    # Loop through the indices the split() method returns
    for index, (train_indices, test_indices) in enumerate(skf.split(X, y)):
        print("Training on fold " + str(index + 1) + "/" + str(kfold) + "...")

        model_dir = os.path.join(outDir, 'kfold', svtype,
                                 str(index + 1))

        cv_train_and_evaluate(X, y, y_binary, win_ids,
                              train_indices, test_indices, model_dir, svtype)


def cross_validation_by_chrom(training_windows, outDir, npz_mode, svtype, chrlist):
    X, y, win_ids = get_data(training_windows, npz_mode, svtype)
    y_binary = to_categorical(y, num_classes=len(mapclasses.keys()))

    # print(win_ids)
    chrom_num1 = map(lambda x: x.split('_')[0], win_ids)

    chrom_array = np.array([c for c in chrom_num1 if c in chrlist])
    # print(chrom_array)

    cv_dict = {}

    for c in np.unique(chrom_array):
        # print('Considering chromosome: {}'.format(c))

        idx_chr = np.where(chrom_array == c)
        idx_not_chr = np.where(chrom_array != c)

        cv_dict[c] = (idx_not_chr, idx_chr)

    # Loop through the indices the split() method returns
    for chrom in cv_dict.keys():
        train_indices, test_indices = cv_dict[chrom]

        print("Test on chromosome " + chrom + "...")

        model_dir = os.path.join(outDir, 'chrom', svtype, chrom)

        cv_train_and_evaluate(X, y, y_binary, win_ids,
                              train_indices, test_indices, model_dir, svtype)


def train_and_test_model(training_name, test_name, training_windows, test_windows,
                         outDir,
                         npz_mode, svtype):
    X_test, y_test, win_ids_test = get_data(training_windows, npz_mode, svtype)
    X_test, y_test, win_ids_test = get_data(test_windows, npz_mode, svtype)

    # Parameters
    params = {
        'dim': X_train.shape[1],
        'batch_size': model_params['batch_size'],
        'epochs': model_params['epochs'],
        'val_split': model_params['validation_split'],
        'n_classes': len(mapclasses.keys()),
        'n_channels': X_train.shape[2],
        'shuffle': True
    }

    y_train_binary = to_categorical(y_train, num_classes=params['n_classes'])
    y_test_binary = to_categorical(y_test, num_classes=params['n_classes'])

    model_dir = os.path.join(outDir, 'trained_on_' +
                             training_name + '_tested_on_' + test_name)
    os.makedirs(model_dir, exist_ok=True)
    model_fn = os.path.join(model_dir, 'model.hdf5')

    print('Training model on {}...'.format(training_name))
    model, history, train_set_size, validation_set_size = train(
        model_fn, params,
        X_train, y_train, y_train_binary)

    results = pd.DataFrame()

    intermediate_results, metrics = evaluate_model(model, X_test, y_test_binary, win_ids_test,
                                                   results, mapclasses, model_dir, svtype)

    results = results.append(intermediate_results)
    results.to_csv(os.path.join(model_dir, 'metrics.csv'), sep='\t')


def main():

    default_win = 25
    default_path = os.path.join('./cnn/win' + str(default_win), 'split_reads')
    def_windows_file = os.path.join(
        default_path, 'windows', 'DEL', 'windows_en.npz')
    parser = argparse.ArgumentParser(description='Train and test model')
    parser.add_argument('-p',
                        '--outputpath',
                        type=str,
                        default=default_path,
                        help="Specify output path")
    parser.add_argument('-t',
                        '--training_windows',
                        type=str,
                        default=def_windows_file + ',' + def_windows_file,
                        help="Comma separated list of training data")
    parser.add_argument('-x',
                        '--test_windows',
                        type=str,
                        default=def_windows_file + ',' + def_windows_file,
                        help="Specify training sample")
    parser.add_argument('-tn',
                        '--training_sample_name',
                        type=str,
                        default='git-data',
                        help="Specify training sample")
    parser.add_argument('-xn',
                        '--test_sample_name',
                        type=str,
                        default='git-data',
                        help="Specify training sample")
    parser.add_argument('-l',
                        '--logfile',
                        default='training.log',
                        help='File in which to write logs.')
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    parser.add_argument('-cv',
                        '--cv',
                        type=str,
                        default='kfold',
                        help="kfold/chrom")
    parser.add_argument('-c',
                        '--chrlist',
                        type=str,
                        default='12,22',
                        help="Comma separated list of chromosomes to consider")
    parser.add_argument('-k',
                        '--kfold',
                        type=int,
                        default=2,
                        help="Specify [k]-fold cross validation")
    parser.add_argument('-e',
                        '--epochs',
                        type=int,
                        default=1,
                        help="Number of epochs")
    parser.add_argument('-b',
                        '--batch_size',
                        type=int,
                        default=32,
                        help="Batch size")
    parser.add_argument('-val',
                        '--validation_split',
                        type=float,
                        default=0.2,
                        help="Percent of training set to use for validation")
    parser.add_argument('-npz',
                        '--load_npz',
                        type=bool,
                        default=True,
                        help="load npz?")
    parser.add_argument('-cnn_layers',
                        '--cnn_layers',
                        type=int,
                        default=4,
                        help="Number of convolutional layers")
    parser.add_argument('-cnn_filters',
                        '--cnn_filters',
                        type=int,
                        default=8,
                        help="Number of convolutional filters")
    parser.add_argument('-kernel_size',
                        '--kernel_size',
                        type=int,
                        default=7,
                        help="Number of convolutional filters")
    parser.add_argument('-fc_nodes',
                        '--fc_nodes',
                        type=int,
                        default=16,
                        help="Number of neurons in the dense layer")
    parser.add_argument('-learning_rate_exp',
                        '--learning_rate_exp',
                        type=int,
                        default=4,
                        help="Learning rate = 10 ** (-learning_rate_exp)")
    parser.add_argument('-regularization_rate_exp',
                        '--regularization_rate_exp',
                        type=int,
                        default=1,
                        help="Regularization rate = 10 ** (-regularization_rate_exp)")
    args = parser.parse_args()
    global mapclasses
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    global model_params
    model_params = {
        'batch_size': args.batch_size,
        'epochs': args.epochs,
        'validation_split': args.validation_split,
        'cnn_layers': args.cnn_layers,
        'cnn_filters': args.cnn_filters,
        'kernel_size': args.kernel_size,
        'fc_nodes': args.fc_nodes,
        'learning_rate_exp': args.learning_rate_exp,
        'regularization_rate_exp': args.regularization_rate_exp
    }
    output_dir = args.outputpath
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)

    print('Writing log file to {}'.format(logfilename))
    t0 = time()
    training_windows_list = args.training_windows.split(',')
    test_windows_list = args.test_windows.split(',')

    for t in training_windows_list:
        assert os.path.exists(t)

    for t in test_windows_list:
        assert os.path.exists(t)

    if len(set(test_windows_list) & set(training_windows_list)) == 0:
        train_and_test_model(
            training_name=args.training_sample_name,
            test_name=args.test_sample_name,
            training_windows=training_windows_list,
            test_windows=args.test_windowsr,
            outDir=output_dir,
            npz_mode=args.load_npz,
            svtype=args.svtype
        )
    else:
        if args.cv == 'kfold':
            cross_validation(training_windows=training_windows_list,
                             outDir=output_dir,
                             npz_mode=args.load_npz,
                             svtype=args.svtype,
                             kfold=args.kfold)
        elif args.cv == 'chrom':
            cross_validation_by_chrom(training_windows=training_windows_list,
                                      outDir=output_dir,
                                      npz_mode=args.load_npz,
                                      svtype=args.svtype,
                                      chrlist=args.chrlist)
    logging.info('Elapsed time training and testing = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
