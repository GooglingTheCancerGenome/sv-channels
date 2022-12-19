"""
Leave one chromosome out procedure
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
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import (EarlyStopping, ModelCheckpoint,
                                        TensorBoard)
from tensorflow.keras.layers import (Activation, BatchNormalization,
                                     Convolution1D, Dense, Flatten)
from tensorflow.keras.models import Sequential, load_model


def load_windows(win_file, lab_file):
    X = zarr.load(win_file)
    with gzip.GzipFile(lab_file, 'r') as fin:
        y = json.loads(fin.read().decode('utf-8'))
    return X, y


def create_model(X, outputdim, learning_rate, regularization_rate,
                 filters, layers, kernel_size, fc_nodes):

    weightinit = 'lecun_uniform'  # weight initialization

    model = Sequential()

    model.add(BatchNormalization(input_shape=(X.shape[1], X.shape[2])))

    filters_list = [filters] * layers

    for filter_number in filters_list:
        model.add(
            Convolution1D(filter_number,
                          kernel_size=(kernel_size,),
                          padding='same',
                          kernel_regularizer=l2(regularization_rate),
                          kernel_initializer=weightinit))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

    model.add(Flatten())

    model.add(
        Dense(units=fc_nodes,
              kernel_regularizer=l2(regularization_rate),
              kernel_initializer=weightinit))  # Fully connected layer
    model.add(Activation('relu'))  # Relu activation

    model.add(Dense(units=outputdim, kernel_initializer=weightinit))
    model.add(BatchNormalization())
    model.add(Activation("sigmoid"))  # Final classification layer

    model.compile(loss='categorical_crossentropy',
                  optimizer=Adam(lr=learning_rate),
                  metrics=['accuracy'])

    return model


def train(args):

    global train_X, val_X, train_y, val_y, class_weights, batch_size, max_epoch, path_best_model

    randomState = 46
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    batch_size = args.batch_size
    max_epoch = args.epochs
    path_best_model = args.model

    windows = args.windows.split(',')
    labels = args.labels.split(',')

    X = []
    y = []
    for w, l in zip(windows, labels):
        partial_X, partial_y = load_windows(w, l)
        X.extend(partial_X)
        y.extend(partial_y.values())

    train_X = np.stack(X, axis=0)

    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    classes = np.array(np.unique(y))
    y_lab = np.asarray(y)
    class_weights = compute_class_weight('balanced', classes, y_lab)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    train_y = to_categorical(y, num_classes=2)

    val_X, y = load_windows(args.validation_windows, args.validation_labels)
    y = y.values()
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    classes = np.array(np.unique(y))
    y_lab = np.asarray(y)
    class_weights = compute_class_weight('balanced', classes, y_lab)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    val_y = to_categorical(y, num_classes=2)

    hparams = np.load(args.hparams)

    cnn_filters = int(hparams[0])
    cnn_layers = int(hparams[1])
    cnn_kernel_size = int(hparams[2])
    cnn_fc_nodes = int(hparams[3])
    cnn_init_learning_rate = hparams[4]
    cnn_regularization_rate = hparams[5]

    print()
    print('cnn_filters: ', cnn_filters)
    print('cnn_layers: ', cnn_layers)
    print('cnn_kernel_size: ', cnn_kernel_size)
    print('cnn_fc_nodes: ', cnn_fc_nodes)
    print('cnn_init_learning_rate: ', cnn_init_learning_rate)
    print('cnn_regularization_rate: ', cnn_regularization_rate)
    print()

    model = create_model(train_X, 2,
                         learning_rate=cnn_init_learning_rate, regularization_rate=cnn_regularization_rate,
                         filters=cnn_filters, layers=cnn_layers, kernel_size=cnn_kernel_size, fc_nodes=cnn_fc_nodes)
    print(model.summary())

    callback_log = TensorBoard(
        log_dir='log_dir',
        histogram_freq=0,
        batch_size=32,
        write_graph=True,
        write_grads=True,
        write_images=False)

    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0,
                              patience=3,
                              verbose=1,
                              restore_best_weights=True)

    callbacks = [callback_log, earlystop]

    validation_data = (val_X, val_y)

    history = model.fit(x=train_X, y=train_y,
                        epochs=max_epoch, batch_size=batch_size,
                        shuffle=True,
                        validation_data=validation_data,
                        class_weight=class_weights,
                        verbose=0,
                        callbacks=callbacks)

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
                        help="Comma separated list of JSON.GZ file for training labels")
    parser.add_argument('-val_win',
                        '--validation_windows',
                        type=str,
                        default='sv_chan.zarr,sv_chan.zarr',
                        help="Comma separated list of validation data")
    parser.add_argument('-val_lab',
                        '--validation_labels',
                        type=str,
                        default='labels/labels.json.gz,labels/labels.json.gz',
                        help="Comma separated list of JSON.GZ file for validation labels")
    parser.add_argument('-sm',
                        '--samples',
                        type=str,
                        default='SAMPLE1,SAMPLE2',
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
