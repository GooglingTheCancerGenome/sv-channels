import argparse
import logging
import sys
from skopt import gp_minimize
from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
import tensorflow as tf
import numpy as np
from time import time
from sklearn.utils.class_weight import compute_class_weight
from sklearn.model_selection import train_test_split
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.regularizers import l2
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import (EarlyStopping, ModelCheckpoint,
                                        TensorBoard)
from tensorflow.keras.layers import (Activation, BatchNormalization,
                                     Convolution1D, Dense, Flatten)
from tensorflow.keras.models import Sequential, load_model

sys.path.append('../genome_wide/')
from functions import load_windows


dim_cnn_filters = Integer(low=4, high=16, name='cnn_filters')
dim_cnn_layers = Integer(low=1, high=6, name='cnn_layers')
dim_cnn_kernel_size = Integer(low=1, high=8, name='cnn_kernel_size')
dim_cnn_fc_nodes = Integer(low=4, high=10, name='cnn_fc_nodes')
dim_init_learning_rate = Real(low=1e-4, high=1e-1, prior='log-uniform', name='cnn_init_learning_rate')
dim_regularization_rate = Real(low=1e-4, high=1e-1, prior='log-uniform', name='cnn_regularization_rate')

dimensions = [dim_cnn_filters, dim_cnn_layers, dim_cnn_kernel_size, dim_cnn_fc_nodes,
              dim_init_learning_rate, dim_regularization_rate]

default_parameters = [8, 1, 7, 6, 1e-4, 1e-1]

path_best_model = './best_model.keras'
best_accuracy = 0.0


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



@use_named_args(dimensions=dimensions)
def fitness(cnn_filters, cnn_layers, cnn_kernel_size, cnn_fc_nodes,
            cnn_init_learning_rate, cnn_regularization_rate):
    global best_accuracy
    # best_accuracy = 0.0
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

    history = model.fit(x=train_X[:3000], y=train_y[:3000],
                        epochs=max_epoch, batch_size=batch_size,
                        shuffle=True,
                        validation_data=validation_data,
                        class_weight=class_weights,
                        verbose=0,
                        callbacks=callbacks)

    accuracy = history.history['val_accuracy'][-1]
    print()
    print('Accuracy: {0:.2%}'.format(accuracy))
    if accuracy > best_accuracy:
        model.save(path_best_model)
        best_accuracy = accuracy
    del model
    tf.keras.backend.clear_session()
    return -accuracy


def optimize(args):

    global train_X, val_X, train_y, val_y, class_weights, batch_size, max_epoch

    randomState = 46
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    batch_size = args.batch_size
    max_epoch = args.epochs

    X, y = load_windows(args.windows)
    y = y.values()
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    classes = np.array(np.unique(y))
    y_lab = np.asarray(y)
    class_weights = compute_class_weight('balanced', classes, y_lab)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    y = to_categorical(y, num_classes=2)

    train_X, val_X, train_y, val_y = train_test_split(
        X, y, test_size=args.validation_split, random_state=2, stratify=y, shuffle=True)

    search_result = gp_minimize(func=fitness, dimensions=dimensions, acq_func='EI',
                                n_calls=200, x0=default_parameters, random_state=7, n_jobs=-1)

    hyps = np.asarray(search_result.x)
    np.save('./hyperparams.npy', hyps, allow_pickle=False)


def main():

    parser = argparse.ArgumentParser(description='Optimize model')

    parser.add_argument('-w',
                        '--windows',
                        type=str,
                        default='../genome_wide/cnn/win25/split_reads/windows/DEL/windows_en.npz',
                        help="Comma separated list of training data")
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
    parser.add_argument('-n',
                        '--ncalls',
                        type=int,
                        default=50,
                        help="Number of calls of the fitness function")
    parser.add_argument('-val',
                        '--validation_split',
                        type=float,
                        default=0.3,
                        help="Percent of training set to use for validation")
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Type of SV")
    args = parser.parse_args()

    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(format=log_format,
                        filename=args.logfile,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    optimize(args)
    logging.info('Elapsed time = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
