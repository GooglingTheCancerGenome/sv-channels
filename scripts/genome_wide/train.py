from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import gzip
import json
import logging
import os
from time import time
import joblib

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.utils.class_weight import compute_class_weight

from tensorflow.keras.callbacks import (EarlyStopping, ModelCheckpoint,
                                        TensorBoard)
from tensorflow.keras.layers import (Activation, BatchNormalization,
                                     Convolution1D, Dense, Dropout, Flatten,
                                     Lambda, Reshape, TimeDistributed,
                                     concatenate, Input)
from tensorflow.keras.models import Sequential, Model, load_model
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras.utils import to_categorical

from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

from skopt.utils import use_named_args  # decorator to convert a list of parameters to named arguments
from skopt import gp_minimize  # Bayesian optimization using Gaussian Processes

from sklearn.metrics import f1_score

from skopt.space import Real, Integer

from model_functions import (  # create_model_with_mcfly, train_model_with_mcfly
    evaluate_model, get_data)

counter = 0


# onstep and make_objective functions adapted from the notebook "Bayesian Optimization Workshop" by Luca Massaron
# https://colab.research.google.com/github/lmassaron/kaggledays-2019-gbdt/blob/master/skopt_workshop_part2.ipynb#scrollTo=Gc9IfCbORvux


def onstep(res):
    global counter
    x0 = res.x_iters  # List of input points
    y0 = res.func_vals  # Evaluation of input points
    print('Last eval: ', x0[-1],
          ' - Score ', y0[-1])
    print('Current iter: ', counter,
          ' - Score ', res.fun,
          ' - Args: ', res.x)
    joblib.dump((x0, y0), 'checkpoint.pkl')  # Saving a checkpoint to disk
    counter += 1


def make_objective(create_model, X, y, space, cv, scoring, class_weights, callbacks,
                   train_params, verbose=0):
    # This decorator converts your objective function with named arguments into one that
    # accepts a list as argument, while doing the conversion automatically.
    @use_named_args(space)
    def objective(**params):
        arch_params = create_model.__code__.co_varnames

        objective_model_params = {p: params[p] for p in params if p in arch_params}
        fit_params = {p: params[p] for p in params if p not in arch_params}

        model = KerasClassifier(build_fn=create_model2,
                                outputdim=train_params['n_classes'],
                                shape=X.shape,
                                verbose=verbose,
                                **objective_model_params)

        score = list()
        print("Testing parameters:", params)
        print(model.sk_params)

        win1_end = int(X.shape[1] / 2 - 5)
        win2_start = int(X.shape[1] / 2 + 5)

        for j, (train_index, test_index) in enumerate(cv.split(X, y)):
            model.fit([X[train_index, :win1_end], X[train_index, win2_start:]],
                      to_categorical(y[train_index]),
                      epochs=model_params['epochs'],
                      batch_size=model_params['batch_size'],
                      shuffle=True,
                      validation_split=model_params['validation_split'],
                      verbose=0,
                      callbacks=callbacks,
                      class_weight=class_weights,
                      **fit_params)
            val_preds = model.model.predict([X[test_index, :win1_end], X[test_index, win2_start:]])
            val_preds = np.argmax(val_preds, axis=1)
            score.append(scoring(y_true=y[test_index], y_pred=val_preds, average='weighted'))

        print("CV scores:", score)
        return np.mean(score) * -1

    return objective


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


def create_model(shape, outputdim, learning_rate, regularization_rate,
                 filters, layers, kernel_size, fc_nodes):
    weightinit = 'lecun_uniform'  # weight initialization

    model = Sequential()

    model.add(BatchNormalization(input_shape=(shape[1], shape[2])))

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
    model.add(Activation("softmax"))  # Final classification layer

    model.compile(loss='categorical_crossentropy',
                  optimizer=Adam(lr=learning_rate),
                  metrics=['accuracy'])

    return model


def create_model2(shape, outputdim, learning_rate, regularization_rate,
                  filters, layers, kernel_size, fc_nodes):
    weightinit = 'lecun_uniform'  # weight initialization

    dim1 = int(shape[1] / 2 - 5)
    input_A = Input(shape=(dim1, shape[2]), name="window1")
    input_B = Input(shape=(dim1, shape[2]), name="window2")

    filters_list = [filters] * layers

    layers = []
    for input in [input_A, input_B]:
        layer1 = BatchNormalization()(input)
        for filter_number in filters_list:
            layer1 = Convolution1D(filter_number,
                                   kernel_size=(kernel_size,),
                                   padding='same',
                                   kernel_regularizer=l2(regularization_rate),
                                   kernel_initializer=weightinit)(layer1)
            layer1 = BatchNormalization()(layer1)
            layer1 = Activation('relu')(layer1)
        layer1 = Flatten()(layer1)
        layers.append(layer1)

    concat = concatenate([layers[0], layers[1]])

    layer2 = Dense(units=fc_nodes,
                   kernel_regularizer=l2(regularization_rate),
                   kernel_initializer=weightinit)(concat)
    layer2 = Activation('relu')(layer2)  # Relu activation

    layer2 = Dense(units=outputdim, kernel_initializer=weightinit)(layer2)
    layer2 = BatchNormalization()(layer2)

    output = Activation("softmax")(layer2)  # Final classification layer

    model = Model(inputs=[input_A, input_B],
                  outputs=[output])

    model.compile(loss='categorical_crossentropy',
                  optimizer=Adam(lr=learning_rate),
                  metrics=['accuracy'])

    return model


def train(model_fn, model_fn_checkpoint, train_params, X_train, y_train, y_train_binary):
    dimensions = [
        Integer(low=4, high=16, name='filters'),
        Integer(low=1, high=6, name='layers'),
        Integer(low=1, high=8, name='kernel_size'),
        Integer(low=4, high=10, name='fc_nodes'),
        Real(low=1e-4, high=1e-1, prior='log-uniform', name='learning_rate'),
        Real(low=1e-4, high=1e-1, prior='log-uniform', name='regularization_rate')
    ]

    # Setting a 5-fold stratified cross-validation (note: shuffle=True)
    skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=0)

    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0,
                              patience=3,
                              verbose=0,
                              restore_best_weights=True)

    checkpoint = ModelCheckpoint(model_fn_checkpoint,
                                 monitor='val_loss',
                                 mode='min',
                                 save_best_only=True,
                                 verbose=0)

    # csv_logger = CSVLogger(os.path.join(channel_data_dir, 'training.log'))
    #
    # tbCallBack = TensorBoard(log_dir=os.path.join(channel_data_dir, 'Graph'),
    #                          histogram_freq=0,
    #                          write_graph=True,
    #                          write_images=True)

    callbacks = [earlystop, checkpoint]

    class_weights = compute_class_weight('balanced', np.unique(y_train), y_train)
    class_weights = {i: v for i, v in enumerate(class_weights)}

    if model_params['optimize']:
        objective = make_objective(create_model2,
                                   X_train, y_train,
                                   space=dimensions,
                                   cv=skf,
                                   scoring=f1_score,
                                   class_weights=class_weights,
                                   callbacks=callbacks,
                                   train_params=train_params,
                                   verbose=0)

        gp_round = gp_minimize(func=objective,
                               dimensions=dimensions,
                               acq_func='gp_hedge',  # Expected Improvement.
                               n_calls=model_params['gp_calls'],
                               callback=[onstep],
                               random_state=1)

        best_parameters = gp_round.x
        best_result = gp_round.fun
        print(best_parameters, best_result)

        model = create_model2(X_train.shape,
                              outputdim=train_params['n_classes'],
                              learning_rate=best_parameters[4],
                              regularization_rate=best_parameters[5],
                              filters=best_parameters[0],
                              layers=best_parameters[1],
                              kernel_size=best_parameters[2],
                              fc_nodes=best_parameters[3]
                              )

    else:

        model = create_model2(X_train.shape,
                              outputdim=train_params['n_classes'],
                              learning_rate=model_params['learning_rate'],
                              regularization_rate=model_params['regularization_rate'],
                              filters=model_params['cnn_filters'],
                              layers=model_params['cnn_layers'],
                              kernel_size=model_params['kernel_size'],
                              fc_nodes=model_params['fc_nodes']
                              )

    win1_end = int(X_train.shape[1] / 2 - 5)
    win2_start = int(X_train.shape[1] / 2 + 5)

    model.fit(x=[X_train[:, :win1_end, :], X_train[:, win2_start:, :]], y=y_train_binary,
              epochs=model_params['epochs'], batch_size=model_params['batch_size'],
              shuffle=True,
              validation_split=model_params['validation_split'],
              class_weight=class_weights,
              verbose=1,
              callbacks=callbacks)

    return X_train.shape[0], int(X_train.shape[0] *
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
    model_fn_checkpoint = os.path.join(model_dir, 'model.checkpoint.hdf5')

    train_set_size, validation_set_size = train(
        model_fn, model_fn_checkpoint, params, X_train, y_train, y_train_binary)

    model = load_model(model_fn_checkpoint)

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
    model_fn_checkpoint = os.path.join(model_dir, 'model.checkpoint.hdf5')

    print('Training model on {}...'.format(training_name))
    model, history, train_set_size, validation_set_size = train(
        model_fn, model_fn_checkpoint, params,
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
                        default=3,
                        help="Number of convolutional layers")
    parser.add_argument('-cnn_filters',
                        '--cnn_filters',
                        type=int,
                        default=15,
                        help="Number of convolutional filters")
    parser.add_argument('-kernel_size',
                        '--kernel_size',
                        type=int,
                        default=4,
                        help="Number of convolutional filters")
    parser.add_argument('-fc_nodes',
                        '--fc_nodes',
                        type=int,
                        default=6,
                        help="Number of neurons in the dense layer")
    parser.add_argument('-learning_rate',
                        '--learning_rate',
                        type=float,
                        default=3.80247940e-04,
                        help="initial learning rate")
    parser.add_argument('-regularization_rate',
                        '--regularization_rate',
                        type=float,
                        default=2.00180567e-04,
                        help="regularization rate")
    parser.add_argument('-gp_calls',
                        '--gp_calls',
                        type=int,
                        default=10,
                        help="Number of calls of gp_minimize")
    parser.add_argument('-optimize',
                        '--optimize',
                        type=bool,
                        default=False,
                        help="Use gp_minimize?")

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
        'learning_rate': args.learning_rate,
        'regularization_rate': args.regularization_rate,
        'gp_calls': args.gp_calls,
        'optimize': args.optimize
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
