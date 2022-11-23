import os
import argparse
import logging

import matplotlib.pyplot as plt
import pandas as pd
import json
import vcf
from matplotlib import pyplot as plt
import tensorflow as tf
import numpy as np
from time import time
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import auc, precision_recall_curve
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.models import load_model
from collections import Counter

import sys

# setting path
if '__file__' in vars():
    # print("We are running the script non interactively")
    path = os.path.join(os.path.dirname(__file__), os.pardir)
    sys.path.append(path)
else:
    # print('We are running the script interactively')
    sys.path.append("../..")

# importing
from model_functions import create_model

# Search space for the hyperparameters

search_space_dict = {
    'cnn_filters': Integer(low=1, high=8, name='cnn_filters'),
    'cnn_layers': Integer(low=1, high=2, name='cnn_layers'),
    'cnn_filter_size': Integer(low=1, high=8, name='cnn_filter_size'),
    'fc_nodes': Integer(low=3, high=6, name='fc_nodes'),
    'dropout_rate': Real(low=0, high=0.9, name='dropout_rate'),
    'learning_rate': Real(low=1e-4, high=1e-1, prior='log-uniform', name='learning_rate'),
    'regularization_rate': Real(low=1e-4, high=1e-1, prior='log-uniform', name='regularization_rate')
}

dimensions = [search_space_dict[k] for k in search_space_dict.keys()]

default_parameters = [4, 1, 7, 4, 0.2, 1e-4, 1e-1]

# Create a global variable that stores the current best area under the ROC curve (AUC), our metric. Initialize it to 0
best_auc = 0.0

best_hyperparameters = {}


def load_data(inputfile):
    npz_archive = np.load(inputfile)
    X = npz_archive['X']
    y = npz_archive['y']

    logging.info('Data loaded')
    logging.info('Shape of X: {}'.format(X.shape))
    logging.info('Shape of y: {}'.format(y.shape))
    logging.info('Labels => samples: {}'.format(Counter(y[:, 0])))
    logging.info('Labels => 1st chr: {}'.format(Counter(y[:, 1])))
    logging.info('Labels => 2nd chr: {}'.format(Counter(y[:, 3])))

    return X, y


def train_test_val_split_by_sample_and_chrom(X, y, test_sample, test_chrom, val_chrom):
    idx_train = np.where(np.logical_and(y[:, 0] != test_sample,
                                        y[:, 1] != test_chrom,
                                        y[:, 1] != val_chrom))[0]
    idx_test = np.where(np.logical_and(y[:, 0] == test_sample, y[:, 1] == test_chrom))[0]
    idx_val = np.where(np.logical_and(y[:, 0] == test_sample, y[:, 1] == val_chrom))[0]

    idx_len_dict = {'train': len(idx_train), 'test': len(idx_test), 'val': len(idx_val)}

    print('Lengths: {}'.format(idx_len_dict))

    return X[idx_train], X[idx_test], X[idx_val], y[idx_train], y[idx_test], y[idx_val]


def train_test_split_by_sample_and_chrom(X, y, test_sample, test_chrom):
    idx_train = np.where(np.logical_and(y[:, 0] != test_sample,
                                        y[:, 1] != test_chrom))[0]
    idx_test = np.where(np.logical_and(y[:, 0] == test_sample, y[:, 1] == test_chrom))[0]

    idx_len_dict = {'train': len(idx_train), 'test': len(idx_test)}

    print('Lengths: {}'.format(idx_len_dict))

    return X[idx_train], X[idx_test], y[idx_train], y[idx_test]


def select_by_sample_and_chrom(X, y, sample, chrom):
    idx = np.where(np.logical_and(y[:, 0] == sample,
                                  y[:, 1] == chrom))[0]

    idx_len_dict = {'set': len(idx)}

    print('Lengths: {}'.format(idx_len_dict))

    return X[idx], y[idx]


@use_named_args(dimensions=dimensions)
def fitness(cnn_filters, cnn_layers, cnn_filter_size, fc_nodes,
            dropout_rate, learning_rate, regularization_rate):
    # Print the current combination of hyperparameters
    print('cnn_filters: {}'.format(cnn_filters))
    print('cnn_layers: {}'.format(cnn_layers))
    print('cnn_filter_size: {}'.format(cnn_filter_size))
    print('fc_nodes: {}'.format(fc_nodes))
    print('dropout_rate: {}'.format(dropout_rate))
    print('learning_rate: {}'.format(learning_rate))
    print('regularization_rate: {}'.format(regularization_rate))

    model = create_model(inner_X_train, 2,
                         learning_rate=learning_rate,
                         regularization_rate=regularization_rate,
                         filters=cnn_filters,
                         layers=cnn_layers,
                         kernel_size=cnn_filter_size,
                         fc_nodes=fc_nodes,
                         dropout_rate=dropout_rate)
    # print(model.summary())

    history = model.fit(x=inner_X_train, y=y_train,
                        epochs=n_epochs,
                        batch_size=model_batch_size,
                        shuffle=True,
                        validation_data=(inner_X_val, y_val),
                        class_weight=class_weights_train,
                        verbose=0)

    # Get the predicited probability of testing data
    y_score = model.predict(inner_X_test)[:, 0]
    # Data to plot precision - recall curve
    precision, recall, thresholds = precision_recall_curve(y_test, y_score)
    # Use AUC function to calculate the area under the curve of precision recall curve
    auc_precision_recall = auc(recall, precision)
    print('auc_precision_recall = {}'.format(auc_precision_recall))

    # Save the model if it improves on the best found performance which is stored by the global variable best_auc
    global best_auc
    global best_hyperparameters

    # If the AUC of the saved model is greater than the current best performance
    if auc_precision_recall > best_auc:
        # Save the new model
        model.save(path_best_model)

        # Store the best hyperparameters
        best_hyperparameters = {
            'innercv_test_sample': innercv_test_sample,
            'outer_i': outer_i,
            'outer_chrom': outer_c,
            'inner_i': inner_i,
            'inner_chrom': inner_c,
            'cnn_filters': cnn_filters,
            'cnn_layers': cnn_layers,
            'cnn_filter_size': cnn_filter_size,
            'fc_nodes': fc_nodes,
            'dropout_rate': dropout_rate,
            'learning_rate': learning_rate,
            'regularization_rate': regularization_rate,
            'auc_precision_recall': auc_precision_recall
        }
        best_hyperparameters = {k: [v] for k, v in best_hyperparameters.items()}

        # Update the current greatest AUC score
        best_auc = auc_precision_recall

        prefix = '_'.join(['outer', outer_c, 'inner', inner_c, 'cnn-filt', str(cnn_filters),
                           'cnn-lay', str(cnn_layers), 'cnn-filt-size', str(cnn_filter_size),
                           'fc-nodes', str(fc_nodes), 'dropout', str(dropout_rate),
                           'lr', str(learning_rate), 'rr', str(regularization_rate),
                           'pr-auc', str(auc_precision_recall)])
        accuracy_plot = os.path.join(os.path.dirname(path_best_model),
                                     ''.join([prefix, '_accuracy.png']))
        loss_plot = os.path.join(os.path.dirname(path_best_model),
                                 ''.join([prefix, '_loss.png']))

        plt.plot(history.history['accuracy'])
        plt.plot(history.history['val_accuracy'])
        plt.title('model accuracy')
        plt.ylabel('accuracy')
        plt.xlabel('epoch')
        plt.legend(['train', 'val'], loc='upper right')
        plt.savefig(accuracy_plot)

        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'val'], loc='upper left')
        plt.savefig(loss_plot)

        # Else delete the model that just finishing training from meomory
        del model

        # Clear the Keras session, otherwise it will keep adding new
        # models to the same TensorFlow graph each time we create
        # a model with a different set of hyperparameters.
        tf.keras.backend.clear_session()

        # Scikit-optimize does minimization so it tries to
        # find a set the combination of hyperparameters with the lowest fitness value.
        # We want to maximize AUC so we negate this number.
    return -auc_precision_recall


def train(input_args):
    def get_labels(y):

        svtype = 'DEL'
        mapclasses = {svtype: 0, 'no' + svtype: 1}
        y_mapped = np.array([mapclasses[i] for i in y[:, 6]])
        classes = np.array(np.unique(y_mapped))
        y_lab = np.asarray(y_mapped)
        class_weights = compute_class_weight(class_weight='balanced', classes=classes, y=y_lab)
        class_weights = {i: v for i, v in enumerate(class_weights)}
        y_cat = to_categorical(y=y_lab, num_classes=2)

        return y_lab, y_cat, class_weights

    # Store chromosome and position of DELs as keys and difference of posterior probabilities as values
    pos_dict = {}
    model_df = pd.DataFrame()

    X, y = load_data(input_args.input)

    # test sample for the outer cross-validation
    outercv_test_sample = input_args.outercv_test_sample
    # test sample for the inner cross-validation
    global innercv_test_sample
    innercv_test_sample = input_args.innercv_test_sample
    # chromosome for validation
    val_chrom = input_args.val_chrom

    # list of chromosomes
    outer_chr_list = ['chr' + str(i) for i in np.arange(1, 23)]
    outer_chr_list.remove(val_chrom)

    logging.info('Outer CV')
    global outer_i, outer_c

    global best_auc
    best_auc = 0.0
    for outer_i, outer_c in enumerate([input_args.test_chrom]):

        logging.info('Considering outer test chromosome {}'.format(outer_c))

        outer_X_train, outer_X_test, outer_X_val, \
        outer_y_train, outer_y_test, outer_y_val = train_test_val_split_by_sample_and_chrom(
            X, y, outercv_test_sample, outer_c, val_chrom
        )

        assert outercv_test_sample not in outer_y_train[:, 0], "Test sample {} in training set: {}".format(
            outercv_test_sample, np.unique(outer_y_train[:, 0])
        )
        assert outercv_test_sample == np.unique(outer_y_test[:, 0]), "Test sample {} not in test set: {}".format(
            outercv_test_sample, np.unique(outer_y_test[:, 0])
        )
        assert outercv_test_sample == np.unique(outer_y_val[:, 0]), "Test sample {} not in val set: {}".format(
            outercv_test_sample, np.unique(outer_y_val[:, 0])
        )

        assert outer_c not in outer_y_train[:, 1], "Test chrom {} in training set: {}".format(
            outer_c, np.unique(outer_y_train[:, 1])
        )
        assert outer_c == np.unique(outer_y_test[:, 1]), "Test chrom {} not in test set: {}".format(
            outer_c, np.unique(outer_y_test[:, 1])
        )
        assert val_chrom == np.unique(outer_y_val[:, 1]), "Val chrom {} not in val set: {}".format(
            val_chrom, np.unique(outer_y_val[:, 1])
        )

        assert outer_X_train.shape[0] == outer_y_train.shape[0], "outer_X_train shape:{} different from" \
                                                                 "outer_y_train shape:{}".format(
            outer_X_train.shape[0],
            outer_y_train.shape[0]
        )
        assert outer_X_test.shape[0] == outer_y_test.shape[0], "outer_X_test shape:{} different from" \
                                                               "outer_y_test shape:{}".format(
            outer_X_test.shape[0],
            outer_y_test.shape[0]
        )
        assert outer_X_val.shape[0] == outer_y_val.shape[0], "outer_X_val shape:{} different from" \
                                                             "outer_y_val shape:{}".format(
            outer_X_val.shape[0],
            outer_y_val.shape[0]
        )

        inner_chr_list = outer_chr_list
        inner_chr_list.remove(outer_c)

        logging.info('Inner CV')

        global inner_i, inner_c

        for inner_i, inner_c in enumerate(inner_chr_list):
            logging.info('Considering inner test chromosome {}'.format(inner_c))

            global inner_X_train, inner_X_val, inner_X_test, \
                y_train, y_val, y_test, \
                class_weights_train, class_weights_val, class_weights_test, \
                n_epochs, model_batch_size, path_best_model

            n_epochs = input_args.epochs
            model_batch_size = input_args.batch_size
            path_best_model = os.path.join(input_args.output, ''.join([input_args.test_chrom, '_best_model.h5']))

            accuracy_plot_file = os.path.join(input_args.output,
                                              ''.join([input_args.test_chrom, '_model_accuracy.png']))
            loss_plot_file = os.path.join(input_args.output,
                                          ''.join([input_args.test_chrom, '_model_loss.png']))

            inner_X_train, inner_X_test, \
            inner_y_train, inner_y_test, = train_test_split_by_sample_and_chrom(
                outer_X_train, outer_y_train, innercv_test_sample, inner_c
            )

            assert innercv_test_sample not in inner_y_train[:, 0], "Test sample {} in training set: {}".format(
                innercv_test_sample, np.unique(inner_y_train[:, 0])
            )
            assert innercv_test_sample == np.unique(inner_y_test[:, 0]), "Test sample {} not in test set: {}".format(
                innercv_test_sample, np.unique(inner_y_test[:, 0])
            )

            assert inner_c not in inner_y_train[:, 1], "Test chrom {} in training set: {}".format(
                inner_c, np.unique(inner_y_train[:, 1])
            )
            assert inner_c == np.unique(inner_y_test[:, 1]), "Test chrom {} not in test set: {}".format(
                inner_c, np.unique(inner_y_test[:, 1])
            )

            assert inner_X_train.shape[0] == inner_y_train.shape[0], "inner_X_train shape:{} different from" \
                                                                     "inner_y_train shape:{}".format(
                inner_X_train.shape[0],
                inner_y_train.shape[0]
            )
            assert inner_X_test.shape[0] == inner_y_test.shape[0], "inner_X_test shape:{} different from" \
                                                                   "inner_y_test shape:{}".format(
                inner_X_test.shape[0],
                inner_y_test.shape[0]
            )

            inner_X_val, inner_y_val = select_by_sample_and_chrom(
                X, y, innercv_test_sample, val_chrom
            )
            assert innercv_test_sample == np.unique(inner_y_val[:, 0]), "Test sample {} not in val set: {}".format(
                innercv_test_sample, np.unique(inner_y_val[:, 0])
            )
            assert val_chrom == np.unique(inner_y_val[:, 1]), "Test chrom {} not in test set: {}".format(
                val_chrom, np.unique(inner_y_val[:, 1])
            )
            assert inner_X_val.shape[0] == inner_y_val.shape[0], "inner_X_val shape:{} different from" \
                                                                 "inner_y_val shape:{}".format(
                inner_X_val.shape[0],
                inner_y_val.shape[0]
            )

            _, y_train, class_weights_train = get_labels(inner_y_train)
            _, y_val, class_weights_val = get_labels(inner_y_val)
            y_test, _, class_weights_test = get_labels(inner_y_test)

            search_result = gp_minimize(func=fitness,
                                        dimensions=dimensions,
                                        acq_func='EI',
                                        n_calls=input_args.ncalls,
                                        x0=default_parameters,
                                        random_state=42,
                                        n_jobs=-1)

        # out of the inner CV loop
        df = pd.DataFrame.from_dict(best_hyperparameters)
        model_df = model_df.append(df, ignore_index=True)

        model = create_model(outer_X_train, 2,
                             learning_rate=float(best_hyperparameters['learning_rate'][0]),
                             regularization_rate=float(best_hyperparameters['regularization_rate'][0]),
                             filters=int(best_hyperparameters['cnn_filters'][0]),
                             layers=int(best_hyperparameters['cnn_layers'][0]),
                             kernel_size=int(best_hyperparameters['cnn_filter_size'][0]),
                             fc_nodes=int(best_hyperparameters['fc_nodes'][0]),
                             dropout_rate=float(best_hyperparameters['dropout_rate'][0])
                             )
        # print(model.summary())

        _, y_train, class_weights_train = get_labels(outer_y_train)
        _, y_val, class_weights_val = get_labels(outer_y_val)

        history = model.fit(x=outer_X_train, y=y_train,
                            epochs=n_epochs,
                            batch_size=model_batch_size,
                            shuffle=True,
                            validation_data=(outer_X_val, y_val),
                            class_weight=class_weights_train,
                            verbose=1)

        probs = model.predict(outer_X_test, batch_size=100, verbose=True)

        for i in np.arange(outer_y_val.shape[0]):
            chrom1 = outer_y_val[i, 1]
            pos1a = outer_y_val[i, 2]
            pos_dict[chrom1 + '_' + str(int(pos1a) + 1)] = str(probs[i][0] - probs[i][1])

    model_df_file = os.path.join(input_args.output, ''.join([input_args.test_chrom, '_model_df.csv']))
    model_df.to_csv(model_df_file)

    return pos_dict


def main():
    randomState = 42
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    input_dir = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels' + \
                '/sv-channels_manuscript/UMCU_hpc/'
    output_dir = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels' + \
                 '/sv-channels_manuscript/UMCU_hpc/nested_locso-cv'

    svchans = os.path.join(input_dir, '8_samples_channels.npz')

    parser = argparse.ArgumentParser(description='Optimize the hyperparameters of the model using a nested'
                                                 'cross-validation approach')

    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default=svchans,
                        help="File with the channels (X) and labels (y) in npz format")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='.',
                        help="Output folder")
    parser.add_argument('-ots',
                        '--outercv_test_sample',
                        type=str,
                        default='HG00420',
                        help="Test sample for outer CV")
    parser.add_argument('-its',
                        '--innercv_test_sample',
                        type=str,
                        default='HG01053',
                        help="Test sample for inner CV")
    parser.add_argument('-t',
                        '--test_chrom',
                        type=str,
                        default='chr1',
                        help="Chromosome for testing in outerCV")
    parser.add_argument('-v',
                        '--val_chrom',
                        type=str,
                        default='chr22',
                        help="Chromosome for validation in outerCV and innerCV depending on sample")
    parser.add_argument('-l',
                        '--logfile',
                        default='optimize.log',
                        help='File for writing log information')
    parser.add_argument('-e',
                        '--epochs',
                        type=int,
                        default=1,
                        help="Number of epochs")
    parser.add_argument('-n',
                        '--ncalls',
                        type=int,
                        default=12,
                        help="Number of calls of the fitness function")
    parser.add_argument('-b',
                        '--batch_size',
                        type=int,
                        default=32,
                        help="Batch size")
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(format=log_format,
                        filename=os.path.join(args.output, args.logfile),
                        filemode='w',
                        level=logging.INFO)

    t0 = time()

    pos_dict = train(args)

    pos_dict_file = os.path.join(args.output, ''.join([args.test_chrom, '_pos_dict.json']))
    with open(pos_dict_file, 'w') as fp:
        json.dump(pos_dict, fp)

    logging.info('Elapsed time = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
