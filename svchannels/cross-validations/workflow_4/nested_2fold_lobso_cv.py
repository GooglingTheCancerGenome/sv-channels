import os
import argparse
import logging

import matplotlib.pyplot as plt
import pandas as pd
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
from sklearn.metrics import average_precision_score, accuracy_score, balanced_accuracy_score
from tensorflow.keras.utils import to_categorical
from collections import Counter

import sys

# setting path
if '__file__' in vars():
    # print("We are running the script non interactively")
    module_path = os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
    print("Adding {} to sys".format(module_path))
    sys.path.append(module_path)
else:
    # print('We are running the script interactively')
    module_path = "../.."
    print("Adding {} to sys".format(module_path))
    sys.path.append(module_path)

# importing
from model import create_model

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
best_performance = 0.0

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
    logging.info('Labels => bins: {}'.format(Counter(y[:, 7])))

    return X, y


def train_test_val_split_by_sample_and_bin_list(X, y, test_sample, test_bin_list, val_chrom):
    idx_train = np.where(np.logical_and(y[:, 0] != test_sample,
                                        np.logical_and.reduce([y[:, 7] != test_bin
                                                               for test_bin in test_bin_list]),
                                        y[:, 1] != val_chrom))[0]
    idx_test = np.where(np.logical_and(y[:, 0] == test_sample,
                                       np.logical_or.reduce([y[:, 7] == test_bin
                                                             for test_bin in test_bin_list]),
                                       ))[0]
    idx_val = np.where(np.logical_and(y[:, 0] == test_sample, y[:, 1] == val_chrom))[0]

    idx_len_dict = {'train': len(idx_train), 'test': len(idx_test), 'val': len(idx_val)}

    logging.info('Lengths: {}'.format(idx_len_dict))

    return X[idx_train], X[idx_test], X[idx_val], y[idx_train], y[idx_test], y[idx_val]


def train_test_split_by_sample_and_bin(X, y, test_sample, test_bin):
    idx_train = np.where(np.logical_and(y[:, 0] != test_sample,
                                        y[:, 7] != test_bin))[0]
    idx_test = np.where(np.logical_and(y[:, 0] == test_sample, y[:, 7] == test_bin))[0]

    idx_len_dict = {'train': len(idx_train), 'test': len(idx_test)}

    logging.info('Lengths: {}'.format(idx_len_dict))

    return X[idx_train], X[idx_test], y[idx_train], y[idx_test]


def select_by_sample_and_chrom(X, y, sample, chrom):
    idx = np.where(np.logical_and(y[:, 0] == sample,
                                  y[:, 1] == chrom))[0]

    idx_len_dict = {'set': len(idx)}

    logging.info('Lengths: {}'.format(idx_len_dict))

    return X[idx], y[idx]


@use_named_args(dimensions=dimensions)
def fitness(cnn_filters, cnn_layers, cnn_filter_size, fc_nodes,
            dropout_rate, learning_rate, regularization_rate):

    # Save the model if it improves on the best found performance which is stored by the global variable
    # best_performance
    global best_performance
    global best_hyperparameters

    # Print the current combination of hyperparameters
    current_hyperparameters = {
        'cnn_filters': cnn_filters,
        'cnn_layers': cnn_layers,
        'cnn_filter_size': cnn_filter_size,
        'fc_nodes': fc_nodes,
        'dropout_rate': dropout_rate,
        'learning_rate': learning_rate,
        'regularization_rate': regularization_rate
    }
    logging.info(current_hyperparameters)

    logging.info('Inner CV')

    auc_precision_recall = []
    accuracy = []
    balanced_accuracy = []
    validation_auc_pr = []

    # for inner_i, inner_c in enumerate(inner_chr_list):
    # only for test purposes
    for inner_i, inner_c in enumerate([inner_chr_list[0]]):

        logging.info('Considering inner test chromosome {}'.format(inner_c))

        inner_X_train, inner_X_test, \
        inner_y_train, inner_y_test, = train_test_split_by_sample_and_bin(
            outer_X_train, outer_y_train, innercv_test_sample, inner_c
        )

        assert innercv_test_sample not in inner_y_train[:, 0], "Test sample {} in training set: {}".format(
            innercv_test_sample, np.unique(inner_y_train[:, 0])
        )
        assert innercv_test_sample == np.unique(inner_y_test[:, 0]), "Test sample {} not in test set: {}".format(
            innercv_test_sample, np.unique(inner_y_test[:, 0])
        )

        assert inner_c not in inner_y_train[:, 7], "Test bin {} in training set: {}".format(
            inner_c, np.unique(inner_y_train[:, 7])
        )
        assert inner_c == np.unique(inner_y_test[:, 7]), "Test bin {} not in test set: {}".format(
            inner_c, np.unique(inner_y_test[:, 7])
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
                            batch_size=batch_size,
                            shuffle=True,
                            validation_data=(inner_X_val, y_val),
                            class_weight=class_weights_train,
                            verbose=0)

        # Get the predicted probability of testing data
        y_probs = model.predict(inner_X_test)
        y_score = y_probs[:, 0]

        # Data to plot precision - recall curve
        # precision, recall, thresholds = precision_recall_curve(y_test, y_score)
        # Use AUC function to calculate the area under the curve of precision recall curve
        # auc_precision_recall.append(auc(recall, precision))
        # calculate average precision score
        avg_prec = average_precision_score(y_test, y_score)
        auc_precision_recall.append(avg_prec)

        y_predicted = y_probs.argmax(axis=1)
        acc_score = accuracy_score(y_test, y_predicted)
        accuracy.append(acc_score)

        bal_acc_score = balanced_accuracy_score(y_test, y_predicted)
        balanced_accuracy.append(bal_acc_score)

        val_auc = history.history['val_auc'][-1]
        validation_auc_pr.append(val_auc)

        # Delete the model that just finishing training from memory
        del model

        # Clear the Keras session, otherwise it will keep adding new
        # models to the same TensorFlow graph each time we create
        # a model with a different set of hyperparameters.
        tf.keras.backend.clear_session()

    mean_auc = np.mean(auc_precision_recall)
    mean_accuracy = np.mean(accuracy)
    mean_balanced_accuracy = np.mean(balanced_accuracy)
    mean_validation_auc_pr = np.mean(validation_auc_pr)

    logging.info('mean_auc_precision_recall = {}'.format(mean_auc))
    logging.info('mean_accuracy = {}'.format(mean_accuracy))
    logging.info('mean_balanced_accuracy = {}'.format(mean_balanced_accuracy))
    logging.info('mean_validation_auc_pr = {}'.format(mean_validation_auc_pr))

    mean_performance = mean_balanced_accuracy

    # If the AUC of the saved model is greater than the current best performance
    if mean_performance > best_performance:
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
            'mean_auc_precision_recall': mean_auc
        }
        best_hyperparameters = {k: [v] for k, v in best_hyperparameters.items()}

        prefix = '_'.join(['outer', outer_c, 'inner', inner_c, 'cnn-filt', str(cnn_filters),
                           'cnn-lay', str(cnn_layers), 'cnn-filt-size', str(cnn_filter_size),
                           'fc-nodes', str(fc_nodes), 'dropout', str(dropout_rate),
                           'lr', str(learning_rate), 'rr', str(regularization_rate),
                           'pr-auc', str(best_auc)])
        accuracy_plot = os.path.join(output,
                                     ''.join([prefix, '_auc_pr.png']))
        loss_plot = os.path.join(output,
                                 ''.join([prefix, '_loss.png']))

        # print(history.history)
        plt.plot(history.history['auc'])
        plt.plot(history.history['val_auc'])
        plt.title('model auc PR')
        plt.ylabel('AUC PR')
        plt.xlabel('epoch')
        plt.legend(['train', 'val'], loc='upper right')
        plt.tight_layout()
        plt.savefig(accuracy_plot)
        plt.close()

        plt.plot(history.history['loss'])
        plt.plot(history.history['val_loss'])
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'val'], loc='upper left')
        plt.tight_layout()
        plt.savefig(loss_plot)
        plt.close()

        # Update the current greatest AUC score
        best_performance = mean_performance

    # Scikit-optimize does minimization so it tries to
    # find a set the combination of hyperparameters with the lowest fitness value.
    # We want to maximize the model performance so we negate this number.
    return -mean_performance


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


def train(input_args):

    # Store chromosome and position of DELs as keys and difference of posterior probabilities as values
    pos_dict = {}
    model_df = pd.DataFrame()

    global X, y
    X, y = load_data(input_args.input)

    # test sample for the outer cross-validation
    outercv_test_sample = input_args.outercv_test_sample
    # test sample for the inner cross-validation
    global innercv_test_sample
    innercv_test_sample = input_args.innercv_test_sample
    # chromosome for validation
    global val_chrom
    val_chrom = input_args.val_chrom

    # list of bins
    # bins are from 1 to 66
    outer_chr_list = [str(i) for i in np.sort(np.unique(y[:, 7]))]
    # remove bins of the validation chromosome
    val_bins = np.unique(y[np.where(y[:, 3] == val_chrom)[0], 7])
    logging.info('Bins for validation chromosome {} are {}'.format(val_chrom, val_bins))
    outer_chr_list = [b for b in outer_chr_list if b not in val_bins]

    first_list = np.random.choice(outer_chr_list, int(len(outer_chr_list)/2), replace=False)
    second_list = [c for c in outer_chr_list if c not in first_list]
    logging.info("First list:{}\nSecond list:{}".format(first_list, second_list))

    logging.info('Outer CV')
    global outer_i, outer_c

    global best_auc
    best_auc = 0.0

    for outer_i, outer_c_list in enumerate([first_list, second_list]):

        logging.info('Considering outer test chromosome list {}'.format(outer_c_list))

        global outer_X_train, outer_y_train

        outer_X_train, outer_X_test, outer_X_val, \
        outer_y_train, outer_y_test, outer_y_val = train_test_val_split_by_sample_and_bin_list(
            X, y, outercv_test_sample, outer_c_list, val_chrom
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

        for outer_c in outer_c_list:

            assert outer_c not in outer_y_train[:, 7], "Test bin {} in training set: {}".format(
                outer_c, np.unique(outer_y_train[:, 7])
            )
            assert outer_c in np.unique(outer_y_test[:, 7]), "Test bin {} not in test set: {}".format(
                outer_c, np.unique(outer_y_test[:, 7])
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

        global inner_chr_list
        inner_chr_list = [c for c in outer_chr_list if c not in outer_c_list]

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
                            batch_size=batch_size,
                            shuffle=True,
                            validation_data=(outer_X_val, y_val),
                            class_weight=class_weights_train,
                            verbose=0)

        probs = model.predict(outer_X_test, batch_size=100, verbose=False)

        for i in np.arange(outer_y_test.shape[0]):
            chrom1 = outer_y_test[i, 1]
            pos1a = outer_y_test[i, 2]
            pos_dict[chrom1 + '_' + str(int(pos1a) + 1)] = str(probs[i][0])

    model_df.to_csv('model_df.csv')

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
    parser.add_argument('-mi',
                        '--manta_vcf_in',
                        type=str,
                        default='/Users/lsantuari/Documents/Projects/GTCG/sv-channels/'
                                'sv-channels_manuscript/1KG_trios/Manta/HG00420/manta.vcf',
                        help="Manta callset in input for the test sample in the outer CV")
    parser.add_argument('-mo',
                        '--manta_vcf_out',
                        type=str,
                        default=os.path.join(output_dir, 'HG00420.sv-channels.vcf'),
                        help="Manta callset in output for the test sample in the outer CV")
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

    assert os.path.exists(args.manta_vcf_in), f"{args.manta_vcf_in} does not exist!"

    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(format=log_format,
                        filename=args.logfile,
                        filemode='w',
                        level=logging.INFO)

    global n_epochs, batch_size, output

    n_epochs = args.epochs
    batch_size = args.batch_size
    output = args.output

    t0 = time()

    pos_dict = train(args)

    reader = vcf.Reader(open(args.manta_vcf_in, 'r'))
    writer = vcf.Writer(open(args.manta_vcf_out, 'w'), reader)

    for record in reader:
        k = record.CHROM + '_' + str(record.POS)
        if k in pos_dict.keys():
            record.QUAL = str(pos_dict[k])
            writer.write_record(record)

    logging.info('Elapsed time = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
