# Imports
import gzip
import os
import argparse
import logging
import pickle

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

HPC_MODE = True

sample_name = 'NA12878'
date = '270219'
label_type = 'Mills2011_nanosv'
datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari' if HPC_MODE else '/Users/lsantuari/Documents'
datapath_training = datapath_prefix + '/Processed/Test/' + \
                    date + '/TestData_' + date + '/' + sample_name + '/TrainingData/'
datapath_test = datapath_prefix + '/Processed/Test/' + \
                date + '/TestData_' + date + '/' + sample_name + '/TestData/'


def create_dir(directory):
    '''
    Create a directory if it does not exist. Raises an exception if the directory exists.
    :param directory: directory to create
    :return: None
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def get_channel_labels():
    # Fill labels for legend

    labels = list()
    labels.append("coverage")
    labels.append("mean_read_quality")
    labels.append("#left_clipped_reads")
    labels.append("#right_clipped_reads")
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
                labels.append(direction + '_' + clipped + '_Clipped_' + value)

    labels.append("#left split reads")
    labels.append("#right split reads")

    for clipped in ['L', 'R']:
        for value in ['median']:
            labels.append(clipped + '_SplitRead_' + value)

    labels.append("Mappability")

    for nuc in ['A', 'T', 'C', 'G', 'N']:
        labels.append("One_hot_encoding_" + nuc)

    for k, l in enumerate(labels):
        logging.info(str(k) + ':' + l)

    return labels


def transpose_dataset(X):
    image = []
    for i in range(0, len(X - 1)):
        tr = X[i].transpose()
        image.append(tr)
    return np.array(image)


def data(datapath, channels):
    data_output_file = datapath + sample_name + '_' + label_type + '_channels.npy.gz'
    with gzip.GzipFile(data_output_file, "rb") as f:
        X = np.load(f)

    label_output_file = datapath + sample_name + '_' + label_type + '_labels.npy.gz'
    with gzip.GzipFile(label_output_file, "rb") as f:
        y = np.load(f)
    f.close()

    label_output_file = datapath + sample_name + '_' + label_type + '_labels_binary.npy.gz'
    with gzip.GzipFile(label_output_file, "rb") as f:
        y_binary = np.load(f)
    f.close()

    id_output_file = datapath + sample_name + '_' + label_type + '_ids.npy.gz'
    with gzip.GzipFile(id_output_file, "rb") as f:
        win_ids = np.load(f)
    f.close()

    return X[:, :, channels], y, y_binary, win_ids


def real_data():
    # Create reference test set using Chr4 to ChrX

    data_output_file = os.path.join(datapath_training, sample_name + '_' + label_type + '_training.npz')

    def get_label_dict():
        # Load label dictionary
        dico_file = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                                 date, 'TestData_' + date, sample_name, 'MultiLabelData/labels.pickle.gz')
        with gzip.GzipFile(dico_file, "rb") as f:
            dico = np.load(f)
        f.close()

        return dico

    def save_data(training_data, training_labels, training_id):

        np.savez(data_output_file, training_data=training_data,
                 training_labels=training_labels, training_id=training_id)
        os.system('gzip -f ' + data_output_file)

    if not os.path.exists(data_output_file + '.gz'):

        dico = get_label_dict()

        # Leaving out chromosome Y and MT for the moment
        # chr_list = list(map(str, np.arange(4, 23)))
        chr_list = list(map(str, np.arange(1, 23)))
        chr_list.append('X')

        logging.info(chr_list)

        training_data = []
        training_labels = []
        training_id = []

        datapath = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                                date, 'TestData_' + date, sample_name, 'ChannelData')

        for i in chr_list:
            logging.info('Loading data for Chr%s' % i)
            data_file = os.path.join(datapath, sample_name + '_' + str(i) + '.npy.gz')
            with gzip.GzipFile(data_file, "rb") as f:
                data_mat = np.load(f)
                training_data.extend(data_mat)
            f.close()

            training_labels.extend(dico[label_type][i])
            training_id.extend([d['chromosome'] + '_' + str(d['position']) for d in dico['id'][i]])

        logging.info(Counter(training_labels))

        training_data = np.array(training_data)
        training_labels = np.array(training_labels)
        training_id = np.array(training_id)

        assert len(training_data) == len(training_labels)

        save_data(training_data, training_labels, training_id)

    else:

        logging.info('Loading real data...')

        with gzip.GzipFile(data_output_file + '.gz', 'rb') as f:

            npzfiles = np.load(f)
            training_data = npzfiles['training_data']
            training_labels = npzfiles['training_labels']
            training_id = npzfiles['training_id']

    return training_data, training_labels, training_id


def artificial_data():
    training_data = []
    training_labels = []

    base_dir = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                            date, 'TrainingData_' + date)
    sample = 'G1'

    for svtype in ['INDEL', 'INDEL_HOM']:
        datapath = os.path.join(base_dir, svtype, sample)
        data_file = os.path.join(datapath, 'ChannelData', sample + '.npy.gz')
        label_file = os.path.join(datapath, 'LabelData', sample + '_17_label.npy.gz')

        with gzip.GzipFile(data_file, "rb") as f:
            data_mat = np.load(f)
            training_data.extend(data_mat)
        f.close()

        with gzip.GzipFile(label_file, "rb") as f:
            label_list = np.load(f)
            training_labels.extend(label_list)
        f.close()

    training_data = np.array(training_data)
    training_labels = np.array(training_labels)

    assert len(training_data) == len(training_labels)

    return training_data, training_labels


def mixed_data(output, data_mode):
    # create_dir('Plots')

    filename, file_extension = os.path.splitext(output)

    results = pd.DataFrame()

    labels = get_channel_labels()

    channels = np.arange(len(labels))
    channel_set = 'all'

    logging.info('Running cv with channels ' + channel_set + ':')
    for i in channels:
        logging.info(str(i) + ':' + labels[i])

    # Load the test data
    logging.info('Loading hold-out test set')
    X_test, y_test, y_test_binary, win_ids_test = data(datapath_test, channels)

    def subsample_nosv(data, labels, pc, lab):

        logging.info('subsample noSV:')

        indices_label = np.where(labels == lab)[0]
        # print(indices_label)
        indices_to_remove = indices_label[np.arange(int(round(len(indices_label) * pc)), len(indices_label))]
        # print(indices_to_remove)
        X = np.delete(data, indices_to_remove, axis=0)
        y = np.delete(labels, indices_to_remove, axis=0)

        logging.info('X shape: %s' % str(X.shape))
        logging.info('y shape: %s' % str(y.shape))
        logging.info('y labels: %s' % str(Counter(y)))

        return X, y

    def get_labelled_windows(data_mode):

        logging.info('Loading data...')

        real_training_data, real_training_labels, real_training_id = real_data()

        logging.info('Real data shape: %s' % str(real_training_data.shape))
        logging.info('Real labels shape: %s' % str(real_training_labels.shape))
        logging.info('Real labels: %s' % str(Counter(real_training_labels)))

        art_training_data, art_training_labels = artificial_data()

        logging.info('Artificial data shape: %s' % str(art_training_data.shape))
        logging.info('Artificial labels shape: %s' % str(art_training_labels.shape))
        logging.info('Artificial labels: %s' % str(Counter(art_training_labels)))

        # artificial data only 0
        # real data only 1
        # mixed data 2
        # for data_mode in ['artificial','real','mixed']:

        if data_mode == 'artificial':

            # artificial data only 0
            indices_label = np.where(real_training_labels == 'noSV')[0]
            training_data = np.concatenate((art_training_data,
                                            real_training_data[indices_label]), axis=0)
            training_labels = np.concatenate((art_training_labels,
                                              real_training_labels[indices_label]), axis=0)
        elif data_mode == 'real':

            # real data only 1
            training_data = real_training_data
            training_labels = real_training_labels

        elif data_mode == 'mixed':

            # mixed data 2
            training_data = np.concatenate((real_training_data, art_training_data), axis=0)
            training_labels = np.concatenate((real_training_labels, art_training_labels), axis=0)

        for l in ['UK', 'INS_pos']:
            logging.info('Removing label %s' % l)
            training_data, training_labels = remove_label(training_data, training_labels, label=l)

        training_data = transpose_dataset(training_data)
        logging.info('Training data shape: %s' % str(training_data.shape))
        logging.info('Training labels shape: %s' % str(training_labels.shape))
        logging.info('Training labels: %s' % str(Counter(training_labels)))

        logging.info('Data loaded.')

        return training_data, training_labels

    logging.info('Running with mode ' + data_mode + '...')

    metrics = dict()

    for pc in np.linspace(0.1, 1, num=10):
        # for pc in [0.1]:

        windows, labels = get_labelled_windows(data_mode)

        # print(pc)
        logging.info('Running with proportion ' + str(pc) + '...')

        pc_str = str(round(pc, 1))
        metrics[pc_str] = dict()

        X, y = subsample_nosv(windows, labels, pc, 'noSV')

        del windows
        del labels

        logging.info('X shape: %s' % str(X.shape))
        logging.info('y shape: %s' % str(y.shape))

        mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}
        y_num = np.array([mapclasses[c] for c in y], dtype='int')
        y_binary = to_categorical(y_num)

        intermediate_result, metrics[pc_str] = cross_validation(X, y, y_binary,
                                                                X_test, y_test, y_test_binary,
                                                                channel_set, proportion=round(pc, 1),
                                                                data_mode=data_mode,
                                                                output=filename)
        logging.info(intermediate_result)
        intermediate_result.to_csv(
            filename + '_' + data_mode + '_' + pc_str + file_extension, sep='\t')
        results = results.append(intermediate_result)

        del X, y

    logging.info('Writing metrics...')

    metrics_output_file = filename + '_metrics_' + data_mode + '.pickle.gz'
    with gzip.GzipFile(metrics_output_file, "wb") as f:
        pickle.dump(metrics, f)
    f.close()

    logging.info(results)

    results.to_csv(filename + '_' + data_mode + file_extension, sep='\t')


def remove_label_with_id(training_data, training_labels, training_id, label='UK'):
    # Remove windows labelled as label
    keep = np.where(np.array(training_labels) != label)
    training_data = training_data[keep]
    training_labels = training_labels[keep]
    training_id = training_id[keep]

    return training_data, training_labels, training_id


def remove_label(training_data, training_labels, label='UK'):
    # Remove windows labelled as label
    keep = np.where(np.array(training_labels) != label)
    training_data = training_data[keep]
    training_labels = training_labels[keep]

    return training_data, training_labels


def load_data(datapath, channels):
    data_output_file = datapath + sample_name + '_' + label_type + '_channels.npy.gz'
    with gzip.GzipFile(data_output_file, "rb") as f:
        X = np.load(f)

    label_output_file = datapath + sample_name + '_' + label_type + '_labels.npy.gz'
    with gzip.GzipFile(label_output_file, "rb") as f:
        y = np.load(f)
    f.close()

    label_output_file = datapath + sample_name + '_' + label_type + '_labels_binary.npy.gz'
    with gzip.GzipFile(label_output_file, "rb") as f:
        y_binary = np.load(f)
    f.close()

    id_output_file = datapath + sample_name + '_' + label_type + '_ids.npy.gz'
    with gzip.GzipFile(id_output_file, "rb") as f:
        win_ids = np.load(f)
    f.close()

    return X[:, :, channels], y, y_binary, win_ids


def create_model(X, y_binary):

    models = modelgen.generate_models(X.shape,
                                      y_binary.shape[1],
                                      number_of_models=1,
                                      model_type='CNN',
                                      cnn_min_layers=2,
                                      cnn_max_layers=2,
                                      cnn_min_filters=4,
                                      cnn_max_filters=4,
                                      cnn_min_fc_nodes=6,
                                      cnn_max_fc_nodes=6,
                                      low_lr=2, high_lr=2,
                                      low_reg=1, high_reg=1,
                                      kernel_size=7)

    i = 0
    for model, params, model_types in models:
        logging.info('model ' + str(i))
        i = i + 1
        logging.info(params)
        logging.info(model.summary())

    return models


def cross_validation(X, y, y_binary,
                     X_hold_out_test, y_hold_out_test, y_hold_out_test_binary,
                     channels, proportion, data_mode, output):
    results = pd.DataFrame()

    # From https://medium.com/@literallywords/stratified-k-fold-with-keras-e57c487b1416
    kfold_splits = 10

    metrics = dict()

    # Instantiate the cross validator
    skf = StratifiedKFold(n_splits=kfold_splits, shuffle=True)

    # Loop through the indices the split() method returns
    for index, (train_indices, test_indices) in enumerate(skf.split(X, y)):

        print("Training on fold " + str(index + 1) + "/10...")

        # Generate batches from indices
        xtrain, xtest = X[train_indices], X[test_indices]
        ytrain, ytest = y[train_indices], y[test_indices]
        ytrain_binary, ytest_binary = y_binary[train_indices], y_binary[test_indices]

        # split into train/validation sets
        xtrain_split, xval_split, ytrain_split, yval_split = train_test_split(xtrain, ytrain,
                                                                              test_size=0.2, random_state=2,
                                                                              stratify=ytrain)

        logging.info('Training data shape: %s' % str(xtrain_split.shape))
        logging.info('Training labels shape: %s' % str(ytrain_split.shape))
        logging.info('Training labels: %s' % str(Counter(ytrain_split)))

        logging.info('Validation data shape: %s' % str(xval_split.shape))
        logging.info('Validation labels shape: %s' % str(yval_split.shape))
        logging.info('Validation labels: %s' % str(Counter(yval_split)))

        logging.info('Test data shape: %s' % str(xtest.shape))
        logging.info('Test labels shape: %s' % str(ytest.shape))
        logging.info('Test labels: %s' % str(Counter(ytest)))

        mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}

        ytrain_split_num = np.array([mapclasses[c] for c in ytrain_split], dtype='int')
        ytrain_split_binary = to_categorical(ytrain_split_num)

        yval_split_num = np.array([mapclasses[c] for c in yval_split], dtype='int')
        yval_split_binary = to_categorical(yval_split_num)

        # Create a new model
        model = create_model(xtrain_split, ytrain_split_binary)

        # Debug message I guess
        print("Training new iteration on " + str(xtrain_split.shape[0]) + " training samples, " +
              str(xval_split.shape[0]) + " validation samples, this may take a while...")

        class_weights = class_weight.compute_class_weight('balanced',
                                                          np.unique(ytrain_split),
                                                          ytrain_split)
        class_weight_dict = dict(enumerate(class_weights))

        history, model = train_model(model, xtrain_split, ytrain_split_binary, class_weight_dict,
                                     xval_split, yval_split_binary)

        accuracy_history = history.history['acc']
        val_accuracy_history = history.history['val_acc']
        print("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
            val_accuracy_history[-1]))

        score_test = model.evaluate(xtest, ytest_binary, verbose=False)
        print('Test loss and accuracy of best model: ' + str(score_test))

        results, metrics[str(index + 1)] = evaluate_model(model, xtest, ytest,
                                                          ytest_binary, results, index, channels,
                                                          proportion, data_mode, output,
                                                          train_set_size=xtrain_split.shape[0],
                                                          validation_set_size=xval_split.shape[0])
        # evaluate_model(model, X_test, y_test_binary, results, index, channels)

    return results, metrics


def train_model(model, xtrain, class_weights, ytrain, xval, yval):

    train_set_size = xtrain.shape[0]
    nr_epochs = 1

    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(xtrain, ytrain,
                                                                                      xval, yval,
                                                                                      model, nr_epochs=nr_epochs,
                                                                                      subset_size=train_set_size,
                                                                                      verbose=False)

    best_model_index = np.argmax(val_accuracies)
    best_model, best_params, best_model_types = model[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    history = best_model.fit(xtrain, ytrain,
                             epochs=nr_epochs, validation_data=(xval, yval),
                             verbose=False,
                             class_weight=class_weights)

    return history, best_model


def evaluate_model(model, X_test, y_test, ytest_binary, results, cv_iter, channels, proportion, data_mode, output,
                   train_set_size, validation_set_size):

    mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}

    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    # print(dict_sorted)
    class_labels = [i[0] for i in dict_sorted]

    n_classes = ytest_binary.shape[1]
    # print(y_binarized)
    # print(n_classes)

    probs = model.predict_proba(X_test, batch_size=1000, verbose=False)

    # columns are predicted, rows are truth
    predicted = probs.argmax(axis=1)
    # print(predicted)
    y_index = ytest_binary.argmax(axis=1)

    # print(y_index)
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [class_labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [class_labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in class_labels], fill_value=0)
    confusion_matrix.to_csv(output + '_confusion_matrix_' + data_mode +
                            '_' + str(proportion) + '_' + str(cv_iter + 1) + '.csv', sep='\t')

    # For each class
    precision = dict()
    recall = dict()
    f1_score_metric = dict()
    thresholds = dict()
    average_precision = dict()

    # for i in range(n_classes):
    for k, i in mapclasses.items():
        precision[k], recall[k], thresholds[k] = precision_recall_curve(ytest_binary[:, i],
                                                                        probs[:, i])
        average_precision[k] = average_precision_score(ytest_binary[:, i], probs[:, i], average="weighted")
        f1_score_metric[k] = f1_score(y_index, predicted, average=None)[i]

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(ytest_binary.ravel(),
                                                                    probs.ravel())

    average_precision["weighted"] = average_precision_score(ytest_binary, probs, average="weighted")
    print('Average precision score, weighted over all classes: {0:0.2f}'
          .format(average_precision["weighted"]))

    f1_score_metric["weighted"] = f1_score(y_index, predicted, average="weighted")

    results = results.append({
        "channels": channels,
        "data_mode": data_mode,
        "proportion": proportion,
        "fold": cv_iter + 1,
        "training_set_size": train_set_size,
        "validation_set_size": validation_set_size,
        "test_set_size": X_test.shape[0],
        "average_precision_score": average_precision["weighted"],
        "f1_score": f1_score_metric["weighted"]
    }, ignore_index=True)

    # for iter_class in mapclasses.values():
    #
    #     predicted = probs.argmax(axis=1)
    #     #print(predicted)
    #     y_pred_class = np.array([1 if i == iter_class else 0 for i in predicted])
    #     #print(y_pred_class)
    #
    #     # keep probabilities for the positive outcome only
    #     probs_class = probs[:, iter_class]
    #     #print(probs_class)
    #
    #     #print(y_test)
    #
    #     y_test_class = np.array([1 if i[iter_class] == 1 else 0 for i in ytest_binary])
    #
    #     # calculate precision-recall curve
    #     precision, recall, thresholds = precision_recall_curve(y_test_class, probs_class)
    #     # calculate F1 score
    #     f1 = f1_score(y_test_class, y_pred_class)
    #     # calculate precision-recall AUC
    #     auc_value = auc(recall, precision)
    #     # calculate average precision score
    #     ap = average_precision_score(y_test_class, probs_class)
    #     print('f1=%.3f auc=%.3f average_precision_score=%.3f' % (f1, auc_value , ap))
    #     # plot no skill
    #     plt.plot([0, 1], [0.5, 0.5], linestyle='--')
    #     # plot the roc curve for the model
    #     plt.plot(recall, precision, marker='.')
    #     # show the plot
    #     plt.savefig('PrecRec_' + str(cv_iter) +
    #                 '_'+proportion+'_'+str(cv_iter + 1)+'.png', bbox_inches='tight')

    plot_precision_recall(data_mode, proportion, cv_iter, mapclasses,
                          precision, recall, average_precision, output)

    return results, (average_precision, precision, recall, thresholds, f1_score_metric)


def run_cv():
    labels = get_channel_labels()

    results = pd.DataFrame()

    # # for channel_index in np.arange(0,len(labels)):
    # for channel_set, channels in channel_list.items():

    # channels.append(channel_index)

    channels = np.arange(len(labels))
    channel_set = 'all'

    logging.info('Running cv with channels ' + channel_set + ':')
    # for i in channels:
    #     print(str(i) + ':' + labels[i])

    # Load the data
    X, y, y_binary, win_ids = data(datapath_training, channels)
    X_test, y_test, y_test_binary, win_ids_test = data(datapath_test, channels)

    results = results.append(cross_validation(X, y, y_binary, X_test, y_test,
                                              y_test_binary, channel_set, 1.0))

    logging.info(results)
    results.to_csv("NA12878/CV_results.csv", sep='\t')


def plot_results():
    source = pd.read_csv(filepath_or_buffer='NA12878/CV_results.csv', delimiter='\t')

    import numpy as np
    import matplotlib.pyplot as plt

    means = source.groupby('channels')['average_precision_score'].agg(np.mean).sort_values()
    print('average_precision_score mean: %d' % means)
    std = source.groupby('channels')['average_precision_score'].agg(np.std)
    print('average_precision_score mean: %d' % std)
    ind = np.arange(len(list(means.index)))  # the x locations for the groups
    width = 0.50  # the width of the bars: can also be len(x) sequence

    plt.bar(ind, means, width, yerr=std)

    plt.ylabel('average_precision_score')
    plt.title('average_precision_score per channel set')
    plt.xticks(ind, list(means.index))
    plt.xticks(rotation=45, horizontalalignment='right')
    # plt.yticks(np.arange(0.8, 1))
    # plt.legend((p1[0]), ('Bar'))
    plt.ylim(bottom=0.8)
    plt.tight_layout()

    plt.savefig('Results.png', bbox_inches='tight')
    plt.close()


def plot_precision_recall(data_mode, proportion, cv_iter,
                          mapclasses, precision, recall, average_precision, output):
    from itertools import cycle
    # setup plot details
    colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])

    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

    lines.append(l)
    labels.append('iso-f1 curves')
    l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
    lines.append(l)
    labels.append('weighted-average Precision-recall (area = {0:0.2f})'
                  ''.format(average_precision["weighted"]))

    for i, color in zip(mapclasses.keys(), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                      ''.format(i, average_precision[i]))

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))

    plt.savefig(output + '_PrecRec_' + data_mode +
                '_' + str(proportion) + '_' + str(cv_iter + 1) + '.png', bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Tests multiple artificial/real/mixed training sets')
    parser.add_argument('-m', '--mode', default='artificial',
                        help='Data mode: artificial/real/mixed')
    parser.add_argument('-o', '--output', default='CV_results.csv',
                        help='File in which to write output.')
    parser.add_argument('-l', '--logfile', default='CV_results.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    filename, file_extension = os.path.splitext(args.logfile)
    logfilename = filename + '_' + args.mode + file_extension

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    mixed_data(output=args.output, data_mode=args.mode)


if __name__ == '__main__':
    main()
