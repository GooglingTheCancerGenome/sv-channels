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
    for i in range (0, len(X -1)):
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
        chr_list = list(map(str, np.arange(4, 23)))
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
        # print('subsample')

        indices_label = np.where(labels == lab)[0]
        # print(indices_label)
        indices_to_remove = indices_label[np.arange(int(round(len(indices_label) * pc)), len(indices_label))]
        # print(indices_to_remove)
        X = np.delete(data, indices_to_remove, axis=0)
        y = np.delete(labels, indices_to_remove, axis=0)

        # print(X.shape)
        # print(y.shape)
        # print(Counter(y))

        return X, y

    def get_labelled_windows(data_mode):

        logging.info('Loading data...')

        real_training_data, real_training_labels, real_training_id = real_data()
        art_training_data, art_training_labels = artificial_data()

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

        logging.info('Training data shape: %s' % str(training_data.shape))
        logging.info('Training labels shape: %s' % str(training_labels.shape))

        for l in ['UK', 'INS_pos']:
            logging.info('Removing label %s' % l)
            training_data, training_labels = remove_label(training_data, training_labels, label=l)

        logging.info('Data loaded.')

        return training_data, training_labels

    logging.info('Running with mode ' + data_mode + '...')

    metrics = dict()

    #for pc in np.linspace(0.1, 1, num=10):
    for pc in [0.1]:

        # print(pc)
        logging.info('Running with proportion ' + str(pc) + '...')

        pc_str = str(round(pc, 1))
        metrics[pc_str] = dict()
        training_data, training_labels = get_labelled_windows(data_mode)

        X, y = subsample_nosv(training_data, training_labels, pc, 'noSV')

        del training_data
        del training_labels

        X = transpose_dataset(X)

        logging.info('X shape: %s' % str(X.shape))
        logging.info('y shape: %s' % str(y.shape))

        mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}
        y_num = np.array([mapclasses[c] for c in y], dtype='int')
        y_binary = to_categorical(y_num)

        intermediate_result, metrics[pc_str] = cross_validation(X, y, y_binary,
                                                                X_test, y_test, y_test_binary,
                                                                channel_set, proportion=round(pc, 1),
                                                                data_mode=data_mode)
        logging.info(intermediate_result)
        results.to_csv(filename + '_' + data_mode + '_' + str(round(pc, 1)) + file_extension, sep='\t')
        results = results.append(intermediate_result)
        del X, y

    logging.info('Writing metrics...')

    metrics_output_file = filename + '_metrics_' + data_mode + '.pickle.gz'
    with gzip.GzipFile(metrics_output_file, "wb") as f:
        pickle.dump(f, metrics)
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

    # i = 0
    # for model, params, model_types in models:
    #     print('model ' + str(i))
    #     i = i + 1
    #     print(params)
    #     model.summary()

    return models


def cross_validation(X, y, y_binary, X_hold_out_test,
                     y_hold_out_test, y_hold_out_test_binary,
                     channels, proportion, data_mode):
    results = pd.DataFrame()

    # From https://medium.com/@literallywords/stratified-k-fold-with-keras-e57c487b1416
    kfold_splits = 2

    metrics = dict()

    # Instantiate the cross validator
    skf = StratifiedKFold(n_splits=kfold_splits, shuffle=True)

    # Loop through the indices the split() method returns
    for index, (train_indices, test_indices) in enumerate(skf.split(X, y)):
        print("Training on fold " + str(index + 1) + "/10...")

        # Generate batches from indices
        xtrain, xtest = X[train_indices], X[test_indices]
        # ytrain, ytest = y[train_indices], y[test_indices]
        ytrain_binary, ytest_binary = y_binary[train_indices], y_binary[test_indices]

        # split into train/validation sets
        xtrain, xval, ytrain_binary, yval = train_test_split(xtrain, ytrain_binary,
                                                             test_size=0.2, random_state=2)

        # Create a new model
        model = create_model(X, y_binary)

        # Debug message I guess
        print("Training new iteration on " + str(xtrain.shape[0]) + " training samples, " +
              str(xval.shape[0]) + " validation samples, this may take a while...")

        history, model = train_model(model, xtrain, ytrain_binary, xval, yval)

        accuracy_history = history.history['acc']
        val_accuracy_history = history.history['val_acc']
        print("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
            val_accuracy_history[-1]))

        score_test = model.evaluate(xtest, ytest_binary, verbose=False)
        print('Test loss and accuracy of best model: ' + str(score_test))

        results, metrics[str(index + 1)] = evaluate_model(model, X_hold_out_test, y_hold_out_test,
                                                          y_hold_out_test_binary, results, index, channels,
                                                          proportion, data_mode,
                                                          train_set_size=xtrain.shape[0],
                                                          validation_set_size=xval.shape[0]
                                                          )
        # evaluate_model(model, X_test, y_test_binary, results, index, channels)

    return results, metrics


def train_model(model, xtrain, ytrain, xval, yval):

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
                             verbose=False)

    return history, best_model


def evaluate_model(model, X_test, y_test, ytest_binary, results, cv_iter, channels, proportion, data_mode,
                   train_set_size, validation_set_size):
    mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}
    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    # print(dict_sorted)
    # class_labels = [i[0] for i in dict_sorted]

    n_classes = ytest_binary.shape[1]
    # print(y_binarized)
    # print(n_classes)

    probs = model.predict_proba(X_test, batch_size=1000, verbose=False)

    # For each class
    precision = dict()
    recall = dict()
    thresholds = dict()
    average_precision = dict()

    # for i in range(n_classes):
    for k, i in mapclasses.items():

        precision[k], recall[k], thresholds[k] = precision_recall_curve(ytest_binary[:, i],
                                                                        probs[:, i])
        average_precision[k] = average_precision_score(ytest_binary[:, i], probs[:, i])

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(ytest_binary.ravel(),
                                                                    probs.ravel())
    average_precision["micro"] = average_precision_score(ytest_binary, probs,
                                                         average="micro")
    average_precision["micro"] = average_precision_score(ytest_binary, probs,
                                                         average="micro")
    print('Average precision score, micro-averaged over all classes: {0:0.2f}'
          .format(average_precision["micro"]))

    results = results.append({
        "channels": channels,
        "data_mode": data_mode,
        "proportion": proportion,
        "fold": cv_iter + 1,
        "training_set_size": train_set_size,
        "validation_set_size": validation_set_size,
        "test_set_size": X_test.shape[0],
        "average_precision_score": average_precision["micro"]
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
    #     plt.savefig('Plots/Precision_Recall_multiclass_Iter_'+str(cv_iter)+'_'+channels+'.png', bbox_inches='tight')

    return results, (precision, recall, thresholds)


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
