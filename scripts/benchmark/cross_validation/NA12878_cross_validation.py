# Imports
import gzip
import os
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

HPC_MODE = True
sample_name = 'NA12878'
date = '270219'
label_type = 'Mills2011_nanosv'
datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari' if HPC_MODE else '/Users/lsantuari/Documents'
datapath_training = datapath_prefix + '/Processed/Test/' + \
                    date + '/TestData_' + date + '/' + sample_name + '/TrainingData/'
datapath_test = datapath_prefix + '/Processed/Test/' + \
                date + '/TestData_' + date + '/' + sample_name + '/TestData/'

balancing = 'balanced'


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
        print(str(k) + ':' + l)

    return labels


def data(datapath, channels):

    if datapath == datapath_test:

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

    elif datapath == datapath_training:

        data_output_file = os.path.join(datapath, '_'.join([sample_name, label_type, balancing]))

        with gzip.GzipFile(data_output_file + '.npz.gz', 'rb') as f:
            npzfiles = np.load(f)
            X = npzfiles['X']
            y = npzfiles['y']
            y_binary = npzfiles['y_binary']
            win_ids = npzfiles['z']

    # print(X.shape)
    # print(y.shape)
    # print(y.shape)
    # print(win_ids.shape)

    # idx = np.arange(0,9)
    # idx = np.append(idx, np.arange(33,35))
    # idx = np.append(idx, np.arange(41, 44))
    # idx = np.append(idx,[12,16,20,24,28,32])

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


#def cross_validation(X, y, y_binary, X_hold_out_test, y_hold_out_test, y_hold_out_test_binary, channels):
def cross_validation(X, y, y_binary, win_ids, channels):

    results = pd.DataFrame()

    # From https://medium.com/@literallywords/stratified-k-fold-with-keras-e57c487b1416
    kfold_splits = 10

    # Instantiate the cross validator
    skf = StratifiedKFold(n_splits=kfold_splits, shuffle=True)

    # Loop through the indices the split() method returns
    for index, (train_indices, test_indices) in enumerate(skf.split(X, y)):
        print("Training on fold " + str(index + 1) + "/10...")

        # Generate batches from indices
        xtrain, xtest = X[train_indices], X[test_indices]
        ytrain, ytest = y[train_indices], y[test_indices]
        ytrain_binary, ytest_binary = y_binary[train_indices], y_binary[test_indices]
        win_ids_test = win_ids[test_indices]

        # split into train/validation sets
        xtrain_split, xval_split, ytrain_split, yval_split = train_test_split(xtrain, ytrain,
                                                                              test_size=0.2, random_state=2,
                                                                              stratify=ytrain,
                                                                              shuffle=True)

        print('Training data shape: %s' % str(xtrain_split.shape))
        print('Training labels shape: %s' % str(ytrain_split.shape))
        print('Training labels: %s' % str(Counter(ytrain_split)))

        print('Validation data shape: %s' % str(xval_split.shape))
        print('Validation labels shape: %s' % str(yval_split.shape))
        print('Validation labels: %s' % str(Counter(yval_split)))

        print('Test data shape: %s' % str(xtest.shape))
        print('Test labels shape: %s' % str(ytest.shape))
        print('Test labels: %s' % str(Counter(ytest)))

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
        print("Class weights: %s" % str(class_weight_dict))

        sample_weights = np.array([class_weight_dict[mapclasses[c]] for c in ytrain_split])
        print("Sample weights: %s" % Counter(sample_weights))

        history, model = train_model(model, xtrain_split, ytrain_split_binary,
                                     xval_split, yval_split_binary, sample_weights)

        accuracy_history = history.history['acc']
        val_accuracy_history = history.history['val_acc']
        print("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
            val_accuracy_history[-1]))

        score_test = model.evaluate(xtest, ytest_binary, verbose=False)
        print('Test loss and accuracy of best model: ' + str(score_test))

        results = evaluate_model(model, xtest, ytest, ytest_binary, win_ids_test, results, index, channels,
                                 train_set_size=xtrain_split.shape[0],
                                 validation_set_size=xval_split.shape[0]
                                 )
        # evaluate_model(model, X_test, y_test_binary, results, index, channels)

    return results


def train_model(model, xtrain, ytrain, xval, yval, sample_weights):

    nr_epochs = 1

    train_set_size = xtrain.shape[0]

    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(xtrain, ytrain,
                                                                                      xval, yval,
                                                                                      model, nr_epochs=1,
                                                                                      subset_size=train_set_size,
                                                                                      verbose=False)

    best_model_index = np.argmax(val_accuracies)
    best_model, best_params, best_model_types = model[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    history = best_model.fit(xtrain, ytrain,
                             epochs=nr_epochs, validation_data=(xval, yval),
                             verbose=False,
                             # sample_weight=sample_weights,
                             shuffle=True)

    return history, best_model


def evaluate_model(model, X_test, y_test, ytest_binary, win_ids_test, results, cv_iter, channels,
                   train_set_size, validation_set_size):

    def write_bed(win_ids, pred_class, probs):
        with open('NA12878/' + 'predictions_'+ str(cv_iter + 1) + '.bed','w') as f:
            for i, c, p in zip(win_ids, pred_class, probs):
                line ='\t'.join([i['chromosome'], i['position'],
                                 i['position']+1, c, p])
                f.write(line+'\n')

    mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}
    # mapclasses = {'DEL': 0, 'noDEL': 1}

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

    predicted_class = [class_labels[i] for i in predicted]
    probs_class = [probs[i] for i in predicted]
    write_bed(win_ids_test, predicted_class, probs_class)

    # print(y_index)
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [class_labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [class_labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in class_labels], fill_value=0)
    confusion_matrix.to_csv('NA12878/' + 'NA12878_confusion_matrix' +
                            '_' + str(cv_iter + 1) + '.csv', sep='\t')

    # For each class
    precision = dict()
    recall = dict()
    average_precision = dict()
    for i in range(n_classes):
        precision[i], recall[i], _ = precision_recall_curve(ytest_binary[:, i],
                                                            probs[:, i])
        average_precision[i] = average_precision_score(ytest_binary[:, i], probs[:, i])

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
        "fold": cv_iter + 1,
        "training_set_size": train_set_size,
        "validation_set_size": validation_set_size,
        "test_set_size": X_test.shape[0],
        "average_precision_score": average_precision["micro"]
    }, ignore_index=True)

    plt.figure()
    plt.step(recall['micro'], precision['micro'], color='b', alpha=0.2,
             where='post')
    plt.fill_between(recall["micro"], precision["micro"], alpha=0.2, color='b')

    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
    plt.title(
        'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
            .format(average_precision["micro"]))

    create_dir('NA12878/Plots')
    plt.savefig('NA12878/Plots/Precision_Recall_avg_prec_score_Iter_' + str(cv_iter) +
                '_' + channels + '_' + balancing + '.png',
                bbox_inches='tight')
    plt.close()

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
    labels.append('micro-average Precision-recall (area = {0:0.2f})'
                  ''.format(average_precision["micro"]))

    for i, color in zip(range(n_classes), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                      ''.format(class_labels[i], average_precision[i]))

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))

    plt.savefig('NA12878/Plots/Precision_Recall_avg_prec_score_per_class_Iter_' +
                str(cv_iter) + '_' + channels + '_' + balancing + '.png', bbox_inches='tight')
    plt.close()

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

    return results


def run_cv():

    labels = get_channel_labels()

    create_dir('NA12878')
    create_dir('NA12878/Plots')

    # basic_channels = np.append(np.arange(0, 9), [33, 34])
    # channel_list = {"base": basic_channels,
    #                 "base_PE_allReads": np.append(basic_channels, [19, 31]),
    #                 "base_PE_outliers": np.append(basic_channels, [12, 16, 20, 24, 28, 32]),
    #                 "base_PE_allReads_outliers": np.append(basic_channels, [19, 31, 12, 16, 20, 24, 28, 32]),
    #                 "base_PE_allReads_outliers_splitReads":
    #                     np.append(basic_channels, [19, 31, 12, 16, 20, 24, 28, 32, 37, 40]),
    #                 "base_PE_allReads_outliers_splitReads_GC":
    #                     np.append(basic_channels, [19, 31, 12, 16, 20, 24, 28, 32, 37, 40, 41]),
    #                 "base_PE_allReads_outliers_splitReads_Mappability":
    #                     np.append(basic_channels, [19, 31, 12, 16, 20, 24, 28, 32, 37, 40, 42]),
    #                 "base_PE_allReads_outliers_splitReads_OneHot":
    #                     np.append(basic_channels, [19, 31, 12, 16, 20, 24, 28, 32, 37, 40, 43]),
    #                 }
    # print(channel_list)

    results = pd.DataFrame()

    # # for channel_index in np.arange(0,len(labels)):
    # for channel_set, channels in channel_list.items():

    # channels.append(channel_index)

    channels = np.arange(len(labels))
    channel_set = 'all'

    print('Running cv with channels ' + channel_set + ':')
    # for i in channels:
    #     print(str(i) + ':' + labels[i])

    # Load the data
    X, y, y_binary, win_ids = data(datapath_training, channels)
    #X_test, y_test, y_test_binary, win_ids_test = data(datapath_test, channels)

    def oversample(X, y):

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

        print(Counter(labels_balanced))

        X = np.array(data_balanced)
        y = np.array(labels_balanced)

        return X, y

    print('X shape: %s' % str(X.shape))
    print('y shape: %s' % str(y.shape))

    X, y = oversample(X, y)

    print('X oversampled shape: %s' % str(X.shape))
    print('y oversampled shape: %s' % str(y.shape))

    #results = results.append(cross_validation(X, y, y_binary, X_test, y_test, y_test_binary, channel_set))
    results = results.append(cross_validation(X, y, y_binary, win_ids, channel_set))

    print(results)

    results.to_csv('NA12878/CV_results_' + balancing + '.csv', sep='\t')


def plot_results():

    source = pd.read_csv(filepath_or_buffer='NA12878/CV_results_' + balancing + '.csv', delimiter='\t')

    import numpy as np
    import matplotlib.pyplot as plt

    means = source.groupby('channels')['average_precision_score'].agg(np.mean).sort_values()
    print('average_precision_score mean: %d' % means)
    std = source.groupby('channels')['average_precision_score'].agg(np.std)
    print('average_precision_score mean: %d' % std)
    ind = np.arange(len(list(means.index)))  # the x locations for the groups
    width = 0.50  # the width of the bars: can also be len(x) sequence

    p1 = plt.bar(ind, means, width, yerr=std)

    plt.ylabel('average_precision_score')
    plt.title('average_precision_score per channel set')
    plt.xticks(ind, list(means.index))
    plt.xticks(rotation=45, horizontalalignment='right')
    # plt.yticks(np.arange(0.8, 1))
    # plt.legend((p1[0]), ('Bar'))
    plt.ylim(bottom=0.8)
    plt.tight_layout()

    plt.savefig('NA12878/Plots/Results_' + balancing + '.png', bbox_inches='tight')
    plt.close()


def main():
    # get_channel_labels()
    run_cv()
    plot_results()


if __name__ == '__main__':
    main()
