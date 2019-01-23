# Imports
import gzip
import os

import numpy as np

import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

import math

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

#import altair as alt

# Bokeh import
#from bokeh.io import show, output_file
#from bokeh.plotting import figure

import logging


HPC_MODE = False
sample_name = 'OC'
date = '070119'
label_type = 'bpi'
datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari' if HPC_MODE else '/Users/lsantuari/Documents'
datapath_training =  datapath_prefix+'/Processed/Test/'+\
           date+'/TestData_'+date+'/'+sample_name+'/TrainingData/'
datapath_test =  datapath_prefix+'/Processed/Test/'+\
           date+'/TestData_'+date+'/'+sample_name+'/TestData/'


def get_classes(labels):
    return sorted(list(set(labels)))


def get_channel_labels():
    # Fill labels for legend

    labels = list()
    labels.append("coverage")
    labels.append("#left_clipped_reads")
    labels.append("#right_clipped_reads")
    labels.append("INV_before")
    labels.append("INV_after")
    labels.append("DUP_before")
    labels.append("DUP_after")
    labels.append("TRA_opposite")
    labels.append("TRA_same")

    for direction in ['Forward', 'Reverse']:
        for clipped in ['Left', 'Right', 'Not']:
            for value in ['sum', 'num', 'median', 'outliers']:
                labels.append(direction + '_' + clipped + '_Clipped_' + value)

    labels.append("#left split reads")
    labels.append("#right split reads")

    for clipped in ['L', 'R']:
        for value in ['sum', 'num', 'median']:
            labels.append(clipped + '_SplitRead_' + value)

    labels.append("GC")
    labels.append("Mappability")
    labels.append("One_hot_Ncoding")

    for k, l in enumerate(labels):
         logging.info(str(k) + ':' + l)

    return labels


def get_channel_labels_TN():
    # Fill labels for legend

    labels = list()

    for type in ['Tumor:', 'Normal:']:

        labels.append(type+"coverage")
        labels.append(type+"#left_clipped_reads")
        labels.append(type+"#right_clipped_reads")
        labels.append(type+"INV_before")
        labels.append(type+"INV_after")
        labels.append(type+"DUP_before")
        labels.append(type+"DUP_after")
        labels.append(type+"TRA_opposite")
        labels.append(type+"TRA_same")

        for direction in ['Forward', 'Reverse']:
            for clipped in ['Left', 'Right', 'Not']:
                for value in ['outliers']:
                    labels.append(type + direction + '_' + clipped + '_Clipped_' + value)

        labels.append(type + "#left split reads")
        labels.append(type + "#right split reads")

        for clipped in ['L', 'R']:
            for value in ['sum', 'num', 'median']:
                labels.append(type + clipped + '_SplitRead_' + value)

    # labels.append("GC")
    labels.append("Mappability")

    for nuc in ['A', 'T', 'C', 'G', 'N']:
        labels.append("One_hot_"+nuc+"_encoding")

    for k, l in enumerate(labels):
         logging.info(str(k) + ':' + l)

    return labels


def set_figure_size(plt):

    # plt.tight_layout()

    F = plt.gcf()
    # Now check everything with the defaults:
    DPI = F.get_dpi()
    logging.info(
        "DPI:", DPI)
    DefaultSize = F.get_size_inches()
    logging.info(
        "Default size in Inches", DefaultSize)
    logging.info(
        "Which should result in a %i x %i Image" % (DPI * DefaultSize[0], DPI * DefaultSize[1]))

    F.set_figwidth(DefaultSize[0] * 5)
    F.set_figheight(DefaultSize[1] * 4)
    Size = F.get_size_inches()
    logging.info(
        "Size in Inches", Size)


def plot_channels(X, y, ids):

    output_dir = 'OC/channel_plots'
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    labels = get_channel_labels_TN()
    number_channels = len(labels)

    classes = get_classes(y)

    for c in classes:

        output_dir_class = os.path.join(output_dir, c)
        if not os.path.isdir(output_dir_class):
            os.mkdir(output_dir_class)

        c_idx = np.where(y == c)[0]
        c_idx_random = np.random.choice(c_idx, size=10, replace=False)

        for i in c_idx_random:

            fig_name = y[i] + '_' + 'Chr' + ids[i].replace(':', '-')

            logging.info(y[i], 'id:', 'Chr' + ids[i])
            plt.title('Class: ' + y[i] + ' ' + 'Position: Chr' + ids[i], fontsize=30)
            plt.ylim([0, number_channels])
            plt.yticks(np.arange(number_channels), labels, fontsize=15)
            plt.xticks(fontsize=30)

            plt.vlines(x=0, ymin=0, ymax=number_channels, color='black')

            for j in range(number_channels - 1, -1, -1):

                if sum(X[i,:,j]) != 0:
                    X_win = ((X[i,:,j] - min(X[i,:,j])) / max(X[i,:,j]))
                else:
                    X_win = X[i,:,j]

                Z = [x + j for x in X_win]

                plt.plot(np.arange(-100, 100), Z, label=y[j], linewidth=2)
                plt.fill_between(np.arange(-100, 100), Z, j, alpha=.5, interpolate=True)

            set_figure_size(plt)
            plt.savefig(os.path.join(output_dir_class, fig_name + '.png'))
            plt.clf()
            plt.close()


def data(datapath):

    dataset_type = '_balanced'
    data_input_file = datapath + sample_name + '_' + label_type + dataset_type + '.npz.gz'

    with gzip.GzipFile(data_input_file, "rb") as f:
        npzfiles = np.load(f)
        X = npzfiles['X']
        y = npzfiles['y']
        y_binary = npzfiles['y_binary']
        win_ids = npzfiles['ids']

    # logging.info(X.shape)
    # logging.info(y.shape)
    # logging.info(y.shape)
    # logging.info(win_ids.shape)

    # idx = np.arange(0,9)
    # idx = np.append(idx, np.arange(33,35))
    # idx = np.append(idx, np.arange(41, 44))
    # idx = np.append(idx,[12,16,20,24,28,32])

    return X, y, y_binary, win_ids


def remove_classes(X, y, win_ids, removed_labels):

    for l in removed_labels:

        keep = np.where(np.array(y) != l)
        X = X[keep]
        y = y[keep]
        win_ids = win_ids[keep]

    classes = get_classes(y)

    mapclasses = dict()
    for i, c in enumerate(classes):
        mapclasses[c] = i

    y_num = np.array([mapclasses[c] for c in y], dtype='int')
    y_binary = to_categorical(y_num)

    # print(y_binary)
    # print(X.shape)
    # print(y.shape)
    # print(win_ids.shape)

    return X, y, y_binary, win_ids


def create_model(X, y_binary):

    # models = modelgen.generate_models(X.shape,
    #                                   y_binary.shape[1],
    #                                   number_of_models = 1,
    #                                   model_type = 'CNN',
    #                                   cnn_min_layers=2,
    #                                   cnn_max_layers=2,
    #                                   cnn_min_filters = 4,
    #                                   cnn_max_filters = 4,
    #                                   cnn_min_fc_nodes=6,
    #                                   cnn_max_fc_nodes=6,
    #                                   low_lr=2, high_lr=2,
    #                                   low_reg=1, high_reg=1,
    #                                   kernel_size = 7)

    models = modelgen.generate_models(X.shape,
                                      y_binary.shape[1],
                                      number_of_models = 1,
                                      model_type = 'CNN',
                                      cnn_min_layers=4,
                                      cnn_max_layers=4,
                                      cnn_min_filters = 6,
                                      cnn_max_filters = 6,
                                      cnn_min_fc_nodes=12,
                                      cnn_max_fc_nodes=12,
                                      low_lr=2, high_lr=2,
                                      low_reg=1, high_reg=1,
                                      kernel_size = 7)

    # i = 0
    # for model, params, model_types in models:
    #     logging.info('model ' + str(i))
    #     i = i + 1
    #     logging.info(params)
    #     model.summary()

    return models


def cross_validation(X, y, y_binary, channels):

    results = pd.DataFrame()

    # From https://medium.com/@literallywords/stratified-k-fold-with-keras-e57c487b1416
    kfold_splits = 10

    # Instantiate the cross validator
    skf = StratifiedKFold(n_splits=kfold_splits, shuffle=True)

    # Loop through the indices the split() method returns
    for index, (train_indices, test_indices) in enumerate(skf.split(X, y)):

        logging.info("Training on fold " + str(index + 1) + "/10...")

        # Generate batches from indices
        xtrain, xtest = X[train_indices], X[test_indices]
        ytrain, ytest = y[train_indices], y[test_indices]
        ytrain_binary, ytest_binary = y_binary[train_indices], y_binary[test_indices]

        # split into train/validation sets
        xtrain, xval, ytrain_binary, yval = train_test_split(xtrain, ytrain_binary,
                                                             test_size=0.2, random_state=2)

        # Clear model, and create it
        model = None
        model = create_model(X, y_binary)

        # Debug message I guess
        logging.info ("Training new iteration on " + str(xtrain.shape[0]) + " training samples, " +
         str(xval.shape[0]) + " validation samples, this may take a while...")

        history, model = train_model(model, xtrain, ytrain_binary, xval, yval)

        accuracy_history = history.history['acc']
        val_accuracy_history = history.history['val_acc']
        logging.info("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
            val_accuracy_history[-1]))

        score_test = model.evaluate(xtest, ytest_binary, verbose=False)
        logging.info('Test loss and accuracy of best model: ' + str(score_test))

        results = evaluate_model(model, xtest, ytest, ytest_binary, results, index, channels,
                                 train_set_size=xtrain.shape[0],
                                 validation_set_size = xval.shape[0]
        )
        #evaluate_model(model, X_test, y_test_binary, results, index, channels)

    return results


def train_model(model, xtrain, ytrain, xval, yval):

    train_set_size = xtrain.shape[0]

    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(xtrain, ytrain,
                                                                                      xval, yval,
                                                                                      model, nr_epochs=1,
                                                                                      subset_size=train_set_size,
                                                                                      verbose=False)

    best_model_index = np.argmax(val_accuracies)
    best_model, best_params, best_model_types = model[best_model_index]
    # logging.info(best_model_index, best_model_types, best_params)

    nr_epochs = 10
    history = best_model.fit(xtrain, ytrain,
                             epochs=nr_epochs, validation_data=(xval, yval),
                             verbose=False)

    return history, best_model


def evaluate_model(model, X_test, y_test, ytest_binary, results, cv_iter, channels,
                   train_set_size, validation_set_size):

    #Generate classes
    classes = sorted(list(set(y_test)))
    mapclasses = dict()
    for i, c in enumerate(classes):
        mapclasses[c] = i

    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    # logging.info(dict_sorted)
    class_labels = [i[0] for i in dict_sorted]

    n_classes = ytest_binary.shape[1]
    # logging.info(ytest_binary)
    # logging.info(n_classes)

    probs = model.predict_proba(X_test, batch_size=1, verbose=False)

    # generate confusion matrix
    labels = sorted(list(set(y_test)))
    predicted = probs.argmax(axis=1)
    y_index = ytest_binary.argmax(axis=1)
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in labels], fill_value=0)
    logging.info(confusion_matrix)
    confusion_matrix.to_csv(path_or_buf='OC/OC_confusion_matrix_cv_iter_' + str(cv_iter + 1) + '.csv')

    print(np.diag(confusion_matrix))
    print(confusion_matrix.sum(axis=1))
    print(confusion_matrix)
    # logging.info('Precision: %d' % int(np.diag(confusion_matrix) / confusion_matrix.sum(axis=1) * 100))
    # logging.info('Recall: %d' % int(np.diag(confusion_matrix)/confusion_matrix.sum(axis=0)*100))

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
    logging.info('Average precision score, micro-averaged over all classes: {0:0.2f}'
          .format(average_precision["micro"]))

    results = results.append({
        "channels": channels,
        "fold": cv_iter+1,
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

    plt.savefig('OC/Precision_Recall_avg_prec_score_Iter_'+str(cv_iter)+'_'+channels+'.png', bbox_inches='tight')
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

    plt.savefig('OC/Precision_Recall_avg_prec_score_per_class_Iter_' +
                str(cv_iter) +'_'+channels+'.png', bbox_inches='tight')
    plt.close()

    # for iter_class in mapclasses.values():
    #
    #     predicted = probs.argmax(axis=1)
    #     #logging.info(predicted)
    #     y_pred_class = np.array([1 if i == iter_class else 0 for i in predicted])
    #     #logging.info(y_pred_class)
    #
    #     # keep probabilities for the positive outcome only
    #     probs_class = probs[:, iter_class]
    #     #logging.info(probs_class)
    #
    #     #logging.info(y_test)
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
    #     logging.info('f1=%.3f auc=%.3f average_precision_score=%.3f' % (f1, auc_value , ap))
    #     # plot no skill
    #     plt.plot([0, 1], [0.5, 0.5], linestyle='--')
    #     # plot the roc curve for the model
    #     plt.plot(recall, precision, marker='.')
    #     # show the plot
    #     plt.savefig('Plots/Precision_Recall_multiclass_Iter_'+str(cv_iter)+'_'+channels+'.png', bbox_inches='tight')

    return results


def run_cv():

    labels = get_channel_labels_TN()

    results = pd.DataFrame()

    channels = 'all'
    logging.info('Running cv with '+channels+' channels:')
    for i, l in enumerate(labels):
        logging.info(str(i) + ':' + l)

    # Load the data
    X, y, y_binary, win_ids = data(datapath_training)
    # BND
    # removed_labels = ['INV_start', 'INV_end', 'DEL_start', 'DEL_end', 'DUP_start', 'DUP_end']
    # X, y, y_binary, win_ids = remove_classes(X, y, win_ids, removed_labels)

    # X_test, y_test, y_test_binary, win_ids_test = data(datapath_test)

    results = results.append(cross_validation(X, y, y_binary, channels))

    logging.info(results)
    results.to_csv("OC/CV_results.csv", sep='\t')


def plot_results():

    source = pd.read_csv(filepath_or_buffer='OC/CV_results.csv', delimiter='\t')

    import numpy as np
    import matplotlib.pyplot as plt

    means = source.groupby('channels')['average_precision_score'].agg(np.mean).sort_values()
    logging.info(means)
    std = source.groupby('channels')['average_precision_score'].agg(np.std)
    logging.info(std)
    ind = np.arange(len(list(means.index)))  # the x locations for the groups
    width = 0.50  # the width of the bars: can also be len(x) sequence

    p1 = plt.bar(ind, means, width, yerr=std)

    plt.ylabel('average_precision_score')
    plt.title('average_precision_score per channel set')
    plt.xticks(ind, list(means.index))
    plt.xticks(rotation=45,  horizontalalignment='right')
    #plt.yticks(np.arange(0.8, 1))
    #plt.legend((p1[0]), ('Bar'))
    plt.ylim(bottom=0.8)
    plt.tight_layout()

    plt.savefig('OC/Results.png', bbox_inches='tight')
    plt.close()


def main():


    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename='OC/logfile.log',
        level=logging.INFO)

    # get_channel_labels_TN()
    run_cv()
    # plot_results()

    # logging.info('Loading data...')
    # X, y, y_binary, win_ids = data(datapath_training)
    # logging.info(X.shape)
    # logging.info('Plotting channels...')
    # plot_channels(X, y, win_ids)


if __name__ == '__main__':

    main()