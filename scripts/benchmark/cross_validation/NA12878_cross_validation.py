# Imports
import gzip
import os

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


HPC_MODE = False
sample_name = 'NA12878'
date = '231118'
label_type = 'Mills2011_nanosv'
datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari' if HPC_MODE else '/Users/lsantuari/Documents'
datapath_training =  datapath_prefix+'/Processed/Test/'+\
           date+'/TestData_'+date+'/'+sample_name+'/TrainingData/'
datapath_test =  datapath_prefix+'/Processed/Test/'+\
           date+'/TestData_'+date+'/'+sample_name+'/TestData/'

def data(datapath):

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

    # print(X.shape)
    # print(y.shape)
    # print(y.shape)
    # print(win_ids.shape)

    return X, y, y_binary, win_ids


def create_model(X, y_binary):

    models = modelgen.generate_models(X.shape,
                                      y_binary.shape[1],
                                      number_of_models = 1,
                                      model_type = 'CNN',
                                      cnn_min_layers=2,
                                      cnn_max_layers=2,
                                      cnn_min_filters = 4,
                                      cnn_max_filters = 4,
                                      cnn_min_fc_nodes=6,
                                      cnn_max_fc_nodes=6,
                                      low_lr=2, high_lr=2,
                                      low_reg=4, high_reg=4,
                                      kernel_size = 9)

    # i = 0
    # for model, params, model_types in models:
    #     print('model ' + str(i))
    #     i = i + 1
    #     print(params)
    #     model.summary()

    return models


def cross_validation(X, y, y_binary, X_test, y_test, y_test_binary):

    results = defaultdict(list)

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

        # split into train/test sets
        xtrain, xval, ytrain_binary, yval = train_test_split(xtrain, ytrain_binary, test_size=0.2, random_state=2)

        # Clear model, and create it
        model = None
        model = create_model(X, y_binary)

        # Debug message I guess
        # print "Training new iteration on " + str(xtrain.shape[0]) + " training samples, " + str(xval.shape[0]) + " validation samples, this may be a while..."

        history, model = train_model(model, xtrain, ytrain_binary, xval, yval)

        accuracy_history = history.history['acc']
        val_accuracy_history = history.history['val_acc']
        print("Last training accuracy: " + str(accuracy_history[-1]) + ", last validation accuracy: " + str(
            val_accuracy_history[-1]))

        score_test = model.evaluate(xtest, ytest_binary, verbose=False)
        print('Test loss and accuracy of best model: ' + str(score_test))

        evaluate_model(model, xtest, ytest_binary, results)

    print(results)

def train_model(model, xtrain, ytrain, xval, yval):

    train_set_size = xtrain.shape[0]

    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(xtrain, ytrain,
                                                                                      xval, yval,
                                                                                      model, nr_epochs=1,
                                                                                      subset_size=train_set_size,
                                                                                      verbose=False)

    best_model_index = np.argmax(val_accuracies)
    best_model, best_params, best_model_types = model[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    # We make a copy of the model, to start training from fresh
    nr_epochs = 1
    datasize = train_set_size  # Change in `X_train.shape[0]` if training complete data set
    history = best_model.fit(xtrain, ytrain,
                             epochs=nr_epochs, validation_data=(xval, yval),
                             verbose=False)

    return history, best_model


def evaluate_model(model, X_test, ytest_binary, results):

    mapclasses = {'DEL_start': 1, 'DEL_end': 0, 'noSV': 2}
    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    # print(dict_sorted)
    class_labels = [i[0] for i in dict_sorted]

    n_classes = ytest_binary.shape[1]
    # print(y_binarized)
    print(n_classes)

    probs = model.predict_proba(X_test, batch_size=1, verbose=False)

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
    print('Average precision score, micro-averaged over all classes: {0:0.2f}'
          .format(average_precision["micro"]))

    results['average_precision'].append(average_precision["micro"])

    # plt.figure()
    # plt.step(recall['micro'], precision['micro'], color='b', alpha=0.2,
    #          where='post')
    # plt.fill_between(recall["micro"], precision["micro"], alpha=0.2, color='b')
    #
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.ylim([0.0, 1.05])
    # plt.xlim([0.0, 1.0])
    # plt.title(
    #     'Average precision score, micro-averaged over all classes: AP={0:0.2f}'
    #         .format(average_precision["micro"]))
    #
    # plt.show()

    # from itertools import cycle
    # # setup plot details
    # colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])
    #
    # plt.figure(figsize=(7, 8))
    # f_scores = np.linspace(0.2, 0.8, num=4)
    # lines = []
    # labels = []
    # for f_score in f_scores:
    #     x = np.linspace(0.01, 1)
    #     y = f_score * x / (2 * x - f_score)
    #     l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
    #     plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))
    #
    # lines.append(l)
    # labels.append('iso-f1 curves')
    # l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
    # lines.append(l)
    # labels.append('micro-average Precision-recall (area = {0:0.2f})'
    #               ''.format(average_precision["micro"]))
    #
    # for i, color in zip(range(n_classes), colors):
    #     l, = plt.plot(recall[i], precision[i], color=color, lw=2)
    #     lines.append(l)
    #     labels.append('Precision-recall for class {0} (area = {1:0.2f})'
    #                   ''.format(class_labels[i], average_precision[i]))
    #
    # fig = plt.gcf()
    # fig.subplots_adjust(bottom=0.25)
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.title('Extension of Precision-Recall curve to multi-class')
    # plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))
    #
    # plt.show()


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
    #     y_test_class = np.array([1 if i[iter_class] == 1 else 0 for i in y_test])
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
    #     plt.show()

    return results


def main():

    # Load the data
    X, y, y_binary, win_ids = data(datapath_training)

    X_test, y_test, y_test_binary, win_ids_test = data(datapath_test)

    cross_validation(X, y, y_binary, X_test, y_test, y_test_binary)

if __name__ == '__main__':

    main()