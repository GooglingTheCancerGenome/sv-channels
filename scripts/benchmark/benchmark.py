# Imports

import numpy as np
import sys
import os
import gzip
import pickle

from collections import Counter, defaultdict

# Keras imports
from keras.utils.np_utils import to_categorical
from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution1D, Flatten, MaxPooling1D
from keras.optimizers import Adam
from keras.models import load_model
from keras import backend as K

from mcfly import modelgen, find_architecture

# TensorFlow import
import tensorflow as tf

# Pandas import
import pandas as pd

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = False

# sample_name = 'G1'
# context_dict = {'G1': 'train', 'N1': 'real'}

sample_name = 'NA12878'
context_dict = {'NA12878': 'real'}

date = '060818'  # '010618'
mode = 'Test'

if HPC_MODE:
    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/' + date + '/' + mode + 'Data_' + date
else:
    base_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/' + date + '/' + mode + 'Data_' + date


# def check_units(y_true, y_pred):
#     if y_pred.shape[1] != 1:
#         y_pred = y_pred[:, 1:2]
#         y_true = y_true[:, 1:2]
#     return y_true, y_pred
#
#
# def precision(y_true, y_pred):
#     y_true, y_pred = check_units(y_true, y_pred)
#     true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
#     predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
#     precision = true_positives / (predicted_positives + K.epsilon())
#     return precision
#
#
# def recall(y_true, y_pred):
#     y_true, y_pred = check_units(y_true, y_pred)
#     true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
#     possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
#     recall = true_positives / (possible_positives + K.epsilon())
#     return recall
#
#
# def f1(y_true, y_pred):
#     def recall(y_true, y_pred):
#         true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
#         possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
#         recall = true_positives / (possible_positives + K.epsilon())
#         return recall
#
#     def precision(y_true, y_pred):
#         true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
#         predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
#         precision = true_positives / (predicted_positives + K.epsilon())
#         return precision
#
#     y_true, y_pred = check_units(y_true, y_pred)
#     precision = precision(y_true, y_pred)
#     recall = recall(y_true, y_pred)
#     return 2 * ((precision * recall) / (precision + recall + K.epsilon()))


def load_windows(sample_list, base_dir):
    # Load channel data and labels by chromosome

    # Leaving out chromosome Y and MT for the moment
    # chr_list = list(map(str, np.arange(1, 23)))
    # chr_list.append('X')

    # chr_list = list(map(str, np.arange(1, 18)))
    chr_list = ['22']
    # print(chr_list)

    data = []
    # labels = []

    for sample_name in sample_list:

        print('Considering sample %s' % sample_name)

        filename_suffix = '_channel_maker_' + context_dict[sample_name] + '_germline.npy.gz'

        for i in chr_list:
            print('Loading windows for Chr%s' % i)

            data_file = base_dir + '/' + sample_name + '/ChannelData/' + str(i) + filename_suffix

            with gzip.GzipFile(data_file, "rb") as f:
                data_mat = np.load(f)
                # Only loading 5 top channels
                data_mat = data_mat[:, [0, 1, 2, 19, 20], :]

            f.close()

            # if mode == 'Training':
            #     label_file = base_dir + '/' + sample_name + '/LabelData/' + sample_name + '_' + str(i) + '_label.npy.gz'
            # elif mode == 'Test':
            #     label_file = base_dir + '/' + sample_name + '/LabelData/' + str(i) + '_label_ci_full_overlap.npy.gz'
            #
            # with gzip.GzipFile(label_file, "rb") as f:
            #     data_lab = np.load(f)
            # f.close()

            data.extend(data_mat)
            # labels.extend(data_lab)

            # print('Length of data:%d, label length:%d' % (len(data_mat), len(data_lab)))
            # assert len(data_mat) == len(data_lab)

    # print(Counter(labels))
    # assert len(data) == len(labels)

    data = np.array(data)
    # data = data[:,[1,2,3,20,21],:]

    # labels = np.array(labels)

    print(data.shape)

    return data, chr_list


def load_labels(sampleName):

    print(f'Loading labels for {sampleName}')

    if not HPC_MODE:
        channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
    else:
        #Add HPC path to labels
        channel_dir = ''

    output_dir = '/'.join((channel_dir, sampleName, 'label_npy'))

    pickle_file = '/'.join((output_dir, 'labels.pickle.gz'))
    with gzip.GzipFile(pickle_file, "rb") as f:
        labels = pickle.load(f)
    f.close()

    return labels


def remove_label_from_data(data, labels, label_to_remove):

    labels = np.array(labels)
    no_idx = np.where(labels != label_to_remove)
    data = data[no_idx]
    labels = labels[no_idx]

    return data, labels


def relabel_deletions(labels):
    # Relabel deletions
    del_idx = np.where(np.logical_or(np.array(labels) == 'DEL_start',
                                     np.array(labels) == 'DEL_end'))
    # print(del_idx)
    # print(len(del_idx))
    # print(len(labels))
    labels[del_idx] = 'DEL'
    return labels


def balance_data(data, labels):
    cnt_lab = Counter(labels)
    min_v = min([v for k, v in cnt_lab.items()])

    # print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))

    # Balance data

    data_balanced = []
    labels_balanced = []

    for l in cnt_lab.keys():
        # print(l)
        iw = np.where(np.array(labels) == l)
        # print(iw[0])

        # Problem: random sampling of data
        # ii = np.random.choice(a=iw[0], size=min_v, replace=False)
        ii = iw[0][:min_v]

        data_balanced.extend(data[ii])
        labels_balanced.extend(labels[ii])

    # print(Counter(labels_balanced))

    X = np.array(data_balanced)
    y = np.array(labels_balanced)

    del data
    del labels
    del data_balanced
    del labels_balanced

    return X, y


def split_data(data, labels):
    # Split into training, validation and test set 60/20/20
    cnt_lab = Counter(labels)
    print(cnt_lab)
    n_lab = [v for v in cnt_lab.values()][0]
    # print(n_lab)

    i_train = int(n_lab * 0.6)
    i_val = i_train + int(n_lab * 0.2)
    # print(i_train)
    # print(i_val)

    X = defaultdict(list)
    y = defaultdict(list)

    X_np = defaultdict(np.array)
    y_np = defaultdict(np.array)

    for l in cnt_lab.keys():
        iw = np.where(labels == l)
        # print(iw[0])
        # print(iw[0][:i_train])
        X['train'].extend(data[iw[0][:i_train]])
        X['val'].extend(data[iw[0][i_train:i_val]])
        X['test'].extend(data[iw[0][i_val:]])

        y['train'].extend(labels[iw[0][:i_train]])
        y['val'].extend(labels[iw[0][i_train:i_val]])
        y['test'].extend(labels[iw[0][i_val:]])

    X_np['train'] = np.array(X['train'])
    X_np['val'] = np.array(X['val'])
    X_np['test'] = np.array(X['test'])

    y_np['train'] = np.array(y['train'])
    y_np['val'] = np.array(y['val'])
    y_np['test'] = np.array(y['test'])

    X_np = normalize(X_np)

    X_np = differentiate(X_np)

    X_np = transpose_dataset(X_np)

    # print(X_np['train'].shape)
    # print(X_np['val'].shape)
    # print(X_np['test'].shape)
    #
    print(y_np['train'].shape)
    print(y_np['val'].shape)
    print(y_np['test'].shape)

    return X_np, y_np


def transpose_dataset(X):
    for key in X.keys():
        image = []
        for i in range(0, len(X[key] - 1)):
            image.append(X[key][i].transpose())
        X[key] = np.array(image)

    return X


def normalize(X):
    number_channels = X['train'].shape[1]
    coverage = np.median(X['train'][:, 0, :])

    for key in X.keys():
        for i in range(0, len(X[key])):
            for j in range(0, number_channels):
                X[key][i][j] = np.divide(X[key][i][j], coverage) * 100

    return X


def differentiate(X):
    for key in X.keys():
        for i in range(0, X[key].shape[0]):
            for j in range(0, X[key].shape[1]):
                for k in range(0, X[key].shape[2] - 1):
                    X[key][i][j][k] = X[key][i][j][k + 1] - X[key][i][j][k]
        X[key] = X[key][:, :, :-1]

    return X


def run_model(X, y):
    # Prepare classes:

    classlabels = sorted(list(set(y['train'])))
    print(classlabels)
    mapclasses = {classlabels[i]: i for i in range(len(classlabels))}

    y['train'] = np.array([mapclasses[c] for c in y['train']], dtype='int')
    y['val'] = np.array([mapclasses[c] for c in y['val']], dtype='int')
    y['test'] = np.array([mapclasses[c] for c in y['test']], dtype='int')

    y_binary = dict()
    y_binary['train'] = to_categorical(y['train'])
    y_binary['val'] = to_categorical(y['val'])
    y_binary['test'] = to_categorical(y['test'])

    np.random.seed = 321
    num_classes = y_binary['train'].shape[1]

    models = modelgen.generate_models(X['train'].shape,
                                      num_classes,
                                      number_of_models=1,
                                      model_type='CNN',
                                      cnn_min_layers=1,
                                      cnn_max_layers=3,
                                      cnn_min_fc_nodes=4,
                                      cnn_max_fc_nodes=6,
                                      low_lr=4, high_lr=4,
                                      metrics=['accuracy'])  # , precision, recall, f1])

    # for model, params, model_types in models:
    #     print(params)
    #     model.summary()

    # Define directory where the results, e.g. json file, will be stored
    resultpath = os.path.join(base_dir, 'data/models_DEL')
    if not os.path.exists(resultpath):
        os.makedirs(resultpath)
    # print(resultpath)

    outputfile = os.path.join(resultpath, 'modelcomparison.json')
    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(X['train'], y_binary['train'],
                                                                                      X['val'], y_binary['val'],
                                                                                      models, nr_epochs=1,
                                                                                      subset_size=X['train'].shape[0],
                                                                                      verbose=False,
                                                                                      outputfile=outputfile,
                                                                                      metric='accuracy')
    # print('Details of the training process were stored in ', outputfile)

    best_model_index = np.argmax(val_accuracies)

    best_model, best_params, best_model_types = models[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    model_path = os.path.join(resultpath, 'best_model')
    # print(model_path)

    best_model.save(model_path)

    model_reloaded = load_model(model_path)

    probs = model_reloaded.predict_proba(X['test'], batch_size=1)
    get_confusion_matrix(probs, y_binary['test'], classlabels)

    return model_path, y_binary


def test_model(X, y_binary, model_path):
    model_reloaded = load_model(model_path)  # , custom_objects={'precision': precision, 'recall': recall, 'f1': f1})

    score_test = model_reloaded.evaluate(X['test'], y_binary['test'], verbose=False)
    print('Test loss and accuracy of best model: ' + str(score_test))


def get_confusion_matrix(probs, y_binary, labels):
    predicted = probs.argmax(axis=1)
    y_index = y_binary.argmax(axis=1)
    # Rows: true, columns: predicted
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in labels], fill_value=0)

    print('Confusion matrix:')
    print(confusion_matrix)

    for l in labels:
        # print(confusion_matrix.loc[l,:])
        # print(confusion_matrix.loc[:,l])

        # label_correct = confusion_matrix.loc[l, l]
        label_precision = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[l, :]) * 100)
        label_recall = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[:, l]) * 100)
        label_F1 = 2 * (label_precision * label_recall) / (label_precision + label_recall)

        print(f'{l} -> Precision:{label_precision}%, Recall:{label_recall}%, F1:{label_F1}')


def main():

    data, chr_list = load_windows(sample_list=context_dict.keys(),
                        base_dir=base_dir)
    labels_pickle = load_labels(sample_name)

    for key in labels_pickle:

        if key != 'id':

            print(f'Considering labels {key}')

            labels = []
            for chr in chr_list:
                labels.extend(labels_pickle[key][chr])

            data_filtered, labels_filtered = remove_label_from_data(data, labels, 'UK')

            # random sample?
            data_balanced, labels_balanced = balance_data(data_filtered, labels_filtered)

            X, y = split_data(data_balanced, labels_balanced)
            model_path, y_binary = run_model(X, y)
            # print(model_path)

            test_model(X, y_binary, model_path)


if __name__ == '__main__':
    main()
