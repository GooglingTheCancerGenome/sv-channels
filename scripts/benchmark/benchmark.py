# Imports

import numpy as np
import sys
import os, errno
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

# Flags

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = True

# Artificial mode?
ART_MODE = False

# Model creation mode?
CREATE_MODELS = True

# Already differentiated?
DIFFERENTIATED = False

channel_selection = 'AllChannel'

# sample_name = 'G1'
# context_dict = {'G1': 'train', 'N1': 'real'}

# sample_name = 'PATIENT1'
context_dict = {'NA12878': 'real', 'PATIENT1': 'real', 'PATIENT2': 'real', 'G1': 'train', 'N1': 'real'}

date = '130918'  # '060818'
mode = 'Test'

if HPC_MODE:

    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/' + mode + 'Data_' + date

    base_dir_art = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TrainingData_' + date

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


def get_chr_list():
    # Leaving out chromosome Y and MT for the moment
    chr_list = list(map(str, np.arange(1, 23)))
    chr_list.append('X')
    #chr_list = list(map(str, np.arange(1, 11)))
    # DEBUG:
    # chr_list = ['12']
    return chr_list


def load_windows(sample_list, base_dir):
    # Load channel data and labels by chromosome

    chr_list = get_chr_list()

    # chr_list = list(map(str, np.arange(1, 18)))
    # chr_list = ['22']
    # print(chr_list)

    data = []
    # labels = []

    for sample_name in sample_list:

        print('Loading windows for sample %s' % sample_name)

        filename_suffix = '_channel_maker_' + context_dict[sample_name] + '_germline.npy.gz'

        for i in chr_list:
            print('Loading windows for Chr%s' % i)

            data_file = base_dir + '/' + sample_name + '/ChannelData/' + str(i) + filename_suffix

            with gzip.GzipFile(data_file, "rb") as f:
                data_mat = np.load(f)
                # Only loading 5 top channels
                # 90:110
                data_mat = data_mat[:, [0, 1, 2, 19, 20], ]
                # data_mat = data_mat[:, [0, 1, 2, 19, 20, 15, 24, 34], :]
                # 15 median PE track
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


def save_windows(data, sample_name):
    print('Saving windows for sample %s' % sample_name)

    # if not DIFFERENTIATED:
    print('Start processing windows...')
        # print('Start normalize', data.shape)
        # data = normalize_single(data)
        # print('Start differentiate', data.shape)
        # data = differentiate_single(data)

    print('Start transpose', data.shape)
    data = transpose_single(data)
    print('Finished', data.shape)

    dirname = base_dir + '/' + sample_name + '/ZipData'
    create_dir(dirname)

    #if not DIFFERENTIATED:
    data_file = dirname + '/' + sample_name + '_' + channel_selection + '_windows_diff'
    # else:
    #      data_file = dirname + '/' + sample_name + '_' + channel_selection + '_windows'

    np.save(data_file, data)
    os.system('gzip ' + data_file + '.npy')

    return data


def load_windows_from_zip(sample_name):

    #if DIFFERENTIATED:
    data_file = base_dir + '/' + sample_name + '/ZipData/' + sample_name + '_' + channel_selection + \
                    '_windows_diff.npy.gz'
    # else:
    #    data_file = base_dir + '/' + sample_name + '/ZipData/' + sample_name + '_' + channel_selection + \
    #                '_windows.npy.gz'

    with gzip.GzipFile(data_file, "rb") as f:
        data = np.load(f)
    f.close()

    return data


def load_labels(sampleName):
    print(f'Loading labels for {sampleName}')

    if not HPC_MODE:
        channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
    else:
        # Add HPC path to labels
        channel_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/' + mode + 'Data_' + date

    output_dir = '/'.join((channel_dir, sampleName, 'MultiLabelData'))

    pickle_file = '/'.join((output_dir, 'labels.pickle.gz'))
    with gzip.GzipFile(pickle_file, "rb") as f:
        labels = pickle.load(f)
    f.close()

    return labels


def load_artificial_labels(sample_list, base_dir):

    chr_list = get_chr_list()
    labels = []

    for sampleName in sample_list:
        print(f'Loading labels for {sampleName}')
        for i in chr_list:

            print('Loading windows for Chr%s' % i)

            label_file = base_dir + '/' + sampleName + '/LabelData/' + sampleName + '_' + str(i) + '_label.npy.gz'

            with gzip.GzipFile(label_file, "rb") as f:
                data_lab = np.load(f)
            f.close()

            labels.extend(data_lab)

    labels = np.array(labels)

    print(labels.shape)

    return labels


def remove_label_from_data(data, labels, label_to_remove):
    print('remove_label_from_data: Data shape:')
    print(data.shape)

    labels = np.array(labels)
    no_idx = np.where(labels != label_to_remove)

    data = data[no_idx]
    labels = labels[no_idx]

    return data, labels


def relabel_deletions(labels):
    # Relabel deletions
    labels = np.array(labels)
    # logical_or can be used only with two arguments
    del_idx = np.where(np.logical_or(np.array(labels) == 'DEL_start',
                                     np.array(labels) == 'DEL_end'))

    # del_idx = np.where(np.any((np.array(labels) == 'DEL_start',
    #                            np.array(labels) == 'DEL_end',
    #                            np.array(labels) == 'UK'), axis=0))

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

    #Take the minimum number of labels
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

    # if not DIFFERENTIATED:
    #     X_np = normalize(X_np)
    #     X_np = differentiate(X_np)
    #     X_np = transpose_dataset(X_np)

    # print(X_np['train'].shape)
    # print(X_np['val'].shape)
    # print(X_np['test'].shape)
    #
    print(y_np['train'].shape)
    print(y_np['val'].shape)
    print(y_np['test'].shape)

    return X_np, y_np


def transpose_dataset(X):
    print('Transpose')

    # def func(x0):
    #     return x0.transpose()

    for key in X.keys():
        # print(X[key].shape)
        # X[key] = np.apply_along_axis(func, 1, X[key])
        image = []
        for i in range(0, len(X[key])):
            image.append(X[key][i].transpose())
        X[key] = np.array(image)
        # print(X[key].shape)

    return X


def transpose_single(X):
    print('Transpose')

    image = []
    for i in range(0, len(X)):
        image.append(X[i].transpose())
    X = np.array(image)

    # def func(x0):
    #     return x0.transpose()
    #
    # X = np.apply_along_axis(func, 1, X)

    return X


def normalize(X):
    # print('Normalize:')
    # number_channels = X['train'].shape[1]
    coverage = np.median(X['train'][:, 0, :])

    def func(x0):
        return np.divide(x0, coverage) * 100

    for key in X.keys():
        # print(X[key].shape)
        X[key] = np.apply_along_axis(func, 1, X[key])
        # for i in range(0, len(X[key])):
        #     for j in range(0, number_channels):
        #         X[key][i][j] = np.divide(X[key][i][j], coverage) * 100
        # print(X[key].shape)
    return X


def normalize_single(X):
    print('Normalize:')

    coverage = np.median(X[:, 0, :])

    def func(x0):
        return np.divide(x0, coverage) * 100

    # number_channels = X.shape[1]
    #
    # for i in range(0, len(X)):
    #     for j in range(0, number_channels):
    #         X[i][j] = np.divide(X[i][j], coverage) * 100

    X = np.apply_along_axis(func, 1, X)
    return X


def differentiate(X):
    print('Differentiate:')

    def func(x0):
        return np.apply_along_axis(lambda x: x[1:] - x[:-1], 0, x0)

    for key in X.keys():
        # print(X[key].shape)
        print('Shape before differentiate: %s' % str(X[key].shape))
        X[key] = np.apply_along_axis(func, 2, X[key])

        # for i in range(0, X[key].shape[0]):
        #     for j in range(0, X[key].shape[1]):
        #         for k in range(0, X[key].shape[2] - 1):
        #             X[key][i][j][k] = X[key][i][j][k + 1] - X[key][i][j][k]
        # X[key] = X[key][:, :, :-1]
        print('Shape after differentiate: %s' % str(X[key].shape))
        # print(X[key].shape)

    return X


def differentiate_single(X_fast):
    print('Differentiate:')

    def func(x0):
        return np.apply_along_axis(lambda x: x[1:] - x[:-1], 0, x0)

    # X_slow = np.copy(X[:100])
    # X_fast = np.copy(X[:100])

    # print('X fast shape before differentiate: %s' % str(X_fast.shape))
    X_fast = np.apply_along_axis(func, 2, X_fast)
    #X_fast = X_fast[:, :, :-1]
    # print('X fast shape after differentiate: %s' % str(X_fast.shape))

    # for i in range(0, X.shape[0]):
    #     X[i] = np.apply_along_axis(lambda x: x[1:] - x[:-1], 0, X[i])

    # print('X slow shape before differentiate: %s' % str(X_slow.shape))
    # # for i in range(0, X.shape[0]):
    # for i in range(0, 100):
    #     for j in range(0, X_slow.shape[1]):
    #         # X_slow[i][j] = X_slow[i][j][1:] - X_slow[i][j][:(X_slow.shape[2] - 1)]
    #         for k in range(0, X_slow.shape[2] - 1):
    #             X_slow[i][j][k] = X_slow[i][j][k + 1] - X_slow[i][j][k]
    # X_slow = X_slow[:, :, :-1]
    # print('X slow shape after differentiate: %s' % str(X_slow.shape))

    # assert np.array_equal(X_slow, X_fast)

    return X_fast


def get_categorical_labels(y):

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

    return y_binary, classlabels


def run_model(X, y_binary, label_name, sample_name, iteration, classlabels, models):
    # Prepare classes:

    # for model, params, model_types in models:
    #     print(params)
    #     model.summary()

    # Define directory where the results, e.g. json file, will be stored
    resultpath = os.path.join(base_dir, 'ModelData' + channel_selection + '/' + sample_name + '/models_DEL')
    if not os.path.exists(resultpath):
        os.makedirs(resultpath)
    # print(resultpath)

    outputfile = os.path.join(resultpath, 'modelcomparison.json')
    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(X['train'], y_binary['train'],
                                                                                      X['val'], y_binary['val'],
                                                                                      models, nr_epochs=1,
                                                                                      subset_size=X['train'].shape[0],
                                                                                      verbose=True,
                                                                                      outputfile=outputfile,
                                                                                      metric='accuracy')
    # print('Details of the training process were stored in ', outputfile)

    best_model_index = np.argmax(val_accuracies)

    best_model, best_params, best_model_types = models[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    # the model path dependent on label_name
    model_path = os.path.join(resultpath, 'best_model_' + label_name + '_' + str(iteration))
    # print(model_path)

    best_model.save(model_path)

    # model_reloaded = load_model(model_path)

    probs = best_model.predict_proba(X['test'], batch_size=1)

    df = get_confusion_matrix(iteration, sample_name, label_name, probs, y_binary['test'], classlabels)

    # df = get_precision_recall_curve(sample_name, label_name, probs, y_binary['test'], classlabels)

    return model_path, df


def test_model(X, y_binary, model_path):
    model_reloaded = load_model(model_path)  # , custom_objects={'precision': precision, 'recall': recall, 'f1': f1})

    score_test = model_reloaded.evaluate(X['test'], y_binary['test'], verbose=True)
    print('Test loss and accuracy of best model: ' + str(score_test))


def get_confusion_matrix(iteration, sample_name, label_name, probs, y_binary, labels):

    predicted = probs.argmax(axis=1)
    y_index = y_binary.argmax(axis=1)

    # Rows: true, columns: predicted
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in labels], fill_value=0)

    print('Confusion matrix:')
    print(confusion_matrix)

    df_conf = pd.DataFrame()

    labels_list = sorted(list(set(confusion_matrix.index) & set(confusion_matrix.columns)))

    for l in labels_list:
        # print(confusion_matrix.loc[l,:])
        # print(confusion_matrix.loc[:,l])

        # label_correct = confusion_matrix.loc[l, l]
        label_precision = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[:, l]) * 100)
        label_recall = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[l, :]) * 100)
        label_F1 = 2 * (label_precision * label_recall) / (label_precision + label_recall)

        print(f'Iter:{iteration} {l} -> Precision:{label_precision}%, Recall:{label_recall}%, F1:{label_F1}')

        df_intres = pd.DataFrame({'iteration': [iteration], 'sample': [sample_name], 'evidence': [label_name],
                                  'label': [l], 'precision': [label_precision],
                                  'recall': [label_recall], 'F1': [label_F1]})
        df_conf = df_conf.append(df_intres)

    return df_conf


def get_precision_recall_curve(sample_name, label_name, probs, y_binary, labels):
    # setup uniformly separated cutoffs in interval [0.5, 1)
    for i in np.linspace(1.0 / len(labels), 1, num=50, endpoint=False):

        predicted = np.argwhere(probs > i)[:, 1]
        y_index = np.argwhere(y_binary > i)[:, 1]

        # Rows: true, columns: predicted
        confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
        confusion_matrix.index = [labels[i] for i in confusion_matrix.index]
        confusion_matrix.columns = [labels[i] for i in confusion_matrix.columns]
        confusion_matrix.reindex(columns=[l for l in labels], fill_value=0)

        print('Confusion matrix:')
        print(confusion_matrix)

        df_conf = pd.DataFrame()

        for l in labels:
            if l in confusion_matrix.index:
                # print(confusion_matrix.loc[l,:])
                # print(confusion_matrix.loc[:,l])

                # label_correct = confusion_matrix.loc[l, l]
                label_precision = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[:, l]) * 100)
                label_recall = np.around(confusion_matrix.loc[l, l] / sum(confusion_matrix.loc[l, :]) * 100)
                label_F1 = 2 * (label_precision * label_recall) / (label_precision + label_recall)

                print(f'Iter:{i} {l} -> Precision:{label_precision}%, Recall:{label_recall}%, F1:{label_F1}')

                df_intres = pd.DataFrame(
                    {'sample': [sample_name], 'evidence': [label_name], 'iteration': [i], 'label': [l],
                     'precision': [label_precision], 'recall': [label_recall], 'F1': [label_F1]})
                df_conf = df_conf.append(df_intres)

    return df_conf


def predict_windows_with_model():
    chr_list = get_chr_list()

    model_sample_name = 'NA12878'
    # label_name = 'Mills2011_PacBio_Moleculo_nanosv'
    # iteration = 1

    # label_name = 'Mills2011_nanosv_manta'
    # iteration = 3

    classlabels = ['DEL', 'noSV']
    resultpath = os.path.join(base_dir, 'ModelData' + channel_selection + '/' + model_sample_name + '/models_DEL')

    labels_pickle = load_labels(model_sample_name)

    for sample_name in ['PATIENT1', 'PATIENT2']:

        print(f'Considering sample {sample_name}')

        X = load_windows_from_zip(sample_name)

        # if not DIFFERENTIATED:
        #     X = normalize_single(X)
        #     X = differentiate_single(X)
        #     X = transpose_single(X)

        label_id = labels_pickle['id']

        label_id_list = []
        for k in chr_list:
            label_id_list.extend(label_id[k])
        label_id_list = np.array(label_id_list)
        # print(label_id_list)

        for label_name in labels_pickle:

            if label_name != 'id':

                print(f'Considering label {label_name}')

                for iteration in range(1, 11, 1):

                    print(f'Considering iteration {iteration}')

                    model_path = os.path.join(resultpath, 'best_model_' + label_name + '_' + str(iteration))
                    model_reloaded = load_model(model_path)

                    probs = model_reloaded.predict_proba(X, batch_size=1000, verbose=True)

                    predicted = probs.argmax(axis=1)
                    # print(predicted)
                    resulting_labels = [classlabels[i] for i in predicted]

                    # print(Counter(resulting_labels))
                    # print(label_id_list[resulting_labels == "DEL"])

                    # write BED
                    resulting_labels = np.array(resulting_labels)
                    # print(label_id_list[np.where(resulting_labels == "DEL")])
                    lines = []
                    for x in label_id_list[np.where(resulting_labels == "DEL")]:
                        lines.append(bytes(x['chromosome'] + '\t' + str(x['position']) + '\t' \
                                           + str(x['position'] + 1) + '\t' + 'DEL' + '\n', 'utf-8'))

                    outdir = os.path.join(base_dir, 'BED_results', sample_name)
                    create_dir(outdir)
                    outfile = os.path.join(outdir, sample_name + '_' + \
                                           label_name + '_' + str(iteration) + '_DEL.bed.gz')
                    f = gzip.open(outfile, 'wb')
                    try:
                        for l in set(lines):
                            f.write(l)
                    finally:
                        f.close()

                    bedfile = os.path.join(outdir, sample_name + '_' + \
                                           label_name + '_' + str(iteration) + '_DEL.sorted.bed')
                    os.system('sortBed -i ' + outfile + ' > ' + bedfile)


def create_models_with_real_data():

    if CREATE_MODELS:

        # for sample_name in ['NA12878', 'PATIENT1', 'PATIENT2']:
        # for sample_name in ['PATIENT1', 'PATIENT2']:
        for sample_name in ['NA12878']:

            print(f'Considering sample {sample_name}')

            res_df = pd.DataFrame()

            data, chr_list = load_windows(sample_list=[sample_name],
                                          base_dir=base_dir)

            data = save_windows(data, sample_name)

            # data = load_windows_from_zip(sample_name)

            chr_list = get_chr_list()

            print('Loaded data. Shape:')
            print(data.shape)

            labels_pickle = load_labels(sample_name)

            np.random.seed = 321

            for key in labels_pickle:

                if key != 'id':

                    print(f'Considering labels {key}')

                    labels = []
                    for chr in chr_list:
                        labels.extend(labels_pickle[key][chr])
                    print('%d labels' % len(labels))

                    data_filtered, labels_filtered = remove_label_from_data(data, labels, 'UK')
                    # labels_filtered = relabel_deletions(labels_filtered)

                    print('Shape filtered:')
                    print(data_filtered.shape)

                    # labels_filtered = np.array(labels)
                    # data_filtered = data

                    # random sample?
                    data_balanced, labels_balanced = balance_data(data_filtered, labels_filtered)

                    print('Shape balanced:')
                    print(data_filtered.shape)

                    X, y = split_data(data_balanced, labels_balanced)

                    print('Shape split:')
                    print(X['train'].shape)

                    y_binary, classlabels = get_categorical_labels(y)
                    num_classes = y_binary['train'].shape[1]

                    # models = modelgen.generate_models(X['train'].shape,
                    #                                   num_classes,
                    #                                   number_of_models=1,
                    #                                   model_type='CNN',
                    #                                   cnn_min_layers=1,
                    #                                   cnn_max_layers=1,
                    #                                   cnn_min_fc_nodes=6,
                    #                                   cnn_max_fc_nodes=6,
                    #                                   low_lr=4, high_lr=4,
                    #                                   metrics=['accuracy'])  # , precision, recall, f1])

                    models = modelgen.generate_models(X['train'].shape,
                                                      num_classes,
                                                      number_of_models=1,
                                                      model_type='CNN',
                                                      cnn_min_layers=2,
                                                      cnn_max_layers=2,
                                                      cnn_min_filters=4,
                                                      cnn_max_filters=4,
                                                      cnn_min_fc_nodes=6,
                                                      cnn_max_fc_nodes=6,
                                                      low_lr=2, high_lr=2)

                    for i in range(1, 11, 1):
                        model_path, results = run_model(X, y_binary, label_name=key,
                                                        sample_name=sample_name,
                                                        iteration=i,
                                                        classlabels=classlabels,
                                                        models = models)
                        # print(results)
                        res_df = res_df.append(results)
                        # print(model_path)

                    print(res_df)
                    res_df.to_csv('/'.join((base_dir, sample_name + '_performance_results.csv')))

                    test_model(X, y_binary, model_path)

    else:

        predict_windows_with_model()


def create_models_with_artificial_data():

    res_df = pd.DataFrame()
    sample_list = ['G1', 'N1']
    data, chr_list = load_windows(sample_list=sample_list,
                                  base_dir=base_dir_art)

    data = transpose_single(data)

    # data = save_windows(data, sample_name)
    # data = load_windows_from_zip(sample_name)

    chr_list = get_chr_list()

    print('Loaded data. Shape:')
    print(data.shape)

    labels = load_artificial_labels(sample_list, base_dir_art)

    data, labels = remove_label_from_data(data, labels, 'INS_pos')

    np.random.seed = 321

    data, labels = balance_data(data, labels)

    X, y = split_data(data, labels)

    print('Shape split:')
    print(X['train'].shape)

    y_binary, classlabels = get_categorical_labels(y)
    num_classes = y_binary['train'].shape[1]

    models = modelgen.generate_models(X['train'].shape,
                                      num_classes,
                                      number_of_models=1,
                                      model_type='CNN',
                                      cnn_min_layers=2,
                                      cnn_max_layers=2,
                                      cnn_min_filters=4,
                                      cnn_max_filters=4,
                                      cnn_min_fc_nodes=6,
                                      cnn_max_fc_nodes=6,
                                      low_lr=2, high_lr=2)

    for i in range(1, 11, 1):
        model_path, results = run_model(X, y_binary, label_name='label',
                                        sample_name='_'.join(sample_list),
                                        iteration=i,
                                        classlabels=classlabels,
                                        models=models)
        # print(results)
        res_df = res_df.append(results)
        # print(model_path)

    print(res_df)
    res_df.to_csv('/'.join((base_dir_art, 'artData_performance_results.csv')))

    test_model(X, y_binary, model_path)


def main():

    if ART_MODE:
        create_models_with_artificial_data()
    else:
        create_models_with_real_data()


if __name__ == '__main__':
    main()
