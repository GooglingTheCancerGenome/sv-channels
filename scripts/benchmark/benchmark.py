# Imports

import numpy as np
from keras.utils.np_utils import to_categorical
import keras
import gzip
from collections import Counter

sample_name = 'G1'
context_dict = {'G1': 'train', 'N1': 'real'}
date = '230718' # '010618'

base_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/' + date + '/TrainingData_' + date


def load_data(sample_list, base_dir):

    # Load channel data and labels by chromosome

    # Leaving out chromosome Y and MT for the moment
    #chr_list = list(map(str, np.arange(1, 23)))
    #chr_list.append('X')

    chr_list = ['1']
    # print(chr_list)

    data = []
    labels = []

    for sample_name in sample_list:

        print('Considering sample %s' % sample_name)

        filename_suffix = '_channel_maker_' + context_dict[sample_name] + '_germline.npy.gz'

        for i in chr_list:

            print('Loading data and labels for Chr%s' % i)

            data_file = base_dir + '/' + sample_name + '/ChannelData/' + str(i) + filename_suffix

            with gzip.GzipFile(data_file, "rb") as f:
                data_mat = np.load(f)
                data.extend(data_mat)
            f.close()

            label_file = base_dir + '/' + sample_name + '/LabelData/' + sample_name + '_' + str(i) + '_label.npy.gz'

            with gzip.GzipFile(label_file, "rb") as f:
                data_lab = np.load(f)
                labels.extend(data_lab)
            f.close()

            # Do not consider insertions: remove INS_pos
            no_ins = np.where(np.array(labels) != 'INS_pos')
            data = data[no_ins]
            labels = labels[no_ins]

            # Relabel deletions
            del_idx = np.where(np.logical_or(np.array(labels) == 'DEL_start', np.array(labels) == 'DEL_end'))
            print(del_idx)
            print(len(del_idx))
            print(len(labels))
            labels[del_idx] = 'DEL'

            # print('Length of data:%d, label length:%d' % (len(data_mat), len(data_lab)))
            assert len(data_mat) == len(data_lab)

    print(Counter(labels))
    assert len(data) == len(labels)

    return data, labels


def balance_data(data, labels):

    cnt_lab = Counter(labels)
    min_v = min([v for k, v in cnt_lab.items()])

    print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))

    # Balance data

    data_balanced = []
    labels_balanced = []

    for l in cnt_lab.keys():
        print(l)
        iw = np.where(np.array(labels) == l)
        # print(iw[0])
        ii = np.random.choice(a=iw[0], size=min_v, replace=False)
        data_balanced.extend(data[ii])
        labels_balanced.extend(labels[ii])

    print(Counter(labels_balanced))

    X = np.array(data_balanced)
    y = np.array(labels_balanced)

    del data
    del labels
    del data_balanced
    del labels_balanced

    return X, y


def split_data(data, labels)

    # Split into training, validation and test set 60/20/20
    cnt_lab = Counter(labels)
    print(cnt_lab)
    n_lab = [v for v in cnt_lab.values()][0]
    print(n_lab)

    i_train = int(n_lab * 0.6)
    i_val = i_train + int(n_lab * 0.2)
    print(i_train)
    print(i_val)

    X_train = []
    y_train = []
    X_val = []
    y_val = []
    X_test = []
    y_test = []

    for l in cnt_lab.keys():
        iw = np.where(y == l)
        # print(iw[0])
        # print(iw[0][:i_train])
        X_train.extend(X[iw[0][:i_train]])
        X_val.extend(X[iw[0][i_train:i_val]])
        X_test.extend(X[iw[0][i_val:]])

        y_train.extend(y[iw[0][:i_train]])
        y_val.extend(y[iw[0][i_train:i_val]])
        y_test.extend(y[iw[0][i_val:]])

    X_train = np.array(X_train)
    X_val = np.array(X_val)
    X_test = np.array(X_test)

    y_train = np.array(y_train)
    y_val = np.array(y_val)
    y_test = np.array(y_test)

    print(X_train.shape)
    print(X_val.shape)
    print(X_test.shape)

    print(y_train.shape)
    print(y_val.shape)
    print(y_test.shape)


def main():

    data, labels =load_data(sample_list=['G1', 'N1'], base_dir=base_dir)


if __name__ == '__main__':
    main()