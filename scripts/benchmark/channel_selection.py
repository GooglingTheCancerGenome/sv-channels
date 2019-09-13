# Dependencies
import os
import numpy as np
from keras.utils.np_utils import to_categorical
import keras
import gzip
from collections import Counter
from sklearn.model_selection import StratifiedKFold

parameters = {
    'sample_name': 'NA12878',
    'date': '231118',
    'label_type':  'Mills2011_nanosv',
    'base_dir': '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/'
}


def get_channel_names(data):

    number_channels = data.shape[1]
    channel_name = ["None"] * number_channels

    channel_name[0] = "coverage"
    channel_name[1] = "#left clipped reads"
    channel_name[2] = "#right clipped reads"
    channel_name[3] = "INV_before"
    channel_name[4] = "INV_after"
    channel_name[5] = "DUP_before"
    channel_name[6] = "DUP_after"
    channel_name[7] = "TRA_opposite"
    channel_name[8] = "TRA_same"

    i = 9
    for direction in ['Forward', 'Reverse']:
        for clipped in ['Left', 'Right', 'Not']:
            for value in ['sum', 'num', 'median', 'outliers']:
                channel_name[i] = direction + '_' + clipped + '_Clipped_' + value
                i = i + 1

    channel_name[i] = "#left split reads"
    i = i + 1
    channel_name[i] = "#right split reads"
    i = i + 1

    for clipped in ['L', 'R']:
        for value in ['sum', 'num', 'median']:
            channel_name[i] = clipped + '_SplitRead_' + value
            i = i + 1

    channel_name[i] = "GC"
    i = i + 1
    channel_name[i] = "Mappability"
    i = i + 1
    channel_name[i] = "One_hot_Ncoding"

    return channel_name


def get_label_dict():

    dict_file = parameters['base_dir'] + \
                parameters['date']+'/TestData_'+parameters['date']+\
                '/'+parameters['sample_name']+'/MultiLabelData/labels.pickle.gz'
    with gzip.GzipFile(dict_file, "rb") as f:
        dict = np.load(f)
    f.close()
    return dict


def balance_dataset(data, labels):

    cnt_lab = Counter(labels)
    min_v = min([v for k, v in cnt_lab.items()])
    max_v = max([v for k, v in cnt_lab.items()])

    print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))
    print('Maximum number of labels = ' + str(max_v))

    data_balanced = []
    labels_balanced = []

    for l in cnt_lab.keys():

        iw = np.where(labels == l)
        ii = iw[0][:min_v]
        data_balanced.extend(data[ii])
        labels_balanced.extend(labels[ii])

    print(Counter(labels_balanced))

    X = np.array(data_balanced)
    y = np.array(labels_balanced)

    return X, y


def get_data(label_dict):

    chr_list = list(map(str, np.arange(4, 23)))
    chr_list.append('X')

    data = []
    labels = []

    label_type = parameters['label_type']
    datapath = parameters['base_dir'] + \
               parameters['date'] + '/TestData_' + parameters['date'] + '/' + \
               parameters['sample_name'] + '/ChannelData/'

    for c in chr_list:

        print('Loading data for Chr%s' % c)
        data_file = datapath + parameters['sample_name'] + '_' + str(c) + '.npy.gz'
        with gzip.GzipFile(data_file, "rb") as f:
            data_mat = np.load(f)
            data.extend(data_mat)
        f.close()

        labels.extend(label_dict[label_type][c])

    print(Counter(labels))
    assert len(data) == len(labels)

    data = np.array(data)
    labels = np.array(labels)

    # Remove unknowns
    keep = np.where(labels != 'UK')
    data = data[keep]
    labels = labels[keep]

    # Remove NANs
    idx = numpy.unique(np.where(np.isnan(data))[0])
    data = np.delete(data, idx, 0)
    labels = np.delete(labels, idx, 0)