import os
import numpy as np
from keras.utils.np_utils import to_categorical
import keras
import gzip
from collections import Counter
import pandas as pd

# Auxiliary functions


def transposeDataset(X):
    image = []
    for i in range (0, len(X -1)):
        tr = X[i].transpose()
        image.append(tr)
    return np.array(image)


sample_name = 'NA12878'
#date = '231118'
date = '270219'
label_type = 'Mills2011_nanosv'


def real_data():
    # Create reference test set using Chr4 to ChrX

    def get_label_dict():
        # Load label dictionary
        dico_file = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                                 date, 'TestData_' + date, sample_name, 'MultiLabelData/labels.pickle.gz')
        with gzip.GzipFile(dico_file, "rb") as f:
            dico = np.load(f)
        f.close()

        return dico

    dico = get_label_dict()

    # Leaving out chromosome Y and MT for the moment
    chr_list = list(map(str, np.arange(4, 23)))
    chr_list.append('X')

    print(chr_list)

    training_data = []
    training_labels = []
    training_id = []

    datapath = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                             date, 'TestData_'+date, sample_name, 'ChannelData')

    for i in chr_list:
        print('Loading data for Chr%s' % i)
        data_file = datapath + sample_name + '_' + str(i) + '.npy.gz'
        with gzip.GzipFile(data_file, "rb") as f:
            data_mat = np.load(f)
            training_data.extend(data_mat)
        f.close()

        training_labels.extend(dico[label_type][i])
        training_id.extend([d['chromosome'] + '_' + str(d['position']) for d in dico['id'][i]])

    print(Counter(training_labels))

    training_data = np.array(training_data)
    training_labels = np.array(training_labels)
    training_id = np.array(training_id)

    assert len(training_data) == len(training_labels)

    return training_data, training_labels, training_id


def artificial_data():

    training_data = []
    training_labels = []

    base_dir = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                            date, 'TrainingData_'+date)
    sample = 'G1'

    for svtype in ['INDEL', 'INDEL_HOM']:

        datapath = os.path.join(base_dir, svtype, sample)
        data_file = os.path.join(datapath, 'ChannelData', sample+'.npy.gz')
        label_file = os.path.join(datapath, 'LabelData', sample+'_17_label.npy.gz')

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


def remove_label(training_data, training_labels, training_id, label = 'UK'):

    # Remove windows labelled as label
    keep = np.where(np.array(training_labels) != label)
    training_data = training_data[keep]
    training_labels = training_labels[keep]
    training_id = training_id[keep]

    return training_data, training_labels, training_id


def balance_dataset(training_data, training_labels, training_id):

    cnt_lab = Counter(training_labels)
    min_v = min([v for k, v in cnt_lab.items()])
    print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))

    data_balanced = []
    labels_balanced = []
    id_balanced = []

    for l in cnt_lab.keys():
        # print(l)
        iw = np.where(training_labels == l)
        ii = iw[0][:min_v]
        data_balanced.extend(training_data[ii])
        labels_balanced.extend(training_labels[ii])
        id_balanced.extend(training_id[ii])

    print(Counter(labels_balanced))

    X = np.array(data_balanced)
    y = np.array(labels_balanced)
    z = np.array(id_balanced)

    return X, y, z

