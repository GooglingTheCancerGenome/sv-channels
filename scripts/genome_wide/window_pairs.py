import os
import numpy as np
from keras.utils.np_utils import to_categorical
import gzip
from collections import Counter, defaultdict
import argparse
from time import time
import bz2file
import pickle
import logging
from itertools import chain
import json

with open('parameters.json', 'r') as f:
    config = json.load(f)

HPC_MODE = config["DEFAULT"]["HPC_MODE"]
# Window size
win_hlen = config["DEFAULT"]["WIN_HLEN"]
win_len = config["DEFAULT"]["WIN_HLEN"] * 2

date = '060219'

chr_list = list(map(str, np.arange(1, 23)))
chr_list.append('X')


def get_channel_dir(sample_name):
    if HPC_MODE:

        datapath_prefix = '/hpc/cog_bioinf/ridder/users/lsantuari'
        channel_dir = datapath_prefix + '/Git/DeepSV_runs/' + date + '/CNN/scripts'
    else:
        channel_dir = '/Users/lsantuari/Documents/Processed/Test/test_060219'

    return channel_dir


def transposeDataset(X):
    image = []
    for i in range(0, len(X - 1)):
        tr = X[i].transpose()
        image.append(tr)
    return np.array(image)


def load_labels(sample_name):
    # Load label dictionary
    labels_file = os.path.join(get_channel_dir(sample_name), sample_name, 'label_npy_win'+win_len,
                               'labels.pickle.gz')
    with gzip.GzipFile(labels_file, "rb") as f:
        labels_dict = np.load(f)
    f.close()
    return labels_dict


def load_split_read_positions(sample_name):
    positions = dict()
    locations = dict()

    vec_type = 'split_read_pos'

    for chrName in chr_list:
        logging.info('Loading SR positions for Chr%s' % chrName)
        # Load files

        fn = os.path.join(get_channel_dir(sample_name), sample_name, vec_type, chrName + '_' + vec_type + '.pbz2')

        with bz2file.BZ2File(fn, 'rb') as f:
            positions[chrName], locations[chrName] = pickle.load(f)

    return positions, locations


def load_windows(sample_name, labels_dict):
    training_data = []
    training_labels = defaultdict(list)
    training_id = []

    for i in chr_list:
        logging.info('Loading data for Chr%s' % i)
        data_file = os.path.join(get_channel_dir(sample_name), sample_name, 'channel_maker_real_germline_win'+win_len,
                                 sample_name + '_' + str(i) + '.npy.gz')
        with gzip.GzipFile(data_file, "rb") as f:
            data_mat = np.load(f)
            for label_type in labels_dict.keys():
                assert len(data_mat) == len(labels_dict[label_type][i])
            training_data.extend(data_mat)
        f.close()
        for label_type in labels_dict.keys():
            training_labels[label_type].extend(labels_dict[label_type][i])
        training_id.extend([d for d in labels_dict['id'][i]])

    for k in training_labels.keys():
        logging.info('Labels: {}'.format(k))
        logging.info(Counter(training_labels[k]))
        assert len(training_data) == len(training_labels[k])

    training_data = np.array(training_data)
    training_id = np.array(training_id)

    return training_data, training_labels, training_id


def save_window_pairs(sample_name, X, y, y_binary, z):
    data_output_file = os.path.join(get_channel_dir(sample_name), '_'.join([sample_name, 'pairs']))
    np.savez_compressed(data_output_file, X=X)
    np.savez_compressed(data_output_file + '_labels', y=y, y_binary=y_binary, z=z)
    # os.system('gzip -f ' + data_output_file + '.npz')


def load_window_pairs(sample_name):
    data_output_file = os.path.join(get_channel_dir(sample_name), '_'.join([sample_name, 'pairs']))

    with data_output_file + '.npz' as f:
        npzfiles = np.load(f)
        X = npzfiles['X']

    with data_output_file + '_labels.npz' as f:
        npzfiles = np.load(f)
        y = npzfiles['y']
        y_binary = npzfiles['y_binary']
        z = npzfiles['z']

    print(X.shape)
    for k in y.keys():
        print(y[k].shape)

    return X, y, y_binary, z


def make_window_pairs(sample_name):

    labels_dict = load_labels(sample_name)
    training_data, training_labels, training_id = load_windows(sample_name, labels_dict)

    positions, locations = load_split_read_positions(sample_name)
    list_of_locations = list(chain.from_iterable(locations.values()))

    win_id_dict = {value['chromosome'] + '_' + str(value['position']): counter
                   for counter, value in enumerate(training_id)}

    lab_dict = dict()
    for i, k in enumerate(training_labels.keys()):
        lab_dict[k] = {value['chromosome'] + '_' + str(value['position']): l
                       for value, l in zip(training_id, training_labels[k])}

    padding = np.zeros(shape=(training_data.shape[1], 10), dtype=np.uint32)

    training_data_pairs = []
    training_labels_pairs = defaultdict(list)
    training_pos_pairs = []

    for loc in list_of_locations:

        chr1, pos1, chr2, pos2 = loc
        bp1 = chr1 + '_' + str(pos1)
        bp2 = chr2 + '_' + str(pos2)

        if bp1 in win_id_dict.keys() and bp2 in win_id_dict.keys():
            # print('%s -- %s' % (bp1, bp2))
            training_pos_pairs.append(bp1 + ':' + bp2)

            # windows side by side
            np_win_pair = np.concatenate((training_data[win_id_dict[bp1]],
                                          padding,
                                          training_data[win_id_dict[bp2]]),
                                         axis=1)

            # windows on top of one another
            # np_win_pair = np.concatenate((training_data[win_id_dict[bp1]],
            #                                  training_data[win_id_dict[bp2]]),
            #                                 axis=0)

            # print(np_win_pair.shape)
            training_data_pairs.append(np_win_pair)

            for k in lab_dict.keys():
                if lab_dict[k][bp1] == 'DEL_start' and lab_dict[k][bp2] == 'DEL_end':
                    training_labels_pairs[k].append('DEL')
                else:
                    training_labels_pairs[k].append('noSV')

    training_data_pairs = np.array(training_data_pairs)

    training_pos_pairs = np.array(training_pos_pairs)

    logging.info(training_data_pairs.shape)
    for k in training_labels_pairs.keys():
        logging.info('Label set:{}, counter:{}'.format(k, Counter(training_labels_pairs[k])))
    logging.info(training_pos_pairs.shape)

    logging.info('Transposing...')
    X = transposeDataset(training_data_pairs)
    z = training_pos_pairs

    mapclasses = {'DEL': 0, 'noSV': 1}

    y = training_labels_pairs
    y_binary = dict()
    for k in y.keys():
        logging.info(Counter(y[k]))
        y_num = np.array([mapclasses[c] for c in y[k]], dtype='int')
        y_binary[k] = to_categorical(y_num)

    logging.info('Saving window pairs...')
    save_window_pairs(sample_name, X, y, y_binary, z)


def main():
    # Parse the arguments of the script
    parser = argparse.ArgumentParser(description='Make window pairs')
    parser.add_argument('-s', '--samplename', type=str, default='NA12878',
                        help="Specify sample name")
    parser.add_argument('-l', '--logfile', default='window_pairs.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    make_window_pairs(sample_name=args.samplename)
    tn = (time() - t0)
    logging.info('Time: split read positions on sample %s and label set %s: %f seconds, %f minutes' % (
        args.samplename, args.labeltype, tn, tn / 60))


if __name__ == '__main__':
    main()
