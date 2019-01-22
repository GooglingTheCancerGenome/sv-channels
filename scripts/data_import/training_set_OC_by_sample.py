import os
import numpy as np
import gzip
from collections import Counter
import logging
import argparse
from keras.utils.np_utils import to_categorical


chr_list = list(map(str, np.arange(1, 23)))
chr_list.extend(['X', 'Y'])

date = '070119'
# label_type = 'bpi'


def transposeDataset(X):

    image = []
    for i in range (0, len(X -1)):
        tr = X[i].transpose()
        image.append(tr)
    return np.array(image)


def balance_data(training_data, training_labels, training_id):

    cnt_lab = Counter(training_labels)
    min_v = min([v for k, v in cnt_lab.items()])
    print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))

    data_balanced = []
    labels_balanced = []
    id_balanced = []

    for l in cnt_lab.keys():
        #print(l)
        iw = np.where(training_labels==l)
        ii = iw[0][:min_v]
        data_balanced.extend(training_data[ii])
        labels_balanced.extend(training_labels[ii])
        id_balanced.extend(training_id[ii])

    print(Counter(labels_balanced))

    X = np.array(data_balanced)
    y = np.array(labels_balanced)
    z = np.array(id_balanced)

    return X, y, z


def remove_nan(X, y, z):

    # Remove windows with nan if present
    # print(np.where(np.isnan(X)))
    idx = np.unique(np.where(np.isnan(X))[0])
    print(idx)
    # print(X[np.where(np.isnan(X))])
    print(z[idx])
    idx = np.unique(np.where(np.isnan(X))[0])
    X = np.delete(X, idx, 0)
    y = np.delete(y, idx, 0)
    z = np.delete(z, idx, 0)

    return X, y, z


def data(sample_name, label_type):

    base_dir = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
               'intermediate_data', label_type, sample_name)

    if not os.path.isdir(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    data_output_file = os.path.join(base_dir, 'data.npy')
    label_output_file = os.path.join(base_dir, 'labels.npy')
    id_output_file = os.path.join(base_dir, 'ids.npy')

    if not os.path.isfile(data_output_file + '.gz') and \
            not os.path.isfile(label_output_file + '.gz') and \
            not os.path.isfile(id_output_file + '.gz'):

        logging.info('Loading sample: %s...' % sample_name)
        datapath = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + \
                   date + '/' + sample_name

        partial_data = []
        partial_labels = []
        partial_id = []

        dico_file = datapath + '/MultiLabelData/labels.pickle.gz'
        with gzip.GzipFile(dico_file, "rb") as f:
            dico = np.load(f)
        f.close()

        logging.info(dico[label_type].keys())
        # logging.info(chr_list)

        # for i in dico[label_type].keys():
        for i in chr_list:
            if i in dico[label_type].keys():
                # if not (sample_name == 'O16_B16' and i == '2'):
                logging.info('Loading data for Chr%s' % i)

                labels_chr = dico[label_type][i]
                partial_labels.extend(labels_chr)
                id_chr = [d['chromosome'] + '_' + str(d['position']) for d in dico['id'][i]]
                partial_id.extend(id_chr)

                assert len(id_chr) == len(labels_chr)

                data_file = datapath + '/ChannelData/' + str(i) + '_channel_maker_real_germline.npy.gz'
                with gzip.GzipFile(data_file, "rb") as f:
                    data_mat = np.load(f)
                    logging.info(data_mat.shape)
                    assert data_mat.shape[0] == len(labels_chr)
                    partial_data.extend(data_mat)
                    del data_mat
                f.close()

        partial_data = np.array(partial_data)
        logging.info('partial data shape: %s' % str(partial_data.shape))
        partial_labels = np.array(partial_labels)
        logging.info('partial labels shape: %s' % str(partial_labels.shape))
        partial_id = np.array(partial_id)
        logging.info('partial ids shape: %s' % str(partial_id.shape))

        i_nosv = np.where(partial_labels == 'noSV')[0]
        logging.info('Number of noSV labels: %d' % len(i_nosv))
        # logging.info(i_nosv)

        i_nosv_idx = np.random.choice(a=i_nosv,
                                      # size=int(np.round(i_nosv.shape[0]/100)),
                                      size=100,
                                      replace=False)
        logging.info('Sampled noSV labels: %d' % len(i_nosv_idx))

        i_sv = np.where(partial_labels != 'noSV')[0]
        logging.info('Number of !noSV labels: %d' % len(i_sv))

        partial_data = np.append(partial_data[i_sv, :, :], partial_data[i_nosv_idx, :, :], axis=0)
        logging.info('partial data shape after append: %s' % str(partial_data.shape))
        partial_labels = np.append(partial_labels[i_sv], partial_labels[i_nosv_idx], axis=0)
        logging.info('partial labels shape after append: %s' % str(partial_labels.shape))
        partial_id = np.append(partial_id[i_sv], partial_id[i_nosv_idx], axis=0)
        logging.info('partial ids shape after append: %s' % str(partial_id.shape))

        intermediate_data_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
                                'intermediate_data/'
        if not os.path.exists(intermediate_data_dir):
            os.mkdir(intermediate_data_dir)

        # save intermediate data
        np.save(data_output_file, partial_data)
        os.system('gzip -f ' + data_output_file)

        np.save(label_output_file, partial_labels)
        os.system('gzip -f ' + label_output_file)

        np.save(id_output_file, partial_id)
        os.system('gzip -f ' + id_output_file)

    else:

        logging.info('Loading sample: %s...' % sample_name)

        with gzip.GzipFile(data_output_file + '.gz', "rb") as f:
            partial_data = np.load(f)
        with gzip.GzipFile(label_output_file + '.gz', "rb") as f:
            partial_labels = np.load(f)
        with gzip.GzipFile(id_output_file + '.gz', "rb") as f:
            partial_id = np.load(f)

    assert partial_data.shape[0] == len(partial_labels)
    assert len(partial_labels) == len(partial_id)


def combine_data(label_type):

    data = []
    labels = []
    ids = []

    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + date + '/'
    comparisons = os.listdir(base_dir)

    for sample_name in comparisons:

        base_dir = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
                                'intermediate_data', label_type, sample_name)

        data_output_file = os.path.join(base_dir, 'data.npy')
        label_output_file = os.path.join(base_dir, 'labels.npy')
        id_output_file = os.path.join(base_dir, 'ids.npy')

        logging.info('Loading sample: %s...' % sample_name)

        with gzip.GzipFile(data_output_file + '.gz', "rb") as f:
            partial_data = np.load(f)
        with gzip.GzipFile(label_output_file + '.gz', "rb") as f:
            partial_labels = np.load(f)
        with gzip.GzipFile(id_output_file + '.gz', "rb") as f:
            partial_id = np.load(f)

        logging.info(partial_data.shape)
        logging.info(partial_labels.shape)
        logging.info(partial_id.shape)

        data.extend(partial_data)
        labels.extend(partial_labels)
        ids.extend(partial_id)

        logging.info(Counter(labels))
        assert len(data) == len(labels)
        assert len(ids) == len(labels)

        training_data = np.array(data)
        training_labels = np.array(labels)
        training_id = np.array(ids)

        out_dir = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test', date, 'TestData', label_type)
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        data_output_file = os.path.join(out_dir, 'data.npy')
        np.save(data_output_file, training_data)
        os.system('gzip -f ' + data_output_file)

        label_output_file = os.path.join(out_dir, 'labels.npy')
        np.save(label_output_file, training_labels)
        os.system('gzip -f ' + label_output_file)

        id_output_file = os.path.join(out_dir, 'ids.npy')
        np.save(id_output_file, training_id)
        os.system('gzip -f ' + id_output_file)

        for l in ['UK', 'INS_start']:
            # Remove windows labelled as unknown ('UK') or INS_start (too few labels)
            keep = np.where(np.array(training_labels) != l)
            training_data = training_data[keep]
            training_labels = training_labels[keep]
            training_id = training_id[keep]

        X = np.array(training_data)
        y = np.array(training_labels)
        z = np.array(training_id)

        X = transposeDataset(X)

        # Derive mapclasses
        classes = sorted(list(set(y)))
        mapclasses = dict()
        for i, c in enumerate(classes):
            mapclasses[c] = i

        logging.info('Mapclasses: %s' % mapclasses)
        y_num = np.array([mapclasses[c] for c in y], dtype='int')
        y_binary = to_categorical(y_num)

        logging.info('X shape: %s' % X.shape)
        logging.info('y shape: %s' % y.shape)
        logging.info('y_binary shape: %s' % y_binary.shape)

        datapath_training = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test',
                                         date, 'TrainingData', label_type)

        if not os.path.isdir(datapath_training):
            os.makedirs(base_dir, datapath_training=True)

        data_output_file = datapath_training + 'OC' + '_' + label_type + '.npz'
        np.savez(data_output_file, X=X, y=y, y_binary=y_binary, ids=z)
        os.system('gzip -f ' + data_output_file)


def main():

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-s', '--sample', type=str, default='COMBINE',
                        help="Specify sample")
    parser.add_argument('-t', '--label_type', type=str, default='bpi',
                        help="Specify sample")
    parser.add_argument('-l', '--logfile', default='my.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    # channel_list = get_channels()
    sample_name = args.sample
    if sample_name == 'COMBINE':
        data(sample_name=sample_name, label_type=args.label_type)
    else:
        combine_data(label_type = args.label_type)


if __name__ == '__main__':

    main()