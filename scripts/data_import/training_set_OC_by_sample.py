import os
import numpy as np
import gzip
from collections import Counter
import logging
import argparse

chr_list = list(map(str, np.arange(1, 23)))
chr_list.extend(['X', 'Y'])

date = '070119'
label_type = 'bpi'


def data(sample_name):

    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
               'intermediate_data/' + sample_name
    data_output_file = base_dir + '_data.npy'
    label_output_file = base_dir + '_labels.npy'
    id_output_file = base_dir + '_ids.npy'

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

    data.extend(partial_data)
    labels.extend(partial_labels)
    ids.extend(partial_id)

    logging.info(Counter(labels))
    assert len(data) == len(labels)

    training_data = np.array(data)
    training_labels = np.array(labels)
    training_id = np.array(ids)

    data_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/data.npy'
    np.save(data_output_file, training_data)
    os.system('gzip -f ' + data_output_file)

    label_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/labels.npy'
    np.save(label_output_file, training_labels)
    os.system('gzip -f ' + label_output_file)

    id_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/ids.npy'
    np.save(id_output_file, training_id)
    os.system('gzip -f ' + id_output_file)


def combine_data():

    data = []
    labels = []
    ids = []

    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + date + '/'
    comparisons = os.listdir(base_dir)

    for sample_name in comparisons:

        base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
                   'intermediate_data/' + sample_name
        data_output_file = base_dir + '_data.npy'
        label_output_file = base_dir + '_labels.npy'
        id_output_file = base_dir + '_ids.npy'

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

        data_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/data.npy'
        np.save(data_output_file, training_data)
        os.system('gzip -f ' + data_output_file)

        label_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/labels.npy'
        np.save(label_output_file, training_labels)
        os.system('gzip -f ' + label_output_file)

        id_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/ids.npy'
        np.save(id_output_file, training_id)
        os.system('gzip -f ' + id_output_file)


def main():

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-s', '--sample', type=str, default='NA12878',
                        help="Specify sample")
    parser.add_argument('-l', '--logfile', default='channel_maker.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    # channel_list = get_channels()
    data(sample_name = args.sample)


if __name__ == '__main__':

    main()