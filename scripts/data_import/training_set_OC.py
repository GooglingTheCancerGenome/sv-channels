import os
import numpy as np
import gzip
from collections import Counter


# def get_channels():
#
#     channel_list = np.append(np.arange(0, 9),
#                              [33, 34])
#     channel_list = np.append(channel_list,
#                              [12, 16, 20, 24, 28, 32])
#
#     channel_list = np.append(channel_list,
#                              np.arange(41, 50))
#     channel_list = np.append(channel_list,
#                              [74, 75])
#     channel_list = np.append(channel_list,
#                              [53, 57, 61, 65, 69, 73])
#
#     channel_list = np.append(channel_list,
#                              np.arange(82, 85))
#
#     print(len(channel_list))
#
#     return channel_list


def main():

    # channel_list = get_channels()

    chr_list = list(map(str, np.arange(1, 23)))
    chr_list.extend(['X', 'Y'])

    date = '070119'
    label_type = 'bpi'
    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + date + '/'
    comparisons = os.listdir(base_dir)

    label_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/all_labels.npy'

    if not os.path.isfile(label_output_file+'.gz'):

        ids = []
        labels = []
        for sample_name in comparisons:
            print('Loading %s labels...' % sample_name)
            dico_file = base_dir + sample_name + '/MultiLabelData/labels.pickle.gz'
            with gzip.GzipFile(dico_file, "rb") as f:
                dico = np.load(f)
            f.close()
            for chrom_name in chr_list:
                if chrom_name in dico[label_type].keys():
                    assert len(dico[label_type][chrom_name]) == len(dico['id'][chrom_name])
                    print('Chr%s:%d labels' % (chrom_name, len(dico[label_type][chrom_name])))
                    labels.extend(dico[label_type][chrom_name])
                    ids.extend(dico['id'][chrom_name])

        labels = np.array(labels)
        ids = np.array(ids)

        print(Counter(labels))

        np.save(label_output_file, labels)
        os.system('gzip -f ' + label_output_file)

        id_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/all_ids.npy'
        np.save(id_output_file, ids)
        os.system('gzip -f ' + id_output_file)

    data = []
    labels = []
    ids = []

    for sample_name in comparisons:

        base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' + \
                           'intermediate_data/' + sample_name
        data_output_file =  base_dir + '_data.npy'
        label_output_file = base_dir + '_labels.npy'
        id_output_file = base_dir + '_ids.npy'

        if not os.path.isfile(data_output_file + '.gz') and \
                not os.path.isfile(label_output_file + '.gz') and \
                not os.path.isfile(id_output_file + '.gz'):

            print('Loading sample: %s...' % sample_name)
            datapath = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + \
                       date + '/' + sample_name

            partial_data = []
            partial_labels = []
            partial_id = []

            dico_file = datapath + '/MultiLabelData/labels.pickle.gz'
            with gzip.GzipFile(dico_file, "rb") as f:
                dico = np.load(f)
            f.close()

            print(dico[label_type].keys())
            # print(chr_list)

            # for i in dico[label_type].keys():
            for i in chr_list:
                if i in dico[label_type].keys():
                    # if not (sample_name == 'O16_B16' and i == '2'):
                        print('Loading data for Chr%s' % i)

                        labels_chr = dico[label_type][i]
                        partial_labels.extend(labels_chr)
                        id_chr = [d['chromosome'] + '_' + str(d['position']) for d in dico['id'][i]]
                        partial_id.extend(id_chr)

                        assert len(id_chr) == len(labels_chr)

                        data_file = datapath + '/ChannelData/' + str(i) + '_channel_maker_real_germline.npy.gz'
                        with gzip.GzipFile(data_file, "rb") as f:
                            data_mat = np.load(f)
                            print(data_mat.shape)
                            assert data_mat.shape[0] == len(labels_chr)
                            partial_data.extend(data_mat)
                            del data_mat
                        f.close()

            partial_data = np.array(partial_data)
            print('partial data shape: %s' % partial_data.shape)
            partial_labels = np.array(partial_labels)
            print('partial labels shape: %s' % partial_labels.shape)
            partial_id = np.array(partial_id)
            print('partial ids shape: %s' % partial_id.shape)

            i_nosv = np.where(partial_labels == 'noSV')[0]
            print('Number of noSV labels: %d' % len(i_nosv))
            # print(i_nosv)

            i_nosv_idx = np.random.choice(a=i_nosv,
                                          # size=int(np.round(i_nosv.shape[0]/100)),
                                          size=100,
                                          replace=False)
            print('Sampled noSV labels: %d' % len(i_nosv_idx))

            i_sv = np.where(partial_labels != 'noSV')[0]
            print('Number of !noSV labels: %d' % len(i_nosv))

            partial_data = np.append(partial_data[i_sv,:,:], partial_data[i_nosv_idx,:,:], axis = 0)
            print('partial data shape after append: %s' % partial_data.shape)
            partial_labels = np.append(partial_labels[i_sv], partial_labels[i_nosv_idx], axis = 0)
            print('partial labels shape after append: %s' % partial_labels.shape)
            partial_id = np.append(partial_id[i_sv], partial_id[i_nosv_idx], axis = 0)
            print('partial ids shape after append: %s' % partial_id.shape)

            intermediate_data_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/' +\
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

            print('Loading sample: %s...' % sample_name)

            with gzip.GzipFile(data_output_file+'.gz', "rb") as f:
                partial_data = np.load(f)
            with gzip.GzipFile(label_output_file+'.gz', "rb") as f:
                partial_labels = np.load(f)
            with gzip.GzipFile(id_output_file+'.gz', "rb") as f:
                partial_id = np.load(f)

            print(partial_data.shape)
            print(partial_labels.shape)
            print(partial_id.shape)

            assert partial_data.shape[0] == len(partial_labels)
            assert len(partial_labels) == len(partial_id)

        data.extend(partial_data)
        labels.extend(partial_labels)
        ids.extend(partial_id)

    print(Counter(labels))
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


if __name__ == '__main__':

    main()