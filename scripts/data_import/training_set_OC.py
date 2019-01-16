import os
import numpy as np
import gzip
from collections import Counter


def main():

    chr_list = list(map(str, np.arange(1, 23)))
    chr_list.extend(['X', 'Y'])

    date = '070119'
    label_type = 'bpi'
    base_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + date + '/'
    comparisons = os.listdir(base_dir)

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
                labels.extend(dico[label_type][chrom_name])
                ids.extend(dico['id'][chrom_name])

    labels = np.array(labels)
    ids = np.array(ids)

    print(Counter(labels))

    label_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/all_labels.npy'
    np.save(label_output_file, labels)
    os.system('gzip ' + label_output_file)

    id_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/all_ids.npy'
    np.save(id_output_file, ids)
    os.system('gzip ' + id_output_file)

    data = []
    labels = []
    ids = []

    for sample_name in comparisons:

        print('Loading sample: %s...' % sample_name)
        datapath = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData_' + date + '/' + sample_name

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
                if not (sample_name == 'O16_B16' and i == '2'):
                    print('Loading data for Chr%s' % i)

                    partial_labels.extend(dico[label_type][i])
                    partial_id.extend([d['chromosome'] + '_' + str(d['position']) for d in dico['id'][i]])

                    data_file = datapath + '/ChannelData/' + sample_name + '_' + str(i) + '.npy.gz'
                    with gzip.GzipFile(data_file, "rb") as f:
                        data_mat = np.load(f)
                        partial_data.extend(data_mat)
                    f.close()

        partial_labels = np.array(partial_labels)
        i_nosv = np.where(partial_labels == 'noSV')[0]

        # print(i_nosv)

        i_nosv_idx = np.random.choice(a=i_nosv,
                                      # size=int(np.round(i_nosv.shape[0]/100)),
                                      size=100,
                                      replace=False)
        i_sv = np.where(partial_labels != 'noSV')[0]

        partial_data = np.array(partial_data)
        partial_data = np.append(partial_data[i_sv], partial_data[i_nosv_idx])

        partial_labels = np.array(partial_labels)
        partial_data = np.append(partial_labels[i_sv], partial_labels[i_nosv_idx])

        partial_id = np.array(partial_id)
        partial_id = np.append(partial_id[i_sv], partial_id[i_nosv_idx])

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
    os.system('gzip ' + data_output_file)

    label_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/labels.npy'
    np.save(label_output_file, training_labels)
    os.system('gzip ' + label_output_file)

    id_output_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/' + date + '/TestData/ids.npy'
    np.save(id_output_file, training_id)
    os.system('gzip ' + id_output_file)


if __name__ == '__main__':
    main()