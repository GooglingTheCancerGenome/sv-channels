import argparse
import zarr
import gzip
import json
import numpy as np
from tensorflow.keras.models import load_model
import vcf


def predict(windows_list, labels_list, samples_list, model_names_list,
            models, manta_list, output_list):

    def load_windows(win_file, lab_file):
        X = zarr.load(win_file)
        with gzip.GzipFile(lab_file, 'r') as fin:
            y = json.loads(fin.read().decode('utf-8'))
        return X, y

    manta_dict = {s: m for s,m in zip(samples_list, manta_list)}
    output_dict = {s: m for s, m in zip(samples_list, output_list)}

    X = []
    y = []
    win_ids = []
    samples = []

    for w, l, s in zip(windows_list, labels_list, samples_list):

        partial_X, partial_y = load_windows(w, l)
        X.extend(partial_X)
        y.extend(partial_y.values())
        win_ids.extend(partial_y.keys())
        samples.extend([s]*len(partial_y))

    X = np.stack(X, axis=0)
    first_chrom = [w.split('/')[0] for w in win_ids]

    pos_dict = {}

    for s in set(samples):
        for c in model_names_list:

            chrom_idx = [i for i, k in enumerate(zip(first_chrom, samples)) if k[0] == c and k[1] == s]
            chrom_idx = np.asarray(chrom_idx)

            X_chrom = X[chrom_idx]
            win_ids_chrom = [win_ids[i] for i in chrom_idx]

            print('Predicting {} with model trained on {}...'.format(
                s, c))

            model = load_model(models[c])
            probs = model.predict(X_chrom, batch_size=100, verbose=True)
            predicted = probs.argmax(axis=1)

            n = 0
            for i, line in enumerate(win_ids_chrom):

                chrom1, pos1a, chrom2, pos2a, *_ = line.split('/')
                # TODO: not sure if we need probs[i][0] here.
                pos_dict[chrom1 + '_' + str(pos1a)] = str(probs[i][0] - probs[i][1])
                n += 1

            assert n == len(predicted), "number of variants and channel-sets should match"

            # print(pos_dict)

        reader = vcf.Reader(open(manta_dict[s], 'r'))
        writer = vcf.Writer(open(output_dict[s], 'w'), reader)

        for record in reader:
            k = record.CHROM + '_' + str(record.POS)
            if k in pos_dict.keys():
                record.QUAL = str(pos_dict[k])
                writer.write_record(record)


def main():

    parser = argparse.ArgumentParser(description='Use model(s) to predict')
    parser.add_argument('-m',
                        '--models',
                        type=str,
                        default='best_model1.keras,best_model2.keras',
                        help="Comma separated list of TensorFlow models in HDF5 format"
                        )
    parser.add_argument('-n',
                        '--model_names',
                        type=str,
                        default='chr1,chr2',
                        help="Comma separated list of chromosome names for the models (test chromosome)"
                        )
    parser.add_argument('-w',
                        '--windows',
                        type=str,
                        default='sv_chan.zarr,sv_chan.zarr',
                        help="Comma separated list of training data")
    parser.add_argument('-lab',
                        '--labels',
                        type=str,
                        default='labels/labels.json.gz,labels/labels.json.gz',
                        help="Comma separated list of JSON.GZ file for labels")
    parser.add_argument('-sm',
                        '--samples',
                        type=str,
                        default='SAMPLE1,SAMPLE2',
                        help="Comma separated list of sample names")
    parser.add_argument('-ma',
                        '--manta',
                        type=str,
                        default='manta1.vcf,manta2.vcf',
                        help="Comma separated list of input Manta VCF files")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='manta1_out.vcf,manta2_out.vcf',
                        help="Comma separated list of output VCF files")

    args = parser.parse_args()

    windows_list = args.windows.split(',')
    labels_list = args.labels.split(',')
    samples_list = args.samples.split(',')

    manta_list = args.manta.split(',')
    output_list = args.output.split(',')

    models_list = args.models.split(',')
    model_names_list = args.model_names.split(',')
    models = {k: v for k, v in zip(model_names_list, models_list)}

    predict(windows_list, labels_list, samples_list, model_names_list,
            models, manta_list, output_list)


if __name__ == '__main__':
    main()
