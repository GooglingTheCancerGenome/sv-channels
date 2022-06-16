import argparse
import os
import subprocess
import zarr
import gzip
import json
import numpy as np
import pandas as pd
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import to_categorical

from model_functions import evaluate_model, get_data


def predict(windows_list, labels_list, samples_list, model_names_list, svtype, models, output_dir):

    def load_windows(win_file, lab_file):
        X = zarr.load(win_file)
        with gzip.GzipFile(lab_file, 'r') as fin:
            y = json.loads(fin.read().decode('utf-8'))
        return X, y

    os.makedirs(output_dir, exist_ok=True)

    mapclasses = {svtype: 0, 'no' + svtype: 1}

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

    # X, y, win_ids = get_data(input_data, input_labels, svtype)

    y = np.array([mapclasses[i] for i in y])
    y_binary = to_categorical(y, num_classes=params['n_classes'])

    first_chrom = [w.split('/')[0] for w in win_ids]

    for s in set(samples):
        for c in model_names_list:

            chrom_idx = [i for i, k in enumerate(zip(first_chrom, samples)) if k[0] == c and k[1] == s]
            chrom_idx = np.asarray(chrom_idx)

            X_chrom = X[chrom_idx]
            y_binary_chrom = y_binary[chrom_idx]
            win_ids_chrom = [win_ids[i] for i in chrom_idx]

            print('Predicting {} with model trained on {}...'.format(
                s, c))

            model = load_model(models[c])

            results = pd.DataFrame()

            intermediate_results, metrics = evaluate_model(model,
                                                           X_chrom, y_binary_chrom, win_ids_chrom,
                                                           results,
                                                           params['mapclasses'],
                                                           output_dir+'_'+s+'_'+c, svtype)

            results = results.append(intermediate_results)
            results.to_csv(os.path.join(output_dir+'_'+s+'_'+c, 'metrics.csv'), sep='\t')


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
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    parser.add_argument('-fe',
                        '--encode_blacklist',
                        type=str,
                        default='../../data/ENCFF001TDO.bed',
                        help="ENCODE blacklist")
    parser.add_argument('-fn',
                        '--n_regions',
                        type=str,
                        default='../../data/reference_N_regions.bed',
                        help="Regions in the genome containing Ns")
    parser.add_argument('-tb',
                        '--twobit',
                        type=str,
                        default='../data/test.2bit',
                        help="TwoBit reference genome")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='results',
                        help="Output folder")
    parser.add_argument('-svc',
                        '--sv_channels',
                        type=str,
                        default='../',
                        help="sv-channels folder")

    args = parser.parse_args()

    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    # Parameters
    global params
    params = {
        'mapclasses': mapclasses,
        'n_classes': len(mapclasses.keys())
    }

    windows_list = args.windows.split(',')
    labels_list = args.labels.split(',')
    samples_list = args.samples.split(',')

    models_list = args.models.split(',')
    model_names_list = args.model_names.split(',')
    models = {k: v for k, v in zip(model_names_list, models_list)}

    predict(windows_list, labels_list, samples_list, model_names_list, args.svtype,
            models, args.output)

    for s in samples_list:

        os.makedirs(args.output + '_' + s + '/predictions', exist_ok=True)

        concat_chroms = 'cat ' + args.output + '_' + s + '_*/predictions/correct.bedpe > ' + args.output + '_' + s + \
                        '/predictions/correct.bedpe'

        print(concat_chroms)
        cmd_out = subprocess.run(concat_chroms, shell=True, check=True)
        print(cmd_out)

        # Merge calls
        merge_sv_calls = ' '.join([
            "cd ", os.path.join(args.sv_channels, "scripts/R") + "; ",
            "Rscript merge_sv_calls.R",
            "-i", os.path.join(args.output + '_' + s),
            "-f", args.encode_blacklist,
            "-n", args.n_regions,
            "-m split_reads",
            "-o", os.path.join(args.output + '_' + s, "sv-channels" + "." + s)
        ])

        print(merge_sv_calls)
        cmd_out = subprocess.run(merge_sv_calls, shell=True, check=True)
        print(cmd_out)

        assert os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py")
        assert os.path.join(args.sv_channels, "scripts/genome_wide", args.output + '.bedpe')

        # Convert the BEDPE output into VCF format
        bedpe_to_vcf = ' '.join([
            "conda activate sv-channels; pwd; python ",
            os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py"),
            "-i", os.path.join(args.output + '_' + s, "sv-channels" + "." + s + ".DEL.bedpe"),
            "-b", args.twobit,
            "-s", s,
            "-o", os.path.join(args.output + '_' + s, "sv-channels" + "." + s + ".DEL.vcf")
        ])
        print(bedpe_to_vcf)
        cmd_out = subprocess.run(bedpe_to_vcf, shell=True, check=True)
        print(cmd_out)


if __name__ == '__main__':
    main()
