import argparse
import os
import subprocess

import pandas as pd
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import to_categorical

from model_functions import evaluate_model, get_data


def predict(input_data, input_labels, sample_name, svtype, model_fn, model_name, output_dir):

    os.makedirs(output_dir, exist_ok=True)
    X, y, win_ids = get_data(input_data, input_labels, svtype)
    y_binary = to_categorical(y, num_classes=params['n_classes'])
    print('Predicting {} with model trained on {}...'.format(
        sample_name, model_name))
    model = load_model(model_fn)
    results = pd.DataFrame()

    intermediate_results, metrics = evaluate_model(model,
                                                   X, y_binary, win_ids,
                                                   results,
                                                   params['mapclasses'],
                                                   output_dir, svtype)

    results = results.append(intermediate_results)
    results.to_csv(os.path.join(output_dir, 'metrics.csv'), sep='\t')


def main():

    parser = argparse.ArgumentParser(description='Use model to predict')
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='best_model.keras',
                        help="TensorFlow model in HDF5 format"
                        )
    parser.add_argument('-n',
                        '--model_name',
                        type=str,
                        default='manta_model',
                        help="Name of the model"
                        )
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default='sv_chan.zarr',
                        help="Specify list of windows"
                        )
    parser.add_argument('-lab',
                        '--labels',
                        type=str,
                        default='labels/labels.json.gz',
                        help="Specify list of labels"
                        )
    parser.add_argument('-sn',
                        '--sample_name',
                        type=str,
                        default='test',
                        help="Specify sample name"
                        )
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

    windows_list = args.input.split(',')
    labels_list = args.labels.split(',')

    predict(windows_list, labels_list, args.sample_name, args.svtype,
            args.model, args.model_name, args.output)

    out_prefix = os.path.join(args.output, "sv-channels")

    merge_sv_calls = ' '.join([
        "Rscript", os.path.join(args.sv_channels, "scripts", "R", "merge_sv_calls.R"),
        "-i", os.path.join("../../svchannels", args.output), # BSP-TODO: what do do here?
        "-f", args.encode_blacklist,
        "-n", args.n_regions,
        "-m split_reads",
        "-o", os.path.join(args.output, "sv-channels")
    ])

    cmd_out = subprocess.run(merge_sv_calls, shell=True, check=True)

    assert os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py")
    assert os.path.join(args.sv_channels, "scripts/genome_wide", args.output + '.bedpe')

    bedpe_to_vcf = ' '.join([
        os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py"),
        "-i", os.path.join(args.output, "sv-channels.DEL.bedpe"),
        "-b", args.twobit,
        "-s", args.sample_name,
        "-o", os.path.join(args.output, "sv-channels.DEL.vcf")
    ])
    print(bedpe_to_vcf)
    cmd_out = subprocess.run(bedpe_to_vcf, shell=True, check=True)
    print(cmd_out)


if __name__ == '__main__':
    main()
