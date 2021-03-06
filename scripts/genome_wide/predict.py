import argparse
import os
import subprocess

import pandas as pd
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import to_categorical

from model_functions import evaluate_model, get_data


def predict(input_data, sample_name, svtype, model_fn, model_name, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    X, y, win_ids = get_data(input_data, True, svtype)
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
    default_win = 200
    default_path = os.path.join('./cnn/win'+str(default_win), 'split_reads')
    def_windows_file = os.path.join(
        default_path, 'windows', 'DEL', 'windows_en.npz')
    parser = argparse.ArgumentParser(description='Use model to predict')
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='./cnn/win200/split_reads/cv/DEL/1/model.hdf5',
                        help="TensorFlow model in HDF5 format"
                        )
    parser.add_argument('-n',
                        '--model_name',
                        type=str,
                        default='cnn',
                        help="Name of the model"
                        )
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default=def_windows_file,
                        help="Specify list of windows"
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
                        default='../../data/test.2bit',
                        help="TwoBit reference genome")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='results',
                        help="Output folder")
    args = parser.parse_args()
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    # Parameters
    global params
    params = {
        'mapclasses': mapclasses,
        'n_classes': len(mapclasses.keys())
    }
    windows_list = args.input.split(',')
    predict(windows_list, args.sample_name, args.svtype,
            args.model, args.model_name, os.path.join(args.output, args.svtype))
    out_prefix = os.path.join(args.output, "sv-channels")
    merge_sv_calls = ' '.join([
        "cd ../R; "
        "Rscript merge_sv_calls.R",
        "-i", os.path.join("../genome_wide", args.output),
        "-f", args.encode_blacklist,
        "-n", args.n_regions,
        "-m split_reads",
        "-o", os.path.join("../genome_wide", out_prefix)
    ])
    print(merge_sv_calls)
    cmd_out = subprocess.run(merge_sv_calls, shell=True, check=True)
    print(cmd_out)

    assert os.path.join("../utils/bedpe_to_vcf.py")
    assert os.path.join("../genome_wide", out_prefix+'.bedpe')

    bedpe_to_vcf = ' '.join([
        "source activate sv-channels; python ../utils/bedpe_to_vcf.py",
        "-i", os.path.join("../genome_wide", out_prefix+'.bedpe'),
        "-b", args.twobit,
        "-s", args.sample_name,
        "-o", os.path.join("../genome_wide", out_prefix +
                           '.'+args.sample_name+'.vcf')
    ])
    print(bedpe_to_vcf)
    cmd_out = subprocess.run(bedpe_to_vcf, shell=True, check=True)
    print(cmd_out)


if __name__ == '__main__':
    main()
