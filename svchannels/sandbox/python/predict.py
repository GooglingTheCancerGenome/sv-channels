import argparse
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import zarr
import sys
from tensorflow.keras.models import load_model
from model_functions import unfold_win_id


def evaluate(model, X_test, win_ids_test, output_dir):

    def write_predictions(probs, predicted, win_ids_test):

        outdir = os.path.join(output_dir, 'predictions')
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, 'correct.bedpe')
        print(f"[svchannels] writing correct to {outfile}", file=sys.stderr)
        lines = []
        j = 1
        for prob, p, w in zip(probs, predicted, win_ids_test):

            sv_score = prob[0]

            if unfold_win_id(w) is not None:

                chr1, pos1, chr2, pos2, strand_info = unfold_win_id(w)

                lines.append('\t'.join([
                    str(chr1),
                    str(pos1),
                    str(int(pos1) + 1),
                    str(chr2),
                    str(pos2),
                    str(int(pos2) + 1),
                    'PRED_' + class_labels[p] + '_' + str(j),
                    str(sv_score),
                    strand_info[0],
                    strand_info[1]
                ]) + '\n')
                j += 1

        with open(outfile, 'w') as f:
            # use set to make lines unique
            for ln in lines:
                f.write(ln)

    dict_sorted = sorted(params["mapclasses"].items(), key=lambda x: x[1])
    class_labels = [i[0] for i in dict_sorted]

    probs = model.predict(X_test, batch_size=1000, verbose=False)
    # columns are predicted, rows are truth
    predicted = probs.argmax(axis=1)

    write_predictions(probs, predicted, win_ids_test)


def predict(input_data, input_bedpe, sample_name, model_fn, model_name, output_dir):

    os.makedirs(output_dir, exist_ok=True)

    # Load input data for test sample
    X = []
    for t in input_data:
        print('Loading data from {}...'.format(t))
        X_partial = zarr.load(t)
        X.extend(X_partial)
    X = np.stack(X, axis=0)

    # Load Manta calls for test sample
    win_ids = []
    with (open(input_bedpe, 'r')) as bed:

        for line in bed:
            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])
            svtype = columns[-1]
            win_ids.append('_'.join((chrom1, pos1_start, pos1_end,
                                    chrom2, pos2_start, pos2_end,
                                    "**")))
    win_ids = np.array(win_ids)

    print('Predicting {} with model trained on {}...'.format(
        sample_name, model_name))

    model = load_model(model_fn)
    evaluate(model, X, win_ids,output_dir)


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
    parser.add_argument('-b',
                        '--bedpe',
                        type=str,
                        default='manta.bedpe',
                        help="Specify Manta calls of the test sample",
                        required=True
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
    parser.add_argument('--exclude',
                        type=str,
                        help='BED file or regions to exclude')
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

    predict(windows_list, args.bedpe, args.sample_name,
            args.model, args.model_name, args.output)

    out_prefix = os.path.join(args.output, "sv-channels")

    merge_sv_calls = ' '.join([
        "Rscript", os.path.join(args.sv_channels, "scripts", "R", "merge_sv_calls.R"),
        "-i", args.output,
        "-f", args.exclude,
        "-m split_reads",
        "-o", os.path.join(args.output, "sv-channels")
    ])

    cmd_out = subprocess.run(merge_sv_calls, shell=True, check=True,
            stderr=sys.stderr)
    print(cmd_out)

    #assert os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py")
    #assert os.path.join(args.sv_channels, "scripts/genome_wide", args.output + '.bedpe')

    bedpe_to_vcf = ' '.join([
        os.path.join(args.sv_channels, "scripts/utils/bedpe_to_vcf.py"),
        "-i", os.path.join(args.output, "sv-channels.DEL.bedpe"),
        "-b", args.twobit,
        "-s", args.sample_name,
        "-o", os.path.join(args.output, "sv-channels.DEL.vcf")
    ])
    print(bedpe_to_vcf)
    cmd_out = subprocess.run(bedpe_to_vcf, shell=True, check=True, stderr=sys.stderr)
    print(cmd_out)


if __name__ == '__main__':
    main()
