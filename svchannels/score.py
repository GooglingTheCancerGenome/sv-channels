import argparse
import os
import zarr
import sys
from tensorflow.keras.models import load_model
import vcf


def predict(d, model_fn):
    zarr_zip = f'{d}/channels.zarr.zip'
    assert os.path.exists(zarr_zip), f'unable to access channels file in directory {d}'
    X = zarr.load(zarr_zip)
    model = load_model(model_fn)
    probs = model.predict(X, batch_size=100, verbose=True)
    return probs


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Use model to score a set of variants')
    parser.add_argument('channels_directory',
                        type=str,
                        help="directory containing channels.zarr.zip and sv_positions.bedpe",
                        )
    parser.add_argument('model',
                        type=str,
                        default='best_model.keras',
                        help="TensorFlow model in HDF5 format"
                        )
    parser.add_argument('manta_vcf',
                        type=str,
                        default='manta.vcf',
                        help="Manta VCF file"
                        )
    parser.add_argument('manta_vcf_out',
                        type=str,
                        default='manta_out.vcf',
                        help="Manta VCF file for output"
                        )

    args = parser.parse_args(argv)

    probs = predict(args.channels_directory, args.model)
    predicted = probs.argmax(axis=1)

    pos_dict = {}
    n = 0
    for i, line in enumerate(open(f'{args.channels_directory}/sv_positions.bedpe')):
        p = predicted[i]
        chrom1, pos1a, pos1b, chrom2, pos2a, pos2b, svtype = line.split('\t')
        pos_dict[chrom1 + '_' + str(int(pos1a) + 1)] = str(probs[i][0])
        n += 1
    assert n == len(predicted), "number of variants and channel-sets should match"

    reader = vcf.Reader(open(args.manta_vcf, 'r'))
    writer = vcf.Writer(open(args.manta_vcf_out, 'w'), reader)

    for record in reader:
        k = record.CHROM + '_' + str(record.POS)
        if k in pos_dict.keys():
            record.QUAL = str(pos_dict[k])
            writer.write_record(record)


if __name__ == '__main__':
    main()
