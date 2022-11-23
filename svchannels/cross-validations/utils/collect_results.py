import argparse
import vcf
import glob
import json


def main():

    parser = argparse.ArgumentParser(description='Collect results from dictionaries created by'
                                                 '1_nested_locso_cv_per_chrom and create the new Manta VCF callset'
                                                 'with sv-channels QUALs')

    parser.add_argument('-id',
                        '--input_dir',
                        type=str,
                        default='.',
                        help="results directory with JSON pos_dict files")
    parser.add_argument('-i',
                        '--manta_vcf_in',
                        type=str,
                        default='~/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/'
                                '1KG_trios/Manta/HG00420/manta.vcf',
                        help="Manta callset in input for the test sample in the outer CV")
    parser.add_argument('-o',
                        '--manta_vcf_out',
                        type=str,
                        default='HG00420.sv-channels.vcf',
                        help="Manta callset in output for the test sample in the outer CV")
    args = parser.parse_args()

    # absolute path to search all text files inside a specific folder
    path = ''.join([r'', args.input_dir, '/**/*pos_dict.json'])
    files = glob.glob(path, recursive=True)
    print('Considering JSON files: {}'.format(files))

    # Load JSON data into a dictionary
    pos_dict = {}
    for f in files:
        with open(f, 'r') as file:
            data = json.load(file)
            pos_dict.update(data)

    reader = vcf.Reader(open(args.manta_vcf_in, 'r'))
    writer = vcf.Writer(open(args.manta_vcf_out, 'w'), reader)

    for record in reader:
        k = record.CHROM + '_' + str(record.POS)
        if k in pos_dict.keys():
            record.QUAL = str(pos_dict[k])
            writer.write_record(record)


if __name__ == '__main__':
    main()
