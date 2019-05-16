from pysam import VariantFile
import argparse
from time import time


def extract_chr(vcf_input_file, vcf_output_file, sv_caller):
    print('SV caller: %s' % sv_caller)

    vcf_in = VariantFile(vcf_input_file)
    vcf_out = VariantFile(vcf_output_file, 'w', header=vcf_in.header)

    n_total = 0
    n_written = 0

    for record in vcf_in:
        if record.chrom in ['1', '2', '3']: # and record.info['SVTYPE'] == 'DEL':
            vcf_out.write(record)
            n_written += 1
        n_total += 1

    print('Total:%d records' % n_total)
    print('Written:%d records' % n_written)


def main():

    context = 'trio/NA12878'

    # working director

    input_file = '/Users/lsantuari/Documents/Data/germline/' + context + '/SV/Filtered/gridss.sym.vcf'
    output_file = '/Users/lsantuari/Documents/Data/germline/' + context + '/SV/Filtered/Test_set/gridss.test_set.vcf'

    parser = argparse.ArgumentParser(description='Extract Chr1, Chr2, Chr3')
    parser.add_argument('-i', '--input', type=str, default=input_file,
                        help="Specify input file (VCF)")
    parser.add_argument('-o', '--output', type=str, default=output_file,
                        help="Specify output (VCF)")
    parser.add_argument('-c', '--caller', type=str, default='gridss',
                        help="Specify SV caller")

    args = parser.parse_args()
    #
    # for sv_caller in ['lumpy', 'gridss', 'delly', 'manta', 'last_nanosv.sorted']:
    #
    #     input_file = '/Users/lsantuari/Documents/Data/germline/' + context + '/SV/Filtered/'+sv_caller+'.sym.vcf'
    #     output_file = '/Users/lsantuari/Documents/Data/germline/' + context + '/SV/Filtered/Test_set/'+\
    #                  sv_caller+'.test_set.vcf'

    for sv_caller in ['Mills']:

        input_file = '/Users/lsantuari/Documents/External_GitHub/' + \
                     'sv_benchmark/input.na12878/lumpy-Mills2012-call-set.chr.vcf'
        output_file = '/Users/lsantuari/Documents/External_GitHub/' + \
                      'sv_benchmark/input.na12878/lumpy-Mills2011-call-set.test_set.vcf'
        t0 = time()
        extract_chr(vcf_input_file=input_file, vcf_output_file=output_file, sv_caller=sv_caller)
        print(time() - t0)


if __name__ == '__main__':
    main()
