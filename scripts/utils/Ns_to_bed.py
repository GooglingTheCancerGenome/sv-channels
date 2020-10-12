import twobitreader as twobit
import argparse


def get_regions(twobitfile, chrlist, outbed):

    genome = twobit.TwoBitFile(twobitfile)

    with (open(outbed, 'w')) as fout:

        for c in chrlist:
            seq = str(genome[c])

            s = -1
            for i, b in enumerate(seq):
                if b == 'N':
                    if s == -1:
                        s = i
                        continue
                else:
                    if s != -1:
                        e = i
                        # print('{}:{}-{} => {}'.format(c,s,e,genome[c][s:e]))
                        line = '\t'.join([c, str(s), str(e)]) + '\n'
                        fout.write(line)
                        s = -1


def main():

    parser = argparse.ArgumentParser(
        description='Extract regions with Ns')
    parser.add_argument('-b',
                        '--bed',
                        type=str,
                        default='../../data/Ns.bed',
                        help="Specify output file (BED)")
    parser.add_argument('-t',
                        '--twobit',
                        type=str,
                        default='/Users/lsantuari/Documents/Projects/GTCG/2bit/hs37d5.2bit',
                        help="Specify input file (2bit)")
    parser.add_argument('-c',
                        '--chrlist',
                        type=str,
                        default='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y',
                        help="Comma separated list of chromosomes to consider")

    args = parser.parse_args()

    get_regions(args.twobit, args.chrlist.split(','), args.bed)


if __name__ == '__main__':
    main()
