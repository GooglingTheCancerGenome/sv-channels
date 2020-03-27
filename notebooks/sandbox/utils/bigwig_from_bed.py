import pyBigWig
import argparse


def extract_regions(bwfile, bedfile, chsizefile, bwoutfile):
    bw = pyBigWig.open(bwfile)

    bed_el = []
    intervals = []

    for line in open(bedfile):
        cols = line.strip().split()
        intervals.append(bw.intervals(cols[0], int(cols[1]), int(cols[2])))
        bed_el.append((cols[0], int(cols[1]), int(cols[2])))

    bwout = pyBigWig.open(bwoutfile, 'w')

    # add header
    header = []
    for line in open(chsizefile):
        l = line.strip().split()
        header.append((l[0], int(l[1])))
    bwout.addHeader(header)

    for i, l in enumerate(intervals):
        chroms = [bed_el[i][0]] * len(l)
        starts = [t[0] for t in l]
        ends = [t[1] for t in l]
        values = [t[2] for t in l]

        bwout.addEntries(chroms, starts, ends, values)

    bw.close()
    bwout.close()


def main():

    parser = argparse.ArgumentParser(description='Extract a set of BED regions from a BigWig file')
    parser.add_argument('-bw', '--bigwig', type=str,
                        default='',
                        help="input BigWig file")
    parser.add_argument('-b', '--bed', type=str,
                        default='',
                        help="input BED file")
    parser.add_argument('-c', '--chromsizes', type=str,
                        default='',
                        help="input TSV file with chromosome size")
    parser.add_argument('-bwout', '--bigwigout', type=str,
                        default='',
                        help="output BigWig file")
    args = parser.parse_args()

    extract_regions(bwfile=args.bigwig, bedfile=args.bed,
                    chsizefile=args.chromsizes, bwoutfile=args.bigwigout)


if __name__ == '__main__':
    main()
