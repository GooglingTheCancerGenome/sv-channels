import pysam
import argparse
import logging
from time import time
import numpy as np

from functions import load_windows, save_windows, is_left_clipped, is_right_clipped


def init_log(logfile):
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfile,
                        filemode='w',
                        level=logging.INFO)


def parse_args():

    parser = argparse.ArgumentParser(
        description='Add window specific channels')

    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-w',
                        '--win',
                        type=int,
                        default=250,
                        help="Window size")
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default='./labels/win200/DEL/split_reads/windows/windows.npz',
                        help="input file")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='./labels/win200/DEL/split_reads/windows/windows_en.npz',
                        help="output file")
    parser.add_argument('-l',
                        '--logfile',
                        default='./labels/win200/DEL/split_reads/windows/windows_en.log',
                        help='File in which to write logs.')
    parser.add_argument('-lp',
                        '--log_every_n_pos',
                        type=int,
                        default=1000,
                        help='Write in log file every log_every_n_pos positions')
    parser.add_argument('-p',
                        '--padding',
                        type=int,
                        default=10,
                        help="Length of the padding in between windows")

    return parser.parse_args()


def get_channels():

    ch = [
        # All reads (clipped or not)
        'F_AR_N', 'R_AR_N',
        # Split reads
        'F_SR_L', 'F_SR_R', 'F_SR_B', 'R_SR_L', 'R_SR_R', 'R_SR_B', 'F_SR_N', 'R_SR_N',
        # Clipped reads
        'F_CR_L', 'F_CR_R', 'R_CR_L', 'R_CR_R', 'F_CR_B', 'R_CR_B', 'F_CR_N', 'R_CR_N',
        # Discordant reads
        'DR_F', 'DR_R',
        # SV type channels
        'DUP_A', 'DUP_B', 'INV_A', 'INV_B', 'TRA_O', 'TRA_S'
    ]
    ch = {k: v for v, k in enumerate(ch)}
    return ch


def update_channel(X, ch, iter, read, win_mid_pos, is_second_win, win_len, padding):

    if is_left_clipped(read) and not is_right_clipped(read):
        clipping = 'L'
    elif not is_left_clipped(read) and is_right_clipped(read):
        clipping = 'R'
    elif is_left_clipped(read) and is_right_clipped(read):
        clipping = 'B'
    else:
        clipping = 'N'

    if read.has_tag('SA'):
        clipped_state = 'SR'
    elif is_left_clipped(read) or is_right_clipped(read):
        clipped_state = 'CR'
    else:
        clipped_state = 'AR'

    orientation = 'R' if read.is_reverse else 'F'

    start_win = win_len + padding if is_second_win else 0
    end_win = win_len*2 + padding if is_second_win else win_len

    abs_start = int(win_mid_pos - win_len / 2)
    abs_end = int(win_mid_pos + win_len / 2)

    start = max(read.reference_start, abs_start)
    end = min(read.reference_end, abs_end)

    # print('reference_start:{}, reference_end:{}'.format(s0, e0))
    rel_start = start_win + start - abs_start
    rel_end = start_win + end - abs_start

    skip = False
    if is_left_clipped(read):
        if (is_second_win and win_len + padding <= rel_start < win_len*2 + padding) or \
                (not is_second_win and 0 <= rel_start < win_len):
            rel_pos = rel_start
        else:
            skip = True
    elif is_right_clipped(read):
        if (is_second_win and win_len + padding <= rel_end < win_len*2 + padding) or \
                (not is_second_win and 0 <= rel_end < win_len):
            rel_pos = rel_end
        else:
            skip = True
    else:
        rel_pos = np.arange(max(start_win, rel_start), min(rel_end, end_win))

    if not skip:

        # print('relative reference_start:{}, relative reference_end:{}'.format(s, e))
        k = '_'.join([orientation, clipped_state, clipping])

        if k in ch.keys():
            X[iter, rel_pos, ch[k]] += 1

        if not read.is_proper_pair:
            k = '_'.join(['DR', orientation])
            if k in ch.keys():
                X[iter, rel_pos, ch[k]] += 1

        if read.is_reverse and not read.mate_is_reverse \
                and read.reference_start < read.next_reference_start:
            k = '_'.join(['DUP', 'A'])
            if k in ch.keys():
                X[iter, rel_pos, ch[k]] += 1

        if not read.is_reverse and read.mate_is_reverse \
                and read.reference_start > read.next_reference_start:
            k = '_'.join(['DUP', 'B'])
            if k in ch.keys():
                X[iter, rel_pos, ch[k]] += 1

        if read.is_reverse == read.mate_is_reverse:
            if read.reference_start < read.next_reference_start:
                k = '_'.join(['INV', 'B'])
                if k in ch.keys():
                    X[iter, rel_pos, ch[k]] += 1
            else:
                k = '_'.join(['INV', 'A'])
                if k in ch.keys():
                    X[iter, rel_pos, ch[k]] += 1

        if read.reference_name != read.next_reference_name:
            if read.is_reverse == read.mate_is_reverse:
                if read.reference_start < read.next_reference_start:
                    k = '_'.join(['TRA', 'S'])
                    if k in ch.keys():
                        X[iter, rel_pos, ch[k]] += 1
                else:
                    k = '_'.join(['TRA', 'O'])
                    if k in ch.keys():
                        X[iter, rel_pos, ch[k]] += 1

    return X


def add_channels(args):

    ibam = args.bam
    win = args.win
    ifile = args.input
    padding = args.padding
    log_every_n_pos = args.log_every_n_pos

    def get_reads(chrom, pos):
        with pysam.AlignmentFile(ibam, "rb") as aln:
            return [read for read in aln.fetch(chrom, pos - win / 2, pos + win / 2)]

    # Load the windows
    logging.info("Loading windows...")
    last_t = time()
    X, y = load_windows(ifile)
    now_t = time()
    logging.info("Windows loaded in {} seconds".format(now_t - last_t))

    # Load the channels
    ch = get_channels()
    # get starting time
    last_t = time()
    # Initialize numpy array
    X_enh = np.zeros(shape=(X.shape[:2] + (len(ch),)), dtype=np.int8)

    for i, p in enumerate(y.keys(), start=0):

        # Every n_r alignments, write log informations
        if not i % log_every_n_pos and i != 0:
            # Record the current time
            now_t = time()
            logging.info("%d positions processed (%f positions / s)" %
                         (i, log_every_n_pos / (now_t - last_t)))
            last_t = time()

        # Get genomic coordinates
        chrom1, pos1, chrom2, pos2, strand_info = p.split('_')
        pos1, pos2 = int(pos1), int(pos2)

        # Fetch reads overlapping each window
        win1_reads = get_reads(chrom1, pos1)
        win2_reads = get_reads(chrom2, pos2)

        # Which reads are in both windows?
        win1_read_names_set = set([read.query_name for read in win1_reads])
        win2_read_names_set = set([read.query_name for read in win2_reads])
        common_read_names = win1_read_names_set & win2_read_names_set

        # Only consider reads common to both windows
        win1_reads = {r for r in win1_reads if r.query_name in common_read_names and not r.is_unmapped}
        win2_reads = {r for r in win2_reads if r.query_name in common_read_names and not r.is_unmapped}

        for r in win1_reads:
            X_enh = update_channel(X_enh, ch, i, r, pos1, False, win, padding)

        for r in win2_reads:
            X_enh = update_channel(X_enh, ch, i, r, pos2, True, win, padding)

    for i in np.arange(X_enh.shape[2]):
        logging.info('win channels array:'+ \
                     'non-zero elements at index {}:{}'.format(i, np.argwhere(X_enh[i, :] != 0).shape[0]))

    X = np.concatenate((X, X_enh), axis=2)
    print(X.shape)

    for i in np.arange(X.shape[2]):
        logging.info('full channels array:' + \
                     'NaN elements at index {}:{}'.format(i, len(np.argwhere(np.isnan(X[i, :])))))

    return X, y


def main():

    # parse arguments
    args = parse_args()
    # initialize log file
    init_log(args.logfile)

    t0 = time()

    X, y = add_channels(args)

    save_windows(X, y, args.output)

    logging.info('Finished in %f seconds' % (time() - t0))


if __name__ == '__main__':
    main()
