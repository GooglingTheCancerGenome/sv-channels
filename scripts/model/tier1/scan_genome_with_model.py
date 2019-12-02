import bcolz
import os
import argparse
import logging
import numpy as np
import json, gzip, errno
from time import time
from keras.models import load_model


def create_dir(directory):
    '''
    Create a directory if it does not exist. Raises an exception if the directory exists.
    :param directory: directory to create
    :return: None
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def load_bcolz_array(channel_data_dir, sampleName, chrom):

    carray_file = os.path.join(channel_data_dir, sampleName, 'chr_array', sampleName + '_' + chrom + '_carray')
    # logging.info('Loading file: {}'.format(carray_file))
    assert os.path.exists(carray_file), carray_file + ' not found'
    chr_array = bcolz.open(rootdir=carray_file)
    # logging.info('Array shape: {}'.format(chr_array[c].shape))

    return chr_array


def run_tier1(sampleName, channeldir, chrom, win, model_fn, outFile):

    model = load_model(model_fn)

    bc_array = load_bcolz_array(channeldir, sampleName, chrom)

    probs_list = []
    step = 10 ** 6
    batch_size = 10 ** 5

    chr_len = bc_array.shape[0] // win * win

    dim1 = (step // win) * (chr_len // step) + (chr_len % step) // win
    dim2 = int(model.outputs[0].shape.dims[1])

    res_array = np.empty(shape=(dim1, dim2))

    for i in np.arange(0, chr_len // step + 1):

        vstart = i * step
        vend = min((i + 1) * step, chr_len)
        d0 = step // win

        print('Scanning chr {} from {} to {} by {}bp'.format(chrom, vstart, vend, win))

        # B = []
        # idx = []
        #
        # for j in np.arange(200):
        #
        #     vstart_in = vstart + j
        #     vend_in = vend + j
        #
        #     if vend_in > bc_array.shape[0]:
        #         break
        #
        #     B.extend(
        #         np.split(
        #             bc_array[vstart_in:vend_in, :], step // win
        #         )
        #     )
        #     idx.extend(list(np.arange(vstart_in, vend_in)))

        # print(B.shape)

        split = (chr_len % step)//win if i == (chr_len//step) else step//win

        B = np.array(
            np.split(
                bc_array[vstart:vend, :], split
                )
        )

        # print(B.shape)
        probs = model.predict_proba(B, batch_size=batch_size, verbose=False)
        # print(probs)
        res_array[i*d0:(i+1)*d0] = probs

    # Write
    np.savez(file=outFile, probs=res_array)


def main():

    parser = argparse.ArgumentParser(description='Run Tier1')

    parser.add_argument('-l', '--logfile', type=str, default='tier1.log',
                        help="Specify log file")
    parser.add_argument('-s', '--sample', type=str, default='NA24385',
                        help="Specify sample")
    parser.add_argument('-m', '--modelpath', type=str, default='model.hdf5',
                        help="Specify model")
    parser.add_argument('-chr', '--chromosome', type=str, default='1',
                        help="Specify chromosome")
    parser.add_argument('-w', '--window', type=int, default=200,
                        help="Specify window size")
    parser.add_argument('-c', '--channeldir', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify channel data dir")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")

    args = parser.parse_args()

    # Log file
    cmd_name = 'tier1_' + str(args.window)
    output_dir = os.path.join(args.outputpath, args.sample, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, 'predictions_'+args.chromosome)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    run_tier1(args.sample,
              args.channeldir,
              args.chromosome,
              args.window,
              args.modelpath,
              output_file
              )

    logging.info('Elapsed time making labels = %f' % (time() - t0))


if __name__ == '__main__':
    main()


