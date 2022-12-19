import zarr
import sys
import argparse
import numpy as np


def shift_discordant(a):

    idx = np.array([23, 24, 25, 26, 31, 32, 33, 34])

    idx_a = a.expand * 2
    idx_b = a.expand * 2 + a.gap

    Z = zarr.load(a.input)

    tmp = Z[:, idx, :idx_a]
    Z[:, idx, :idx_a] = np.zeros(shape=tmp.shape)
    Z[:, idx, :idx_a - a.shift] = tmp[:, :, :idx_a - a.shift]

    tmp = Z[:, idx, idx_b:]
    Z[:, idx, idx_b:] = np.zeros(shape=tmp.shape)
    Z[:, idx, idx_b + a.shift:] = tmp[:, :, a.shift:]

    zarr.save_array(a.output, Z)


def main(args=sys.argv[1:]):

    p = argparse.ArgumentParser()
    p.add_argument(
        '--input',
        type=str,
        default="channels.zarr.zip",
        help="Input Zarr file for the channels"
    )
    p.add_argument(
        '--output',
        type=str,
        default="out.channels.zarr.zip",
        help="Output Zarr file for the channels"
    )
    p.add_argument(
        '--shift',
        type=int,
        default=10,
        help="Specify shift for the discordant channels. The left window will be shifted to the left, while" +
             "the right window will be shifted to the right")
    p.add_argument(
        '--expand',
        type=int,
        default=62,
        help="Specify expansion around breakpoints (default: %(default)s) for creating the channel. " +
             "This will automatically extend 8x the specified value for discordant reads")
    p.add_argument(
        '--gap',
        type=int,
        default=8,
        help="Specify width of gap (default: %(default)s)")

    a = p.parse_args(args)

    shift_discordant(a)


if __name__ == "__main__":
    main()
