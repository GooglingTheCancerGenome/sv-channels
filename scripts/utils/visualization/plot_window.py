import argparse
import os
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors
from sklearn.preprocessing import minmax_scale


def get_data(windows_list):

    X = []
    y = []
    win_ids = []

    for t in windows_list:

        print('Loading data from {}...'.format(t))

        npzfile = np.load(t, allow_pickle=True)

        X.extend(npzfile['data'])
        labels = npzfile['labels']
        labels = labels.item()

        y.extend(labels.values())
        win_ids.extend(labels.keys())

        print('Data from {} loaded'.format(t))

    X = np.stack(X, axis=0)

    print('X shape:{}'.format(X.shape))
    print('y:{}'.format(Counter(y)))

    return X, y, win_ids


def plot_window(X, y, w,  idx, outdir):

    W_i = minmax_scale(X[idx, :, :], feature_range=(0, 1), axis=0, copy=True)

    df = pd.DataFrame(W_i, columns=np.arange(W_i.shape[1]))
    # print(df)
    ax = df.plot(subplots=True,
                 figsize=(20, 5),
                 kind='line',
                 legend=False,
                 color='black')

    for i, a in enumerate(ax):
        if i != W_i.shape[1] - 1:
            a.spines['bottom'].set_color('white')
        if i != 0:
            a.spines['top'].set_color('white')
        a.set_yticks([])
        a.set_ylim(0, 1)

    plt.tight_layout()

    img_type = 'jpg'
    figname = '.'.join([w[idx], y[idx], img_type])
    plt.savefig(os.path.join(outdir, figname), dpi=300, format=img_type)


def main():

    default_win = 200
    default_path = os.path.join(
        '../genome_wide/cnn/win'+str(default_win), 'split_reads')
    def_windows_file = os.path.join(
        default_path, 'windows', 'DEL', 'windows_en.npz')

    parser = argparse.ArgumentParser(description='Use model to predict')

    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default=def_windows_file,
                        help="Specify list of windows"
                        )
    parser.add_argument('-j',
                        '--index',
                        type=int,
                        default=1,
                        help="Specify index of window to plot"
                        )
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    parser.add_argument('-c',
                        '--class_label',
                        type=str,
                        default='DEL',
                        help="Specify class")
    parser.add_argument('-o',
                        '--output_dir',
                        type=str,
                        default='plots',
                        help="Specify output directory")

    args = parser.parse_args()

    # grep 'DEL' test.bedpe | awk '{print $1"\t"$2"\t"$3"\n"$4"\t"$5"\t"$6}' > test.igv.DEL.bed

    os.makedirs(os.path.join(args.output_dir, args.svtype), exist_ok=True)

    X, y, win_ids = get_data([args.input])

    sv_i = [i for i in np.arange(len(y)) if y[i] == args.class_label]

    for i in sv_i:
        if win_ids[i] == '12_1053781_12_1054233_+-':
            plot_window(X, y, win_ids, i, args.output_dir)


if __name__ == "__main__":
    main()
