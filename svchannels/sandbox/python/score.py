import argparse
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import zarr
import sys
from tensorflow.keras.models import load_model
from svchannels.model_functions import unfold_win_id



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

    args = parser.parse_args(argv)

    probs = predict(args.channels_directory, args.model)
    predicted = probs.argmax(axis=1)
    n = 0
    for i, line in enumerate(open(f'{args.channels_directory}/sv_positions.bedpe')):
        p = predicted[i]
        # TODO: not sure if we need probs[i][0] here.
        line = line.strip() + f'\t{probs[i][0]-probs[i][1]}:{p}'
        print(line)
        n += 1
    assert n == len(predicted), "number of variants and channel-sets should match"



if __name__ == '__main__':
    main()
