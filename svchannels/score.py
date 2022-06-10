import argparse
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import zarr
import sys
from tensorflow.keras.models import load_model
from .model_functions import unfold_win_id



def predict(zarr_zip, model_fn):


    X = zarr.load(zarr_zip)
    model = load_model(model_fn)
    probs = model.predict(X_test, batch_size=1000, verbose=False)
    predicted = probs.argmax(axis=1)



def main(argv=sys.argv[1:]):

    parser = argparse.ArgumentParser(description='Use model to score a set of variants')
    parser.add_argument('channels',
                        type=str,
                        default='channels.zarr.zip',
                        help="channels.zarr.zip containing SVs",
                        )
    parser.add_argument('model',
                        type=str,
                        default='best_model.keras',
                        help="TensorFlow model in HDF5 format"
                        )

    args = parser.parse_args(argv)

    windows_list = args.input.split(',')

    probs = predict(args.channels, args.model)
    predicted = probs.argmax(axis=1)




if __name__ == '__main__':
    main()
