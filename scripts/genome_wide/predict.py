import argparse
import os
import pandas as pd

from tensorflow.keras.utils import to_categorical
from tensorflow.keras.models import load_model

from model_functions import \
    evaluate_model, get_data


def predict(input_data, sample_name, svtype, model_fn, model_name, output_dir):

    os.makedirs(output_dir, exist_ok=True)

    X, y, win_ids = get_data(input_data, True, svtype)
    y_binary = to_categorical(y, num_classes=params['n_classes'])

    print('Predicting {} with model trained on {}...'.format(sample_name, model_name))

    model = load_model(model_fn)

    results = pd.DataFrame()

    intermediate_results, metrics = evaluate_model(model,
                                                   X, y_binary, win_ids,
                                                   results,
                                                   params['mapclasses'],
                                                   output_dir, svtype)

    results = results.append(intermediate_results)
    results.to_csv(os.path.join(output_dir, 'metrics.csv'), sep='\t')


def main():

    default_win = 200
    default_path = os.path.join('./cnn/win'+str(default_win), 'split_reads')
    def_windows_file = os.path.join(default_path, 'windows', 'DEL', 'windows_en.npz')

    parser = argparse.ArgumentParser(description='Use model to predict')

    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='./cnn/win200/split_reads/cv/DEL/1/model.hdf5',
                        help="TensorFlow model in HDF5 format"
                        )
    parser.add_argument('-n',
                        '--model_name',
                        type=str,
                        default='cnn',
                        help="Name of the model"
                        )
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        default=def_windows_file,
                        help="Specify training sample"
                        )
    parser.add_argument('-sn',
                        '--sample_name',
                        type=str,
                        default='test',
                        help="Specify sample name"
                        )
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='./predictions',
                        help="Output folder"
                        )

    args = parser.parse_args()

    global params
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}

    # Parameters
    params = {
        'mapclasses': mapclasses,
        'n_classes': len(mapclasses.keys())
    }

    predict(args.input, args.sample_name, args.svtype,
            args.model, args.model_name, args.output)


if __name__ == '__main__':
    main()
