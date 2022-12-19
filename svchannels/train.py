import argparse
import os
import logging
import zarr
import gzip
import json
import matplotlib.pyplot as plt
from pathlib import Path
import tensorflow as tf
import numpy as np
from time import time
from sklearn.utils.class_weight import compute_class_weight
from tensorflow.keras.utils import to_categorical

from svchannels.model import create_model

# specify default parameters:
# 1. number of CNN filters
# 2. number of CNN layers
# 3. size of the CNN filter (kernel)
# 4. number of nodes in the dense layer
# 5. dropout rate
# 6. initial learning rate
# 7. regularization rate
default_parameters_names = ['cnn_filters', 'cnn_layers', 'cnn_filter_size', 'fc_nodes',
                            'dropout_rate', 'learning_rate', 'regularization_rate']
default_parameters = [4, 1, 7, 4, 0.2, 1e-4, 1e-1]
default_parameters_dict = {k: v for k, v in zip(default_parameters_names, default_parameters)}


def load_windows(win_file, lab_file):
    X = zarr.load(win_file)
    with gzip.GzipFile(lab_file, 'r') as fin:
        y = json.loads(fin.read().decode('utf-8'))
    return X, y


def train(args):

    randomState = 42
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    windows = args.windows.split(',')
    labels = args.labels.split(',')
    assert len(windows) == len(labels)

    # Load channels and labels
    X = []
    y = []
    for w, l in zip(windows, labels):
        partial_X, partial_y = load_windows(w, l)
        assert partial_X.shape[0] == len(partial_y), (w, l, partial_X.shape, len(partial_y))
        X.extend(partial_X)
        y.extend(partial_y.values())
    train_X = np.stack(X, axis=0)

    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    classes = np.array(np.unique(y))
    y_lab = np.asarray(y)
    class_weights = compute_class_weight('balanced', classes, y_lab)
    train_class_weights = {i: v for i, v in enumerate(class_weights)}

    train_y = to_categorical(y, num_classes=2)

    val_X, y = load_windows(args.validation_windows, args.validation_labels)
    y = y.values()
    mapclasses = {args.svtype: 0, 'no' + args.svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    val_y = to_categorical(y, num_classes=2)

    model = create_model(train_X, 2,
                         learning_rate=default_parameters_dict['learning_rate'],
                         regularization_rate=default_parameters_dict['regularization_rate'],
                         filters=default_parameters_dict['cnn_filters'],
                         layers=default_parameters_dict['cnn_layers'],
                         kernel_size=default_parameters_dict['cnn_filter_size'],
                         fc_nodes=default_parameters_dict['fc_nodes'],
                         dropout_rate=default_parameters_dict['dropout_rate']
                         )
    # print(model.summary())

    history = model.fit(x=train_X, y=train_y,
                        epochs=args.epochs,
                        batch_size=args.batch_size,
                        shuffle=True,
                        validation_data=(val_X, val_y),
                        class_weight=train_class_weights
                        )

    # save model
    model.save(args.model)

    path = Path(args.model)
    accuracy_plot = os.path.join(path.parent.absolute(),
                                 ''.join(['model_auc_pr.png']))
    loss_plot = os.path.join(path.parent.absolute(),
                             ''.join(['model_loss.png']))

    # print(history.history)
    plt.plot(history.history['auc'])
    plt.plot(history.history['val_auc'])
    plt.title('model auc PR')
    plt.ylabel('AUC PR')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper right')
    plt.tight_layout()
    plt.savefig(accuracy_plot)
    plt.close()

    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.tight_layout()
    plt.savefig(loss_plot)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Optimize model')

    parser.add_argument('-w',
                        '--windows',
                        type=str,
                        default='sv_chan.zarr,sv_chan.zarr',
                        help="Comma separated list of training data")
    parser.add_argument('-lab',
                        '--labels',
                        type=str,
                        default='labels/labels.json.gz,labels/labels.json.gz',
                        help="Comma separated list of JSON.GZ file for labels")
    parser.add_argument('-l',
                        '--logfile',
                        default='optimize.log',
                        help='File in which to write logs.')
    parser.add_argument('-e',
                        '--epochs',
                        type=int,
                        default=50,
                        help="Number of epochs")
    parser.add_argument('-b',
                        '--batch_size',
                        type=int,
                        default=32,
                        help="Batch size")
    parser.add_argument('-n',
                        '--ncalls',
                        type=int,
                        default=50,
                        help="Number of calls of the fitness function")
    parser.add_argument('-val',
                        '--validation_windows',
                        type=str,
                        default='sv_chan.zarr',
                        help="Windows used for validation")
    parser.add_argument('-val_lab',
                        '--validation_labels',
                        type=str,
                        default='labels/labels.json.gz',
                        help="JSON.GZ file for labels")
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Type of SV")
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='best_model.keras',
                        help="Output model")
    args = parser.parse_args()

    for p in (args.model, args.logfile):
        od = Path(p)
        for parent in reversed(od.parents):
            parent.mkdir(mode=0o774, exist_ok=True)

    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(format=log_format,
                        filename=args.logfile,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    train(args)
    logging.info('Elapsed time = %f seconds' %
                 (time() - t0))


if __name__ == '__main__':
    main()
