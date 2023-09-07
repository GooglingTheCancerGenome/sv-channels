import argparse
import os
import sys
import logging
import zarr
import gzip
import json
import matplotlib.pyplot as plt
from pathlib import Path
import tensorflow as tf
import numpy as np
from time import time
from sklearn.metrics import average_precision_score, accuracy_score, balanced_accuracy_score
from skopt import gp_minimize
from skopt.space import Real, Integer
from skopt.utils import use_named_args
from sklearn.model_selection import KFold
import numpy as np
from sklearn.utils.class_weight import compute_class_weight
from tensorflow.keras.utils import to_categorical
from sklearn.model_selection import train_test_split

from netcal.scaling import TemperatureScaling
from netcal.presentation import ReliabilityDiagram
from netcal.metrics import ECE

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

# Search space for the hyperparameters

search_space_dict = {
    'cnn_filters': Integer(low=1, high=8, name='cnn_filters'),
    'cnn_layers': Integer(low=1, high=2, name='cnn_layers'),
    'cnn_filter_size': Integer(low=1, high=8, name='cnn_filter_size'),
    'fc_nodes': Integer(low=3, high=6, name='fc_nodes'),
    'dropout_rate': Real(low=0, high=0.9, name='dropout_rate'),
    'learning_rate': Real(low=1e-4, high=1e-1, prior='log-uniform', name='learning_rate'),
    'regularization_rate': Real(low=1e-4, high=1e-1, prior='log-uniform', name='regularization_rate')
}

dimensions = [search_space_dict[k] for k in search_space_dict.keys()]

best_performance = 0.0
best_hyperparameters = {}
verbose_level = 0


class TemperatureScalingLayer(tf.keras.layers.Layer, temperature):

  def __init__(self, num_outputs, temperature):
    super(TemperatureScalingLayer, self).__init__()
    self.num_outputs = num_outputs
    self.temperature = temperature

  def build(self, input_shape):

    self.kernel = self.add_weight("kernel",
                                  shape=[int(input_shape[-1]),
                                         self.num_outputs])

  def call(self, inputs):

    return self.temperature.transform(input)


def get_labels(y):
    '''

    :param y: labels
    :return: y_lab: labels with 0 (positive class) and 1 (negative class)
             y_cat: class vector (integers) in binary class matrix form
             class_weights: weights for the unbalanced classes
    '''
    svtype = 'DEL'
    mapclasses = {svtype: 0, 'no' + svtype: 1}
    y_mapped = np.array([mapclasses[i] for i in y])
    classes = np.array(np.unique(y_mapped))
    y_lab = np.asarray(y_mapped)
    class_weights = compute_class_weight(class_weight='balanced', classes=classes, y=y_lab)
    class_weights = {i: v for i, v in enumerate(class_weights)}
    y_cat = to_categorical(y=y_lab, num_classes=2)

    return y_lab, y_cat, class_weights


def load_windows(win_file, lab_file):
    X = zarr.load(win_file)
    with gzip.GzipFile(lab_file, 'r') as fin:
        y = json.loads(fin.read().decode('utf-8'))
    return X, y


@use_named_args(dimensions=dimensions)
def fitness(cnn_filters, cnn_layers, cnn_filter_size, fc_nodes,
            dropout_rate, learning_rate, regularization_rate):

    # Save the model if it improves on the best found performance which is stored by the global variable best_auc
    global best_performance
    global best_hyperparameters

    # Print the current combination of hyperparameters
    current_hyperparameters = {
        'cnn_filters': cnn_filters,
        'cnn_layers': cnn_layers,
        'cnn_filter_size': cnn_filter_size,
        'fc_nodes': fc_nodes,
        'dropout_rate': dropout_rate,
        'learning_rate': learning_rate,
        'regularization_rate': regularization_rate
    }
    logging.info(current_hyperparameters)

    logging.info('Inner CV')

    auc_precision_recall = []
    accuracy = []
    balanced_accuracy = []
    validation_auc_pr = []

    kf = KFold(n_splits=kfolds)

    for i, (train_index, test_index) in enumerate(kf.split(X_train)):

        logging.info(f"Fold {i}:")

        inner_X_train = X_train[train_index]
        inner_y_train = y_train[train_index]
        inner_X_test = X_train[test_index]
        inner_y_test = y_train[test_index]

        logging.info(f"inner_X_train shape:{inner_X_train.shape}")
        logging.info(f"inner_y_train shape:{inner_y_train.shape}")
        logging.info(f"inner_X_test shape:{inner_X_test.shape}")
        logging.info(f"inner_y_test shape:{inner_y_test.shape}")

        _, inner_y_train, class_weights_train = get_labels(inner_y_train)
        inner_y_test, _, class_weights_test = get_labels(inner_y_test)

        model = create_model(inner_X_train, 2,
                             learning_rate=learning_rate,
                             regularization_rate=regularization_rate,
                             filters=cnn_filters,
                             layers=cnn_layers,
                             kernel_size=cnn_filter_size,
                             fc_nodes=fc_nodes,
                             dropout_rate=dropout_rate)
        # print(model.summary())

        history = model.fit(x=inner_X_train, y=inner_y_train,
                            epochs=n_epochs,
                            batch_size=batch_size,
                            shuffle=True,
                            validation_split=validation_split,
                            class_weight=class_weights_train,
                            # callbacks=[callback],
                            verbose=verbose_level)

        # Get the predicted probability of testing data
        y_probs = model.predict(inner_X_test)
        y_score = y_probs[:, 0]

        # Data to plot precision - recall curve
        # precision, recall, thresholds = precision_recall_curve(y_test, y_score)
        # Use AUC function to calculate the area under the curve of precision recall curve
        # auc_precision_recall.append(auc(recall, precision))
        # calculate average precision score
        avg_prec = average_precision_score(inner_y_test, y_score)
        auc_precision_recall.append(avg_prec)

        y_predicted = y_probs.argmax(axis=1)
        acc_score = accuracy_score(inner_y_test, y_predicted)
        accuracy.append(acc_score)

        bal_acc_score = balanced_accuracy_score(inner_y_test, y_predicted)
        balanced_accuracy.append(bal_acc_score)

        val_auc = history.history['val_auc'][-1]
        validation_auc_pr.append(val_auc)

        # Delete the model that just finishing training from memory
        del model

        # Clear the Keras session, otherwise it will keep adding new
        # models to the same TensorFlow graph each time we create
        # a model with a different set of hyperparameters.
        tf.keras.backend.clear_session()

    mean_auc = np.mean(auc_precision_recall)
    mean_accuracy = np.mean(accuracy)
    mean_balanced_accuracy = np.mean(balanced_accuracy)
    mean_validation_auc_pr = np.mean(validation_auc_pr)

    logging.info('mean_auc_precision_recall = {}'.format(mean_auc))
    logging.info('mean_accuracy = {}'.format(mean_accuracy))
    logging.info('mean_balanced_accuracy = {}'.format(mean_balanced_accuracy))
    logging.info('mean_validation_auc_pr = {}'.format(mean_validation_auc_pr))

    mean_performance = mean_balanced_accuracy

    # If the AUC of the saved model is greater than the current best performance
    if mean_performance > best_performance:
        # Store the best hyperparameters
        best_hyperparameters = {
            'cnn_filters': cnn_filters,
            'cnn_layers': cnn_layers,
            'cnn_filter_size': cnn_filter_size,
            'fc_nodes': fc_nodes,
            'dropout_rate': dropout_rate,
            'learning_rate': learning_rate,
            'regularization_rate': regularization_rate,
            'mean_auc_precision_recall': mean_auc
        }
        best_hyperparameters = {k: [v] for k, v in best_hyperparameters.items()}

        prefix = '_'.join(['fold', str(i), 'cnn-filt', str(cnn_filters),
                           'cnn-lay', str(cnn_layers), 'cnn-filt-size', str(cnn_filter_size),
                           'fc-nodes', str(fc_nodes), 'dropout', str(dropout_rate),
                           'lr', str(learning_rate), 'rr', str(regularization_rate),
                           'pr-auc', str(best_performance)])
        accuracy_plot = os.path.join(output,
                                     ''.join([prefix, '_auc_pr.png']))
        loss_plot = os.path.join(output,
                                 ''.join([prefix, '_loss.png']))

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

        # Update the current greatest AUC score
        best_performance = mean_performance

    # Scikit-optimize does minimization so it tries to
    # find a set the combination of hyperparameters with the lowest fitness value.
    # We want to maximize the model performance so we negate this number.
    return -mean_performance


def train(args):

    randomState = 42
    np.random.seed(randomState)
    tf.random.set_seed(randomState)

    windows = args.windows.split(',')
    labels = args.labels.split(',')
    assert len(windows) == len(labels)

    global X_train, y_train, n_epochs, batch_size, validation_split, output, kfolds

    n_epochs = args.epochs
    batch_size = args.batch_size
    validation_split = args.validation_split
    output = args.output
    kfolds = args.kfolds

    # Load channels and labels
    X = []
    y_train = []
    for w, l in zip(windows, labels):
        partial_X, partial_y = load_windows(w, l)
        assert partial_X.shape[0] == len(partial_y), (w, l, partial_X.shape, len(partial_y))
        X.extend(partial_X)
        y_train.extend(partial_y.values())
    X_train = np.stack(X, axis=0)
    y_train = np.array(y_train)

    search_result = gp_minimize(func=fitness,
                                dimensions=dimensions,
                                acq_func='EI',
                                n_calls=args.ncalls,
                                x0=default_parameters,
                                random_state=42,
                                n_jobs=-1,
                                verbose=verbose_level)

    model = create_model(X_train, 2,
                         learning_rate=best_hyperparameters['learning_rate'][0],
                         regularization_rate=best_hyperparameters['regularization_rate'][0],
                         filters=best_hyperparameters['cnn_filters'][0],
                         layers=best_hyperparameters['cnn_layers'][0],
                         kernel_size=best_hyperparameters['cnn_filter_size'][0],
                         fc_nodes=best_hyperparameters['fc_nodes'][0],
                         dropout_rate=best_hyperparameters['dropout_rate'][0]
                         )
    # print(model.summary())

    _, y_train, class_weights_train = get_labels(y_train)

    X_train_fit, X_val, y_train_fit, y_val = train_test_split(X_train, y_train,
                                                              test_size=validation_split,
                                                              random_state=42)

    history = model.fit(x=X_train_fit, y=y_train_fit,
                        epochs=n_epochs,
                        batch_size=batch_size,
                        shuffle=True,
                        validation_data=(X_val, y_val),
                        class_weight=class_weights_train,
                        # callbacks=[callback],
                        verbose=verbose_level)

    prefix = '_'.join(['cnn-filt', str(best_hyperparameters['cnn_filters'][0]),
                       'cnn-lay', str(best_hyperparameters['cnn_layers'][0]),
                       'cnn-filt-size', str(best_hyperparameters['cnn_filter_size'][0]),
                       'fc-nodes', str(best_hyperparameters['fc_nodes'][0]),
                       'dropout', str(best_hyperparameters['dropout_rate'][0]),
                       'lr', str(best_hyperparameters['learning_rate'][0]),
                       'rr', str(best_hyperparameters['regularization_rate'][0]),
                       'pr-auc', str(best_performance)])

    # Temperature scaling
    confidences_val = model.predict(X_val, batch_size=100, verbose=True)
    ground_truth = y_val

    temperature = TemperatureScaling()
    temperature.fit(confidences_val, ground_truth)

    confidences_train = model.predict(X_train, batch_size=100, verbose=True)
    calibrated = temperature.transform(confidences_train)

    n_bins = 10

    ece = ECE(n_bins)
    uncalibrated_score = ece.measure(confidences_train, ground_truth)
    calibrated_score = ece.measure(calibrated, ground_truth)
    logging.info('Temperature scaling: ECE uncalibrated_score = {}'.format(uncalibrated_score))
    logging.info('Temperature scaling: ECE calibrated_score = {}'.format(calibrated_score))

    diagram = ReliabilityDiagram(n_bins)

    uncalibrated_plot = os.path.join(args.output,
                                     ''.join([prefix, '_uncalibrated.png']))
    calibrated_plot = os.path.join(args.output,
                                   ''.join([prefix, '_calibrated.png']))

    diagram.plot(confidences_train, ground_truth)  # visualize miscalibration of uncalibrated
    plt.tight_layout()
    plt.savefig(uncalibrated_plot)
    plt.close()

    diagram.plot(calibrated, ground_truth)  # visualize miscalibration of calibrated
    plt.tight_layout()
    plt.savefig(calibrated_plot)
    plt.close()

    accuracy_plot = os.path.join(args.output,
                                 ''.join([prefix, '_auc_pr.png']))
    loss_plot = os.path.join(args.output,
                             ''.join([prefix, '_loss.png']))

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

    # save model
    model.save(args.model)


def main(argv=sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Optimize model')

    parser.add_argument('windows',
                        type=str,
                        default='sv_chan.zarr,sv_chan.zarr',
                        help="Comma separated list of training data")
    parser.add_argument('labels',
                        type=str,
                        default='labels/labels.json.gz,labels/labels.json.gz',
                        help="Comma separated list of JSON.GZ file for labels")
    parser.add_argument('--validation_split',
                        type=float,
                        default=0.3,
                        help="Split used for validation")
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        default='.',
                        help="Output folder")
    parser.add_argument('-l',
                        '--logfile',
                        default='optimize.log',
                        help='File in which to write logs.')
    parser.add_argument('-k',
                        '--kfolds',
                        type=int,
                        default=2,
                        help="Number of folds for cross-validation")
    parser.add_argument('-e',
                        '--epochs',
                        type=int,
                        default=10,
                        help="Number of epochs")
    parser.add_argument('-n',
                        '--ncalls',
                        type=int,
                        default=12,
                        help="Number of calls of the fitness function")
    parser.add_argument('-b',
                        '--batch_size',
                        type=int,
                        default=32,
                        help="Batch size")
    parser.add_argument('-s',
                        'python3 -m pip install netcal',
                        type=str,
                        default='DEL',
                        help="Type of SV")
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='best_model.keras',
                        help="Output model")

    args = parser.parse_args(argv)

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
