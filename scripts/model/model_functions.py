from mcfly import modelgen, find_architecture
import logging
import numpy as np
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import average_precision_score
import os, errno
import pandas as pd
import matplotlib.pyplot as plt

date = '250719'


def unfold_win_id(win_id):

    chr1, pos1, chr2, pos2 = win_id.split('_')

    return chr1, pos1, chr2, pos2


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


def create_model_with_mcfly(X, y_binary):

    models = modelgen.generate_models(X.shape,
                                      y_binary.shape[1],
                                      number_of_models=1,
                                      model_type='CNN',
                                      cnn_min_layers=2,
                                      cnn_max_layers=2,
                                      cnn_min_filters=4,
                                      cnn_max_filters=4,
                                      cnn_min_fc_nodes=6,
                                      cnn_max_fc_nodes=6,
                                      low_lr=4, high_lr=4,
                                      low_reg=1, high_reg=1,
                                      kernel_size=7)

    i = 0
    for model, params, model_types in models:
        logging.info('model ' + str(i))
        i = i + 1
        logging.info(params)
        logging.info(model.summary())

    return models


def train_model_with_mcfly(model, xtrain, ytrain, xval, yval):

    train_set_size = xtrain.shape[0]
    nr_epochs = 5

    histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(xtrain, ytrain,
                                                                                      xval, yval,
                                                                                      model, nr_epochs=nr_epochs,
                                                                                      subset_size=train_set_size,
                                                                                      verbose=True)

    best_model_index = np.argmax(val_accuracies)
    best_model, best_params, best_model_types = model[best_model_index]
    # print(best_model_index, best_model_types, best_params)

    history = best_model.fit(xtrain, ytrain,
                             epochs=nr_epochs,
                             validation_data=(xval, yval),
                             verbose=True,
                             shuffle=True)

    return history, best_model


def evaluate_model(model, X_test, ytest_binary, win_ids_test,
                   results, cv_iter, output, mapclasses, output_dir):

    def write_bed_wrong_predictions(predicted, y_index, win_ids_test, class_labels):

        #print(class_labels)

        outdir = os.path.join(output_dir, 'predictions')
        create_dir(outdir)

        outfile = os.path.join(outdir, output + '_wrong_predictions_' + str(int(cv_iter) + 1) + '.bedpe')

        lines = []

        for p, r, w in zip(predicted, y_index, win_ids_test):

            if class_labels[p] != class_labels[r]:

                chr1, pos1, chr2, pos2 = unfold_win_id(w)
                # print('{0}_{1}:{2}_{3}'.format(chr1, pos1, chr2, pos2))
                lines.append('\t'.join([str(chr1), str(pos1), str(int(pos1)+1),
                                        str(chr2), str(pos2), str(int(pos2)+1),
                                        'PRED:' + class_labels[p] + '_TRUE:' + class_labels[r]])+'\n')

        f = open(outfile, 'w')
        try:
            # use set to make lines unique
            for l in lines:
                f.write(l)
        finally:
            f.close()

    def write_bed_predictions(predicted, y_index, win_ids_test, class_labels):

        # print(class_labels)

        outdir = os.path.join(output_dir, 'predictions')
        create_dir(outdir)

        outfile = os.path.join(outdir, output + '_DEL_predicted_' + str(int(cv_iter) + 1) + '.bedpe')

        lines = []

        for p, r, w in zip(predicted, y_index, win_ids_test):

            if class_labels[p] == 'DEL':
                chr1, pos1, chr2, pos2 = unfold_win_id(w)
                # print('{0}_{1}:{2}_{3}'.format(chr1, pos1, chr2, pos2))
                lines.append('\t'.join([str(chr1), str(pos1), str(int(pos1) + 1),
                                        str(chr2), str(pos2), str(int(pos2) + 1),
                                        'PRED:' + class_labels[p] + '_TRUE:' + class_labels[r]]) + '\n')

        f = open(outfile, 'w')
        try:
            # use set to make lines unique
            for l in lines:
                f.write(l)
        finally:
            f.close()


    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    # print(dict_sorted)
    class_labels = [i[0] for i in dict_sorted]

    n_classes = ytest_binary.shape[1]
    # print(y_binarized)
    # print(n_classes)

    probs = model.predict_proba(X_test, batch_size=10000, verbose=True)

    # save model
    outdir = os.path.join(output_dir, 'models')
    create_dir(outdir)
    model.save(os.path.join(outdir, '{0}_model_{1}.hdf5'.format(output, str(int(cv_iter) + 1))))

    # columns are predicted, rows are truth
    predicted = probs.argmax(axis=1)
    # print(predicted)
    # true
    y_index = ytest_binary.argmax(axis=1)
    # print(y_index)

    # write predictions
    write_bed_wrong_predictions(predicted, y_index, win_ids_test, class_labels)
    write_bed_predictions(predicted, y_index, win_ids_test, class_labels)

    # print(y_index)
    outdir = os.path.join(output_dir, 'confusion_matrix')
    create_dir(outdir)

    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [class_labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [class_labels[i] for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=[l for l in class_labels], fill_value=0)
    confusion_matrix.to_csv(
        os.path.join(outdir, '{0}_confusion_matrix_{1}.csv'.format(output, str(int(cv_iter) + 1))), sep='\t'
    )

    # For each class
    precision = dict()
    recall = dict()
    f1_score_metric = dict()
    thresholds = dict()
    average_precision = dict()

    # for i in range(n_classes):
    for k, i in mapclasses.items():
        precision[k], recall[k], thresholds[k] = precision_recall_curve(ytest_binary[:, i],
                                                                        probs[:, i])
        average_precision[k] = average_precision_score(ytest_binary[:, i], probs[:, i], average="weighted")
        # f1_score_metric[k] = f1_score(y_index, predicted, average=None)[i]

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(ytest_binary.ravel(),
                                                                    probs.ravel())

    average_precision["weighted"] = average_precision_score(ytest_binary, probs, average="weighted")
    print('Average precision score, weighted over all classes: {0:0.2f}'
          .format(average_precision["weighted"]))

    # f1_score_metric["weighted"] = f1_score(y_index, predicted, average="weighted")

    results = results.append({
        "run": str(cv_iter + 1),
        "test_set_size": X_test.shape[0],
        "average_precision_score": average_precision["weighted"]
        # "f1_score": f1_score_metric["weighted"]
    }, ignore_index=True)

    plot_precision_recall(cv_iter, mapclasses,
                          precision, recall, average_precision, output, output_dir)

    return results, (average_precision, precision, recall, thresholds, f1_score_metric)


def plot_precision_recall(cv_iter, mapclasses, precision, recall, average_precision, output, output_dir):

    outdir = os.path.join(output_dir, 'plots')
    create_dir(outdir)

    from itertools import cycle
    # setup plot details
    colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])

    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

    lines.append(l)
    labels.append('iso-f1 curves')
    l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=2)
    lines.append(l)
    labels.append('weighted-average Precision-recall (area = {0:0.2f})'
                  ''.format(average_precision["weighted"]))

    for i, color in zip(mapclasses.keys(), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                      ''.format(i, average_precision[i]))

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Extension of Precision-Recall curve to multi-class')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))

    plt.savefig(os.path.join(outdir, '{0}_PrecRec_{1}.png'.format(output, str(int(cv_iter) + 1))),
                bbox_inches='tight')
    plt.close()
