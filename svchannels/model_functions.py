import logging
import os
import sys
from collections import Counter
from itertools import cycle
import gzip
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import zarr
# from mcfly import modelgen, find_architecture
from sklearn.metrics import (average_precision_score, f1_score,
                             precision_recall_curve)


def unfold_win_id(win_id):

    if len(win_id.split('/')) >= 5:
        # NOTE: we are combining sv-type with strand_info!!!
        chr1, pos1, chr2, pos2, strand_info = win_id.split('/', 4)
        return chr1, pos1, chr2, pos2, strand_info
    else:
        return None


def get_data(windows_list, labels_list, svtype):

    def filter_labels(X, y, win_ids):
        keep = [i for i, v in enumerate(y) if v in [svtype, 'no' + svtype]]
        X = X[np.array(keep)]
        y = [y[i] for i in keep]
        win_ids = [win_ids[i] for i in keep]
        return X, y, win_ids

    X = []
    y = []
    win_ids = []

    for t, l in zip(windows_list, labels_list):
        logging.info('Loading data from {}...'.format(t))
        X_partial = zarr.load(t)
        X.extend(X_partial)

        with gzip.GzipFile(l, 'r') as fin:
            labels = json.loads(fin.read().decode('utf-8'))

        y.extend(labels.values())
        win_ids.extend(labels.keys())
        logging.info('Data from {} loaded'.format(t))

    X = np.stack(X, axis=0)
    logging.info(X.shape)
    logging.info(Counter(y))
    mapclasses = {svtype: 0, 'no' + svtype: 1}
    y = np.array([mapclasses[i] for i in y])
    win_ids = np.array(win_ids)
    return X, y, win_ids


def evaluate_model(model, X_test, ytest_binary, win_ids_test,
                   results, mapclasses, output_dir, svtype):

    def write_wrong_predictions(probs, predicted, y_index, win_ids_test,
                                class_labels):
        outdir = os.path.join(output_dir, 'predictions')
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, 'wrong.bedpe')
        lines = []

        for prob, p, r, w in zip(probs, predicted, y_index, win_ids_test):
            if class_labels[p] != class_labels[r]:
                sv_score = prob[0]

                if unfold_win_id(w) is not None:

                    chr1, pos1, chr2, pos2, strand_info = unfold_win_id(w)

                    # print('{0}_{1}:{2}_{3}'.format(chr1, pos1, chr2, pos2))
                    lines.append('\t'.join([
                        str(chr1),
                        str(pos1),
                        str(int(pos1) + 1),
                        str(chr2),
                        str(pos2),
                        str(int(pos2) + 1), 'PRED:' +
                        class_labels[p] + '_TRUE:' + class_labels[r],
                        str(sv_score),
                        strand_info[0],
                        strand_info[1]
                    ]) + '\n')

        with open(outfile, 'w') as f:
            # use set to make lines unique
            for ln in lines:
                f.write(ln)

    def write_correct_predictions(probs, predicted, y_index, win_ids_test,
                                  class_labels, svtype):
        outdir = os.path.join(output_dir, 'predictions')
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, 'correct.bedpe')
        print(f"[svchannels] writing correct to {outfile}", file=sys.stderr)
        lines = []
        j = 1
        for prob, p, r, w in zip(probs, predicted, y_index, win_ids_test):
            if class_labels[p] == svtype:
                sv_score = prob[0]

                if unfold_win_id(w) is not None:

                    chr1, pos1, chr2, pos2, strand_info = unfold_win_id(w)

                    lines.append('\t'.join([
                        str(chr1),
                        str(pos1),
                        str(int(pos1) + 1),
                        str(chr2),
                        str(pos2),
                        str(int(pos2) + 1),
                        'PRED_' + class_labels[p] + '_TRUE_' +
                        class_labels[r] + '_' + str(j),
                        str(sv_score),
                        strand_info[0],
                        strand_info[1]
                    ]) + '\n')
                    j += 1

        with open(outfile, 'w') as f:
            # use set to make lines unique
            for ln in lines:
                f.write(ln)

    dict_sorted = sorted(mapclasses.items(), key=lambda x: x[1])
    class_labels = [i[0] for i in dict_sorted]
    n_classes = ytest_binary.shape[1]
    probs = model.predict(X_test, batch_size=1000, verbose=False)
    # columns are predicted, rows are truth
    predicted = probs.argmax(axis=1)
    print(f"predicted shape: {predicted.shape}")
    y_index = ytest_binary.argmax(axis=1)
    write_wrong_predictions(probs, predicted, y_index, win_ids_test,
                            class_labels)
    write_correct_predictions(probs, predicted, y_index, win_ids_test,
                              class_labels, svtype)
    confusion_matrix = pd.crosstab(pd.Series(y_index), pd.Series(predicted))
    confusion_matrix.index = [class_labels[i] for i in confusion_matrix.index]
    confusion_matrix.columns = [class_labels[i]
                                for i in confusion_matrix.columns]
    confusion_matrix.reindex(columns=class_labels, fill_value=0)
    confusion_matrix.to_csv(os.path.join(
        output_dir, 'confusion_matrix.csv'), sep='\t')

    # dictionaries for each class
    precision = dict()
    recall = dict()
    f1_score_metric = dict()
    thresholds = dict()
    average_precision = dict()

    # for i in range(n_classes):
    for k, i in mapclasses.items():
        precision[k], recall[k], thresholds[k] = precision_recall_curve(
            ytest_binary[:, i], probs[:, i])
        average_precision[k] = average_precision_score(ytest_binary[:, i],
                                                       probs[:, i],
                                                       average="weighted")
        # f1_score_metric[k] = f1_score(y_index, predicted, average=None)[i]

    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(
        ytest_binary.ravel(), probs.ravel())
    average_precision["weighted"] = average_precision_score(ytest_binary,
                                                            probs,
                                                            average="weighted")
    f1_score_metric["weighted"] = f1_score(y_index,
                                           predicted,
                                           average="weighted")
    results = results.append(
        {
            "test_set_size": X_test.shape[0],
            "average_precision_score": average_precision["weighted"],
            "f1_score": f1_score_metric["weighted"]
        }, ignore_index=True)

    plot_precision_recall(mapclasses, precision, recall,
                          average_precision, output_dir)
    return results, (average_precision, precision, recall, thresholds, f1_score_metric)


def plot_precision_recall(mapclasses, precision, recall,
                          average_precision, output_dir):
    # setup plot details
    colors = cycle(
        ['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])

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

    plt.savefig(os.path.join(
        output_dir, 'precision_vs_recall.png'),
        bbox_inches='tight')
    plt.close()
