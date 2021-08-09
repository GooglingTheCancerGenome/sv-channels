from __future__ import print_function

import tensorflow
print("Tensorflow version:", tensorflow.__version__)
import keras
print("Keras version:", keras.__version__)
import numpy
print("Numpy version:", numpy.__version__)

import os
import numpy as np
from keras.models import load_model
from keras.models import model_from_json

import deeplift
from deeplift.util import compile_func
from deeplift.layers import NonlinearMxtsMode
import deeplift.conversion.kerasapi_conversion as kc

from collections import OrderedDict


def deeplift(args):

    def get_entropy(b):
        return np.apply_along_axis(lambda p: -(np.sum(p * np.log(p))), 1, b)

    model = load_model(args.model)
    probs = model.predict_proba(X, batch_size=1000, verbose=True)
    probs_entropy = get_entropy(probs)

    # serialize model
    keras_model_weights = os.path.join(args.dir, "model.h5")
    keras_model_json = os.path.join(args.dir, "model.json")

    model_json = model.to_json()
    with open(keras_model_json, "w") as json_file:
        json_file.write(model_json)
    # serialize weights to HDF5
    model.save_weights(keras_model_weights)
    print("Saved model to disk")

    keras_model = model_from_json(open(keras_model_json).read())
    keras_model.load_weights(keras_model_weights)

    method_to_model = OrderedDict()
    for method_name, nonlinear_mxts_mode in [
        # The genomics default = rescale on conv layers, revealcancel on fully-connected
        ('rescale_conv_revealcancel_fc', NonlinearMxtsMode.DeepLIFT_GenomicsDefault),
        ('rescale_all_layers', NonlinearMxtsMode.Rescale),
        ('revealcancel_all_layers', NonlinearMxtsMode.RevealCancel),
        ('grad_times_inp', NonlinearMxtsMode.Gradient),
        ('guided_backprop', NonlinearMxtsMode.GuidedBackprop)]:
        method_to_model[method_name] = kc.convert_model_from_saved_files(
            h5_file=keras_model_weights,
            json_file=keras_model_json,
            nonlinear_mxts_mode=nonlinear_mxts_mode)

    # load windows
    with np.load(args.windows) as npzfile:
        input_data = npzfile['data']
        labels = npzfile['labels']
        labels = labels.item()
        input_data_labels = labels.values()

    # make sure predictions are the same as the original model
    model_to_test = method_to_model['rescale_conv_revealcancel_fc']
    deeplift_prediction_func = compile_func([model_to_test.get_layers()[0].get_activation_vars()],
                                            model_to_test.get_layers()[-1].get_activation_vars())
    original_model_predictions = keras_model.predict(input_data, batch_size=200)
    converted_model_predictions = deeplift.util.run_function_in_batches(
        input_data_list=[input_data],
        func=deeplift_prediction_func,
        batch_size=200,
        progress_update=None)
    print("maximum difference in predictions:",
          np.max(np.array(converted_model_predictions) - np.array(original_model_predictions)))
    assert np.max(np.array(converted_model_predictions) - np.array(original_model_predictions)) < 10 ** -5
    predictions = converted_model_predictions

    print("Compiling scoring functions")
    method_to_scoring_func = OrderedDict()
    for method, model in method_to_model.items():
        print("Compiling scoring function for: " + method)
        method_to_scoring_func[method] = model.get_target_contribs_func(find_scores_layer_idx=0,
                                                                        target_layer_idx=-3)

    # To get a function that just gives the gradients, we use the multipliers of the Gradient model
    gradient_func = method_to_model['grad_times_inp'].get_target_multipliers_func(find_scores_layer_idx=0,
                                                                                  target_layer_idx=-3)
    print("Compiling integrated gradients scoring functions")
    integrated_gradients10_func = deeplift.util.get_integrated_gradients_function(
        gradient_computation_function=gradient_func,
        num_intervals=10)
    method_to_scoring_func['integrated_gradients10'] = integrated_gradients10_func

    # Use mean as reference
    bg = np.mean(input_data, axis=0)

    method_to_task_to_scores = OrderedDict()
    for method_name, score_func in method_to_scoring_func.items():
        print("on method", method_name)
        method_to_task_to_scores[method_name] = OrderedDict()
        for task_idx in [0, 1, 2]:
            scores = np.array(score_func(
                task_idx=task_idx,
                input_data_list=[X],
                input_references_list=[bg],
                batch_size=32,
                progress_update=None))
            # print(scores.shape)
            assert scores.shape == X.shape
            # scores = np.sum(scores, axis=2)
            method_to_task_to_scores[method_name][task_idx] = scores

    scores_file = os.path.join(args.dir, 'method_to_task_to_scores.npy')
    np.save(scores_file, method_to_task_to_scores)

    # Load scores
    import numpy as np

    method_to_task_to_scores = np.load(scores_file, allow_pickle=True).item()
    # method_to_task_to_scores_loaded
    print(method_to_task_to_scores.keys())
    for k in method_to_task_to_scores.keys():
        for i in [0, 1, 2]:
            print(method_to_task_to_scores[k][i].shape)


def main():
    parser = argparse.ArgumentParser(description='Apply DeepLift')
    parser.add_argument('-w',
                        '--windows',
                        type=str,
                        default='',
                        help="Specify windows path")
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        default='',
                        help="Specify model path")
    parser.add_argument('-d',
                        '--dir',
                        type=str,
                        default='',
                        help="Specify working directory")

    args = parser.parse_args()
    deeplift(args)


if __name__ == '__main__':
    main()
