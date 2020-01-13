import os
import numpy as np
import tensorflow as tf

from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution1D, Lambda, \
    Convolution2D, Flatten, \
    Reshape, LSTM, Dropout, TimeDistributed, BatchNormalization
from keras.regularizers import l2
from keras.optimizers import Adam
from keras.callbacks import TensorBoard, EarlyStopping, ModelCheckpoint, CSVLogger
from keras.utils import to_categorical

import argparse

# set GPU options
# allow_growth allows fair share of GPU memory across processes
gpu_options = tf.GPUOptions(allow_growth=True)
sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
tf.keras.backend.set_session(sess)


def create_model(dim_length, dim_channels, class_number):

    layers = 4 # 2
    filters = [32] * layers  # 4
    fc_hidden_nodes = 8
    learning_rate = 10 ** (-4)
    regularization_rate = 10 ** (-1)
    kernel_size = 7
    drp_out1 = 0
    drp_out2 = 0

    outputdim = class_number  # number of classes

    weightinit = 'lecun_uniform'  # weight initialization

    model = Sequential()
    model.add(
        BatchNormalization(
            input_shape=(
                dim_length,
                dim_channels)))

    for filter_number in filters:
        model.add(Convolution1D(filter_number,
                                kernel_size=kernel_size,
                                padding='same',
                                kernel_regularizer=l2(regularization_rate),
                                kernel_initializer=weightinit))
        model.add(BatchNormalization())
        model.add(Activation('relu'))

    model.add(Flatten())
    model.add(Dropout(drp_out1))
    model.add(Dense(units=fc_hidden_nodes,
                    kernel_regularizer=l2(regularization_rate),
                    kernel_initializer=weightinit))  # Fully connected layer
    model.add(Activation('relu'))  # Relu activation
    model.add(Dropout(drp_out2))
    model.add(Dense(units=outputdim, kernel_initializer=weightinit))
    model.add(BatchNormalization())
    model.add(Activation("softmax"))  # Final classification layer

    model.compile(loss='categorical_crossentropy',
                  optimizer=Adam(lr=learning_rate),
                  metrics=['accuracy'])

    return model


def data(args, mapclasses):

    X_train = []
    y_train = []

    for pos_set in args.positive:
        with np.load(pos_set) as data:
            X_train.extend(data['start'])
            y_train.extend(['DEL_start'] * data['start'].shape[0])
            X_train.extend(data['end'])
            y_train.extend(['DEL_end'] * data['end'].shape[0])

    for neg_set in args.negative:
        with np.load(neg_set) as data:
            X_train.extend(data['neg'])
            y_train.extend(['noSV'] * data['neg'].shape[0])

    X_train = np.array(X_train)
    y_train = np.array(y_train)

    from collections import Counter

    sampling = 'oversample'

    cnt_lab = Counter(y_train)

    # maximum training samples per class
    max_train = 10 ** 5

    min_v = min([v for k, v in cnt_lab.items()])
    max_v = max([v for k, v in cnt_lab.items()])

    print(cnt_lab)
    print('Minimum number of labels = ' + str(min_v))
    print('Maximum number of labels = ' + str(max_v))

    data_balanced = []
    labels_balanced = []

    for l in cnt_lab.keys():
        # print(l)
        iw = np.where(y_train == l)

        if sampling == 'oversample':
            ii = np.random.choice(a=iw[0], size=min(max_v, max_train), replace=True)
        elif sampling == 'undersample':
            ii = np.random.choice(a=iw[0], size=min_v, replace=False)

        data_balanced.extend(X_train[ii])
        labels_balanced.extend(y_train[ii])

    print(Counter(labels_balanced))

    X_train = np.array(data_balanced)
    y_train = np.array(labels_balanced)

    y_train_num = [mapclasses[c] for c in y_train]
    y_train_num = np.array(y_train_num)
    y_train_binary = to_categorical(y_train_num, num_classes=len(mapclasses.keys()))
    y_train_binary.shape

    return X_train, y_train_binary


def train(args, params, X_train, y_train_binary, model):

    earlystop = EarlyStopping(monitor='val_loss',
                              min_delta=0,
                              patience=3,
                              verbose=1,
                              restore_best_weights=True)

    checkpoint = ModelCheckpoint(args.output,
                                 monitor='val_loss',
                                 mode='min',
                                 save_best_only=True,
                                 verbose=1)

    # csv_logger = CSVLogger(os.path.join('./tmp', 'training.log'))
    #
    # tbCallBack = TensorBoard(log_dir=os.path.join('./tmp/Graph'),
    #                          histogram_freq=0,
    #                          write_graph=True,
    #                          write_images=True)

    print('Fitting model...')
    # Train model on dataset
    history = model.fit(X_train, y_train_binary,
                        validation_split=params['val_split'],
                        batch_size=params['batch_size'],
                        epochs=params['epochs'],
                        shuffle=True,
                        verbose=1,
                        callbacks=[earlystop,
                                   checkpoint]
                        )


def main():

    parser = argparse.ArgumentParser(
        description='Train model',
        usage='''T0_S2_training.py [<args>]
        ''')
    parser.add_argument('-positive', nargs='+',
                            default=['positive.npz'],
                            help="Positive set file(s)")
    parser.add_argument('-negative', nargs='+',
                            default=['negative.npz'],
                            help="Positive set file(s)")
    parser.add_argument('-output', type=str,
                            default='model.hdf5',
                            help="Specify output")

    args = parser.parse_args()

    mapclasses = {'DEL_start': 0, 'DEL_end': 1, 'noSV': 2}
    X_train, y_train_binary = data(args, mapclasses)

    # Parameters
    params = {'dim': X_train.shape[1],
              'batch_size': 32,
              'epochs': 50,
              'val_split': 0.2,
              'n_classes': len(mapclasses.keys()),
              'n_channels': X_train.shape[2],
              'shuffle': True}

    model = create_model(params['dim'], params['n_channels'], params['n_classes'])

    train(args, params, X_train, y_train_binary, model)


if __name__ == '__main__':
    main()
