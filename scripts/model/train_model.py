import dask.array as da
import h5py
import os
import numpy as np
import json
import gzip
from classes import DataGenerator

from keras.models import Sequential
from keras.layers import Dense, Activation, Convolution1D, Lambda, \
    Convolution2D, Flatten, \
    Reshape, LSTM, Dropout, TimeDistributed, BatchNormalization
from keras.regularizers import l2
from keras.optimizers import Adam


def get_candidate_positions():

    #channel_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/'
    channel_dir = '/Users/lsantuari/Documents/Processed/channel_maker_output'

    sampleName = 'T1'
    outFile = os.path.join(channel_dir, sampleName, 'candidate_positions.json.gz')
    with gzip.GzipFile(outFile, 'r') as fin:
        candpos = json.loads(fin.read().decode('utf-8'))

    candpos_list = []
    for k in candpos.keys():
        candpos_list.extend(candpos[k])

    data_id = ['_'.join((chr1, str(pos1), chr2, str(pos2))) for chr1, pos1, chr2, pos2 in candpos_list if
                   chr1 == chr2 and chr1 == '17']

    return data_id


def create_model(dim_length, dim_channels, class_number):

    layers = 2
    filters = [4] * layers
    fc_hidden_nodes = 6
    learning_rate = 4
    regularization_rate = 1
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
        model.add(Convolution1D(filter_number, kernel_size=kernel_size, padding='same',
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


def train():

    def load_chr_array():

        chrlist = list(map(str, range(1, 23)))
        chrlist.extend(['X'])
        chr_array = dict()

        for c in chrlist:

            chrname = 'chr'+c
            hdf5_file = os.path.join(channel_data_dir, 'T1_'+c+'.hdf5')
            f = h5py.File(hdf5_file)
            d = f[chrname]
            chr_array[c] = da.from_array(d, chunks=("auto", -1))

        return chr_array

    # Parameters
    params = {'dim': 410,
              'batch_size': 2048,
              'n_classes': 2,
              'n_channels': 34,
              'shuffle': True}

    data_id = get_candidate_positions()
    chr_array = load_chr_array()

    # channel_data_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/NA12878/'
    channel_data_dir = '/Users/lsantuari/Documents/Processed/channel_maker_output'


    # Datasets
    partition = {'train': data_id[:7000000], 'validation': data_id[7000000:]}
    labels = {k: i for k, i in zip(data_id, np.random.randint(2, size=len(data_id)))}

    # Generators
    training_generator = DataGenerator(chr_array, partition['train'], labels, **params)
    validation_generator = DataGenerator(chr_array, partition['validation'], labels, **params)

    # Design model
    model = create_model(params['dim'], params['n_channels'], params['n_classes'])

    # Train model on dataset
    model.fit_generator(generator=training_generator,
                        validation_data=validation_generator,
                        max_queue_size=2000,
                        epochs= 50,
                        workers=1,  # I don't see multi workers can have any performance benefit without multi threading
                        use_multiprocessing=False,  # HDF5Matrix cannot support multi-threads
                        )

def main():

    train()

if __name__ == '__main__':
    main()