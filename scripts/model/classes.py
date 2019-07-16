import numpy as np
import keras


class DataGenerator(keras.utils.Sequence):
    'Generates data for Keras'

    def __init__(self, chr_array, list_IDs, labels, batch_size=32, dim=410, n_channels=34,
                 n_classes=2, shuffle=True):
        'Initialization'
        self.dim = dim
        self.batch_size = batch_size
        self.chr_array = chr_array
        self.labels = labels
        self.list_IDs = list_IDs
        self.n_channels = n_channels
        self.n_classes = n_classes
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        'Denotes the number of batches per epoch'
        return int(np.floor(len(self.list_IDs) / self.batch_size))

    def __getitem__(self, index):
        'Generate one batch of data'
        # Generate indexes of the batch
        indexes = self.indexes[index * self.batch_size:(index + 1) * self.batch_size]

        # Find list of IDs
        list_IDs_temp = [self.list_IDs[k] for k in indexes]

        # Generate data
        X, y = self.__data_generation(list_IDs_temp)

        return X, y

    def on_epoch_end(self):
        'Updates indexes after each epoch'
        self.indexes = np.arange(len(self.list_IDs))
        if self.shuffle == True:
            np.random.shuffle(self.indexes)

    def __data_generation(self, list_IDs_temp):
        'Generates data containing batch_size samples'  # X : (n_samples, *dim, n_channels)
        # Initialization
        X = np.zeros((self.batch_size, self.dim, self.n_channels))
        y = np.zeros((self.batch_size), dtype=int)

        assert len(list_IDs_temp) > 0

        # Generate data
        for i, ID in enumerate(list_IDs_temp):

            # print(ID)
            # Store sample
            chr1, pos1, chr2, pos2 = ID.split('_')
            pos1 = int(pos1)
            pos2 = int(pos2)

            if chr1 == chr2 and chr1 == '17' and \
                    100 < pos1 < self.chr_array[chr1].shape[0] - 100 and \
                    100 < pos2 < self.chr_array[chr2].shape[0] - 100:
                # print(X[i,].shape)

                X[i,] = np.concatenate(
                    (
                        np.array(self.chr_array[chr1][pos1 - 100:pos1 + 100, :]),
                        np.zeros(shape=(10, self.chr_array[chr1].shape[1]), dtype=np.float32),
                        np.array(self.chr_array[chr2][pos2 - 100:pos2 + 100, :])
                    ),
                    axis=0)

                # Store class
                y[i] = self.labels[ID]
                assert self.labels[ID] < self.n_classes
                # print(Counter(y))
        # print(X)
        # print(y)
        assert np.array_equal(X, np.zeros((self.batch_size, self.dim, self.n_channels))) != True

        return X, keras.utils.to_categorical(y, num_classes=self.n_classes)
