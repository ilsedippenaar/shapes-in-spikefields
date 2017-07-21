import numpy as np
import tensorflow as tf
from scipy.io import loadmat
from config import config


class LfpLSTM:
    def __init__(self):
        data_dir, cache_dir, plot_save_dir = config()

        start_time = (24 * 60 + 10) * 1000
        self.lfp, times, labels = self.load_data(data_dir, start_time)

        self.times = times
        label_dict = {label: i for i, label in enumerate(set(labels))}
        self.labels = np.array([label_dict[label] for label in labels]) # assign ids to each string label

        # Hyper-parameters
        self.num_epochs = 30
        self.learning_rate = 1e-2
        self.momentum = 0.5

        # input sizes
        self.num_features = self.lfp.shape[1]
        self.num_time_points = 1000  # = order of process, number of cells, or "max_time" parameter
        self.batch_size = 1  # this is a formality essentially

        # output sizes
        num_labels = len(label_dict)
        self.num_classes = num_labels + 1  # include blank label
        self.num_units = 50  # = vector dimension of internal states (LSTM size)

    @staticmethod
    def load_data(data_dir, start_time):
        lfps = loadmat(data_dir.joinpath('lfps').as_posix())['lfps']
        event_file = loadmat(data_dir.joinpath('events').as_posix(), chars_as_strings=True, squeeze_me=True)
        times = event_file['times']
        labels = event_file['labels']
        lfp = lfps[start_time:]
        # lfp = np.mean(lfps, 0)[start_time:].reshape(1,-1)
        return lfp, times[times >= start_time] - start_time, labels[times >= start_time]

    def make_batch(self, train=True):
        total_time = self.lfp.shape[0]
        inputs = []
        target_indices = []
        target_values = []
        for i in range(self.batch_size):
            start = np.random.randint(total_time-self.num_time_points+1)
            end = start + self.num_time_points
            inputs.append(self.lfp[start:end,:])
            valid = [start <= t < end for t in self.times]
            target_indices.extend([i, t-start] for t in self.times[valid])
            target_values.extend(self.labels[valid])
        return inputs, (target_indices, target_values)

    def _make_graph(self):
        graph = tf.Graph()
        with graph.as_default():
            inputs = tf.placeholder(dtype=tf.float32, shape=(self.batch_size, self.num_time_points, self.num_features))
            labels = tf.sparse_placeholder(dtype=tf.int32)  # any size
            # these would be used for varying sequence lengths
            sequence_length = tf.constant([self.num_time_points], dtype=tf.int32)

            cell = tf.contrib.rnn.BasicLSTMCell(self.num_units)
            output, _ = tf.nn.dynamic_rnn(cell, inputs, sequence_length=sequence_length, dtype=tf.float32)
            weights = tf.Variable(tf.truncated_normal(shape=(self.batch_size, self.num_units, self.num_classes)))
            bias = tf.Variable(tf.zeros(shape=self.num_classes))

            logits = tf.matmul(output, weights) + bias
            loss = tf.reduce_mean(tf.nn.ctc_loss(labels, logits, sequence_length, time_major=False))
            return graph, inputs, labels, \
                   tf.train.MomentumOptimizer(self.learning_rate, self.momentum).minimize(loss), loss

    def run(self):
        graph, inputs, labels, optimizer, loss = self._make_graph()
        with tf.Session(graph=graph) as sess:
            tf.global_variables_initializer().run()
            for i in range(self.num_epochs):
                batch_inputs = []
                batch_labels = [[],[]]
                while len(batch_labels[1]) < 1:
                    batch_inputs, batch_labels = self.make_batch()
                feed_dict = {
                    inputs: batch_inputs,
                    labels: batch_labels
                }
                _, l = sess.run([optimizer, loss], feed_dict=feed_dict)
                print('Loss at epoch {} = {:5f}'.format(i, l))


if __name__ == '__main__':
    lfpLSTM = LfpLSTM()
    lfpLSTM.run()