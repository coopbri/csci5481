# CSCI 5481 Final Project
# University of Minnesota
# Code by Brian Cooper
# This file processes the input data and prepares it for experimentation.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import expm
from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler

# ============================================================================ #
# Data preparation - PCA dimensionality reduction, etc.                        #
# ============================================================================ #
def prepare_data(data, target, n, training_size, test_size, class_labels):
    # Split data for training and testing
    sample_train, sample_test, label_train, label_test = \
        train_test_split(data, target, test_size=0.3, random_state=12)

    # Standarize for Gaussian around 0 with unit variance
    std_scale = StandardScaler().fit(sample_train)
    sample_train = std_scale.transform(sample_train)
    sample_test = std_scale.transform(sample_test)

    # Reduce number of features to number of qubits
    pca = PCA(n_components=n).fit(sample_train)
    sample_train = pca.transform(sample_train)
    sample_test = pca.transform(sample_test)

    # Scale to range (-1, +1) for convenience
    samples = np.append(sample_train, sample_test, axis=0)
    minmax_scale = MinMaxScaler((-1, 1)).fit(samples)
    sample_train = minmax_scale.transform(sample_train)
    sample_test = minmax_scale.transform(sample_test)

    # Pick training size samples from each distribution
    training_input = {key: (sample_train[label_train == k, :])[:training_size] \
        for k, key in enumerate(class_labels)}
    test_input = {key: (sample_train[label_train == k, :])[training_size:(
        training_size+test_size)] for k, key in enumerate(class_labels)}

    return sample_train, label_train, training_input, test_input

# ============================================================================ #
# Ecoli dataset                                                                #
# ============================================================================ #
def ecoli(training_size, test_size, n, PLOT_DATA):
    class_labels = [r'cp', r'im', r'pp', r'imU', r'om', r'omL', r'imL', r'imS']

    df = pd.read_csv('datasets/ecoli.csv', header=None)
    df = df.replace({'cp':0, 'im':1, 'pp':2, 'imU':3, 'om':4, 'omL':5, \
        'imL':6, 'imS':7})
    data = df.iloc[:,1:8].astype(float)
    target = df.iloc[:,8]

    sample_train, label_train, training_input, test_input = \
        prepare_data(data, target, n, training_size, test_size, class_labels)

    if PLOT_DATA:
        for k in range(8):
            if k == 0:
                label = 'cp'
            elif k == 1:
                label = 'im'
            elif k == 2:
                label = 'pp'
            elif k == 3:
                label = 'imU'
            elif k == 4:
                label = 'om'
            elif k == 5:
                label = 'omL'
            elif k == 6:
                label = 'imL'
            else:
                label = 'imS'

            plt.scatter(sample_train[label_train == k, 0][:training_size],
                        sample_train[label_train == k, 1][:training_size],
                        label=label)

        plt.title("E. Coli Dataset")
        plt.legend()
        plt.show()

    return sample_train, training_input, test_input, class_labels

# ============================================================================ #
# Yeast dataset                                                                #
# ============================================================================ #
def yeast(training_size, test_size, n, PLOT_DATA):
    class_labels = [r'CYT', r'NUC', r'MIT', r'ME3', r'ME2', r'ME1', r'EXC', \
        r'VAC', r'POX', r'ERL']

    df = pd.read_csv('datasets/yeast.csv', header=None)
    df = df.replace({'CYT':0, 'NUC':1, 'MIT':2, 'ME3':3, 'ME2':4, 'ME1':5, \
        'EXC':6, 'VAC':7, 'POX':8, 'ERL':9})
    data = df.iloc[:,1:9]
    target = df.iloc[:,9]

    sample_train, label_train, training_input, test_input = \
        prepare_data(data, target, n, training_size, test_size, class_labels)

    if PLOT_DATA:
        for k in range(10):
            if k == 0:
                label = 'CYT'
            elif k == 1:
                label = 'NUC'
            elif k == 2:
                label = 'MIT'
            elif k == 3:
                label = 'ME3'
            elif k == 4:
                label = 'ME2'
            elif k == 5:
                label = 'ME1'
            elif k == 6:
                label = 'EXC'
            elif k == 7:
                label = 'VAC'
            elif k == 8:
                label = 'POX'
            else:
                label = 'ERL'

            plt.scatter(sample_train[label_train == k, 0][:training_size],
                        sample_train[label_train == k, 1][:training_size],
                        label=label)

        plt.title("Yeast Dataset")
        plt.legend()
        plt.show()

    return sample_train, training_input, test_input, class_labels

# ============================================================================ #
# Mouse dataset                                                                #
# ============================================================================ #
def mouse(training_size, test_size, n, PLOT_DATA):
    class_labels = [r'c-CS-s', r'c-CS-m', r'c-SC-s', r'c-SC-m', \
        r't-CS-s', r't-CS-m', r't-SC-s', r't-SC-m']

    df = pd.read_csv('datasets/mouse.csv', header=None)
    df = df.replace({'c-CS-s':0, 'c-CS-m':1, 'c-SC-s':2, 'c-SC-m':3, \
        't-CS-s':4, 't-CS-m':5, 't-SC-s':6, 't-SC-m':7})
    data = df.iloc[:,1:10]
    target = df.iloc[:,78]

    sample_train, label_train, training_input, test_input = \
        prepare_data(data, target, n, training_size, test_size, class_labels)

    if PLOT_DATA:
        for k in range(8):
            if k == 0:
                label = 'c-CS-s'
            elif k == 1:
                label = 'c-CS-m'
            elif k == 2:
                label = 'c-SC-s'
            elif k == 3:
                label = 'c-SC-m'
            elif k == 4:
                label = 't-CS-s'
            elif k == 5:
                label = 't-CS-m'
            elif k == 6:
                label = 't-SC-s'
            else:
                label = 't-SC-m'

            plt.scatter(sample_train[label_train == k, 0][:training_size],
                        sample_train[label_train == k, 1][:training_size],
                        label=label)

        plt.title("Mouse Dataset")
        plt.legend()
        plt.show()

    return sample_train, training_input, test_input, class_labels

# ============================================================================ #
# RNA-seq dataset                                                              #
# ============================================================================ #
# def rnaseq(training_size, test_size, n, PLOT_DATA):
#     class_labels = [r'PRAD', r'LUAD', r'BRCA', r'KIRC', r'COAD']

#     df = pd.read_csv('datasets/rnaseq.csv', header=None)
#     df = df.replace({'PRAD':0, 'LUAD':1, 'BRCA':2, 'KIRC':3, 'COAD':4})
#     data = df.iloc[:,1:20532]
#     target = df.iloc[:,:20532]

#     sample_train, label_train, training_input, test_input = \
#         prepare_data(data, target, n, training_size, test_size, class_labels)

#     if PLOT_DATA:
#         for k in range(5):
#             if k == 0:
#                 label = 'PRAD'
#             elif k == 1:
#                 label = 'LUAD'
#             elif k == 2:
#                 label = 'BRCA'
#             elif k == 3:
#                 label = 'KIRC'
#             else:
#                 label = 'COAD'

#             plt.scatter(sample_train[label_train == k, 0][:10],
#                         sample_train[label_train == k, 1][:10],
#                         label=label)

#         plt.title("RNA-seq Dataset")
#         plt.legend()
#         plt.show()

#     return sample_train, training_input, test_input, class_labels
