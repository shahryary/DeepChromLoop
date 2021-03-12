#!/usr/bin/env python
__author__ = "Yadollah Shahryary Dizaji"
__description__ = "Deep Learning Project"
__version__ = "1.0.0"
__email__ = "y.shahryary@tum.de"

from configparser import ConfigParser
from keras import backend as K
import time
import h5py
import numpy as np

"""
Configuration & useful functions 
"""
def str2bool(str):
    return str.lower() in ("True", "true", "t", "T", "1")

def config_info_general_data():
    """
    reading GENERAL configuration parameters from config.conf file.
    """
    configParser = ConfigParser()
    configFilePath = 'config.conf'
    configParser.read(configFilePath)
    data_dir = configParser.get('GENERAL', 'data_dir')
    data_samples_to_use = configParser.get('GENERAL', 'data_samples_to_use')
    HiChip_marks = configParser.get('GENERAL', 'HiChip_marks')
    percent_in_test_set = configParser.get('GENERAL', 'percent_in_test_set')
    two_neurons_output = str2bool(configParser.get('GENERAL', 'two_neurons_output'))
    HDF5_samples_num = configParser.get('GENERAL', 'HDF5_samples')
    ATAC_based_negative_samples = str2bool(configParser.get('GENERAL', 'ATAC_based_negative_samples'))

    return data_dir, int(data_samples_to_use), HiChip_marks, int(percent_in_test_set), two_neurons_output, int(HDF5_samples_num), ATAC_based_negative_samples


def config_info_sequence_signals():
    """
    reading 'sequence signals configuration' from config.conf file.
    """
    configParser = ConfigParser()
    configFilePath = 'config.conf'
    configParser.read(configFilePath)

    using_histones = str2bool(configParser.get('SeqSignals', 'histones'))
    using_H3 = str2bool(configParser.get('SeqSignals', 'H3'))
    using_ATAC = str2bool(configParser.get('SeqSignals', 'ATAC'))
    using_fasta_seq = str2bool(configParser.get('SeqSignals', 'fasta_seq'))
    using_methylation_info = str2bool(configParser.get('SeqSignals', 'methylation_info'))

    return using_histones, using_H3, using_ATAC, using_fasta_seq, using_methylation_info


def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))


def data_to_train_test(h5_filename, data_samples_to_use, percent_in_test_set):
    print("start reading h5")
    start_time_start_reading_H5 = time.time()
    print("start_time_start_reading_H5")
    print("--- %s seconds ---" % (time.time() - start_time_start_reading_H5))

    with h5py.File(h5_filename, 'r') as hf:
        re_region_matrices = np.array(hf.get('re_region_matrices'))
        interactions_region_matrices = np.array(hf.get('interactions_region_matrices'))
        labels = np.array(hf.get('labels'))

        re_region_matrices = re_region_matrices[:data_samples_to_use]
        interactions_region_matrices = interactions_region_matrices[:data_samples_to_use]
        labels = labels[:data_samples_to_use]

        # transpose re_region_matrices and interactions_region_matrices
        re_region_matrices = re_region_matrices.transpose(0, 2, 1)
        interactions_region_matrices = interactions_region_matrices.transpose(0, 2, 1)

        num_test_set = int(data_samples_to_use * percent_in_test_set / 100)
        num_train_set = int(data_samples_to_use - num_test_set)

        train_r1 = re_region_matrices[:num_train_set]
        train_r2 = interactions_region_matrices[:num_train_set]
        train_labels = labels[:num_train_set]

        test_r1 = re_region_matrices[-num_test_set:]
        test_r2 = interactions_region_matrices[-num_test_set:]
        test_labels = labels[-num_test_set:]

    print("reading H5 file done")
    print("--- %s seconds ---" % (time.time() - start_time_start_reading_H5))

    return train_r1, train_r2, train_labels, test_r1, test_r2, test_labels
