#!/usr/bin/env python
import scipy.io as sio
import pandas as pd
import h5py
import os

# sample values for directories and metadata file
input_dir = '/Users/alissavalentine/ephys_rotation'
data_file = 'qc_depth_duration.csv'
output_dir = '/Users/alissavalentine/ephys_rotation/hdf5'


def create_hdf5(meta_data_file, input_file_directory, output_file_directory):
    # change directory to where metadata file is
    os.chdir(input_file_directory)

    # load csv file with metadata including path names for files
    meta_data = pd.read_csv(meta_data_file)

    # creating empty list to hold file names
    file_paths = []

    # changing directory to where you want the output file saved
    os.chdir(output_file_directory)

    # iterating through each file saving files in h5py
    for path in meta_data.file_path:
        # loading data from mat lab file
        mat_data = sio.loadmat(path)
        # extracting raw signal data from file
        key = list(mat_data.keys())[3]
        data_structure = mat_data[key]  # returns numpy nd array with values within the data structure A,b,c etc.
        data_structure = data_structure[0][0]
        continuous_channel = data_structure[-3]
        continuous_channel = continuous_channel[0][0]
        data = continuous_channel[-1]
        # creating new path and file name for hdf5 version
        new_path = path.replace("MATfiles", "HDF5files").split('.mat')[0] + '.hdf5'
        new_path = path.split('/')[-1].split('.mat')[0] + '.hdf5'
        # saving path/file name to empty list
        file_paths.append(new_path)
        # writing data to h5py file and closing the file
        f = h5py.File(new_path, "w")
        f.create_dataset('data', data=data[:, 0], dtype='<i2')
        f.close()
    return file_paths


# paths = create_hdf5(data_file, input_dir, output_dir)
# print(paths)
