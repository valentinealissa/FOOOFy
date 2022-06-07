#!/usr/bin/env python
# Last edited: 6/7/22
# Alissa Valentine
# This script requires the following inputs:
# 1. EDF file
# 2. Output path for hdf5 converted files
# 3. Frequency range
# 4. Sample rate
# 5. Directory path for reports PDFs (if running create_reports)
# 6. Directory path for results file (if running save_results)
# 7. File name for results file (if running save_results)

import numpy as np
import h5py
import os
import pyedflib
from neurodsp.spectral import compute_spectrum
import scipy.signal
from fooof import FOOOF
from neurodsp.plts.time_series import plot_time_series
from neurodsp.utils import create_times
from neurodsp.filt import filter_signal
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fm


def edf_to_hdf5(input_file_name, output_path):
    """
    :param input_file_name: the edf file with raw signal data
    :param output_path: the location where you want hdf5 files saved
    :return: a list of the file paths for hdf5 files, and an hdf5 file for each electrode recording
    """
    # making list of electrodes to not include
    ignore = ["Pleth", "PR", "OSAT", "TRIG", "DC1", "DC2", "DC3", "DC4", "DC5", "DC6", "DC7", "DC8", "DC9", "DC10",
              "DC11", "DC12", "DC13", "DC14", "DC15", "DC16"]
    # read EDF file
    f = pyedflib.EdfReader(input_file_name)
    # extract header
    header = f.getSignalHeaders()
    # create empty list to hold file paths
    metadata = []
    # iterate through the header
    for i, h in enumerate(header):
        # print the header to look at the info available (currently commented out)
        #print(h)
        # save electrode name
        electrode = list(h.values())[0]
        # iterate through the worthwhile electrodes
        if electrode not in ignore:
            # extract the raw signal data
            data = f.readSignal(i, 0)
            # save patient file name info
            patient_id = input_file_name.rsplit(".")[0]
            # create output file name with patient info and electrode
            output_file_name = "_".join([patient_id, electrode])
            # add path to the output file name
            output_path_name = output_path + "/" + output_file_name + '.hdf5'
            # save the path to list
            metadata.append(output_path_name)
            # open a h5py file and write the data to it then close it
            output = h5py.File(output_path_name, "w")
            output.create_dataset('data', data=data[:], dtype='<i2')
            output.close()
    f.close()
    # return the file paths
    return metadata


def create_reports(metadata_paths, frequency_range, sample_rate, output_directory):
    """
    :param metadata_paths: list of hdf5 file paths
    :param frequency_range: tuple containing lowest and highest frequency ex:(2,50)
    :param sample_rate: integer expressing the sample rate of the recordings
    :param output_directory: where you want the reports saved
    :return: FOOOF report for each electrode in edf file
    """
    # change directory to where you want the reports saved
    os.chdir(output_directory)

    # iterate through each file in the list of file paths
    for file_path in metadata_paths:
        # edit future file name to remove dots (since it caused error) and sample rate
        file_name = file_path.split('.hdf5')[0].rsplit('/')[-1].split('_64')[0]
        # read the h5py file
        with h5py.File(file_path, "r") as f:
            # save raw data as np.float64
            raw_signal = np.array(f.get('data'), dtype=np.float64)
            # bandpass filter signal using sample rate and f range
            filtered_signal = filter_signal(raw_signal, sample_rate, 'bandpass', frequency_range)
            # remove nan values from filtered signal
            filtered_signal = np.array([x for x in filtered_signal if np.isnan(x) == False])
            # downsample the data
            secs = len(filtered_signal) / sample_rate  # Number of seconds in signal X
            samps = int(secs * 1000)  # Number of samples to downsample
            down_filtered_signal = scipy.signal.resample(filtered_signal, samps)
            # optional plots to see the filtered data (uncomment them to see them)
            # times = create_times((len(filt_signal)/sr), sr)
            # plot_time_series(times, raw_signal, ['Raw'])
            # plot_time_series(times, filt_signal, ['Filtered'])
            # compute power spectrum analysis with welch method and given f range
            freqs, powers = compute_spectrum(down_filtered_signal, 1000, f_range=frequency_range)
            # Initialize a FOOOF object
            # manually change parameters as needed
            fm = FOOOF(peak_width_limits=[2, 7], min_peak_height=0.05, max_n_peaks=3, verbose=True)
            # fit FOOOF model to power spectrum data
            fm.fit(freqs, powers)
            # generate report and save the report to output directory
            fm.report(freqs, powers, [5, 45])
            # update the file path to what you want
            fm.save_report(file_name=file_name, file_path='/Users/alissavalentine/ephys_rotation/hdf5/reports/')


def save_results(metadata_paths, frequency_range, sample_rate, output_directory, output_file_name):
    """
    :param metadata_paths: list of hdf5 file paths
    :param frequency_range: tuple containing lowest and highest frequency ex:(2,50)
    :param sample_rate: integer expressing the sample rate of the recordings
    :param output_directory: where you want the FOOOF results file saved
    :param output_file_name: name of the FOOOF results file
    :return: txt file with parameterized FOOOF values for each electrode in recording/edf file
    """
    # empty list to hold future values from FOOOF
    final_values = []
    # iterate through files
    for file in metadata_paths:
        # read h5py file
        with h5py.File(file, "r") as f:
            # save raw data as np.float64
            raw_signal = np.array(f.get('data'), dtype=np.float64)
            # bandpass filter signal using sample rate and f range
            filtered_signal = filter_signal(raw_signal, sample_rate, 'bandpass', frequency_range)
            # remove nan values from filtered signal
            filtered_signal = np.array([x for x in filtered_signal if np.isnan(x) == False])
            # downsample the data
            secs = len(filtered_signal) / sample_rate  # Number of seconds in signal X
            samps = int(secs * 1000)  # Number of samples to down sample
            down_filtered_signal = scipy.signal.resample(filtered_signal, samps)
            # compute power spectrum analysis with welch method and given f range
            freqs, powers = compute_spectrum(down_filtered_signal, 1000, f_range=frequency_range)
            # setting the frequency ranges for the bands
            bands = Bands({'theta': [4, 8], 'alpha': [8, 13], 'beta': [13, 30], 'low_gamma': [30, 50]})
            # Initialize a FOOOF object
            fm = FOOOF(peak_width_limits=[2, 7], min_peak_height=0.05, max_n_peaks=3, verbose=True)
            # fit FOOOF model to power spectrum data
            fm.fit(freqs, powers)
            # extracting largest peak information from theta range
            theta = get_band_peak_fm(fm, bands.theta)
            theta_CF = theta[0]
            theta_PW = theta[1]
            theta_BW = theta[2]
            # extracting largest peak information from alpha range
            alpha = get_band_peak_fm(fm, bands.alpha)
            alpha_CF = alpha[0]
            alpha_PW = alpha[1]
            alpha_BW = alpha[2]
            # extracting largest peak information from beta range
            beta = get_band_peak_fm(fm, bands.beta)
            beta_CF = beta[0]
            beta_PW = beta[1]
            beta_BW = beta[2]
            # extracting largest peak information from low gamma range
            low_gamma = get_band_peak_fm(fm, bands.low_gamma)
            low_gamma_CF = low_gamma[0]
            low_gamma_PW = low_gamma[1]
            low_gamma_BW = low_gamma[2]
            # extracting full model performance values
            error = fm.get_params('error')
            r_squared = fm.get_params('r_squared')
            # extracting aperiodic fit values
            exponent = fm.get_params('aperiodic_params')[1]
            offset = fm.get_params('aperiodic_params')[0]
            # saving recording name
            recording_name = file.split('.hdf5')[0].rsplit('/')[-1].split('_64')[0]
            # saving all extracted values into a row
            new_data = [str(recording_name), str(theta_CF), str(theta_PW), str(theta_BW), str(alpha_CF), str(alpha_PW),
                        str(alpha_BW), str(beta_CF), str(beta_PW), str(beta_BW), str(low_gamma_CF), str(low_gamma_PW),
                        str(low_gamma_BW), str(r_squared), str(error), str(exponent), str(offset)]
            new_data = ",".join(new_data)
            # adding values to the final values list
            final_values.append(new_data)

    # change directory to where you want the reports saved
    os.chdir(output_directory)

    # save results in text file
    # open file for all results
    all_file = open(output_file_name, "w")
    with all_file as file:
        # creating header for the file
        file.writelines("file_name,theta_CF,theta_PW,theta_BW,"
                        "alpha_CF,alpha_PW,alpha_BW,"
                        "beta_CF,beta_PW,beta_BW,"
                        "lgamma_CF,lgamma_PW,lgamma_BW,"
                        "Rsquare,error,exponent,offset\n")
        # adding each line of data
        for line in final_values:
            file.writelines(str(line))
            file.writelines("\n")
    all_file.close()


# setting values, update to your needs
in_file = "MS007_3PreMed2.edf"
path = '/Users/alissavalentine/ephys_rotation/hdf5'
f_range = (2, 50)
sr = float(1024)
report_dir = '/Users/alissavalentine/ephys_rotation/hdf5/reports/'
results_dir = '/Users/alissavalentine/ephys_rotation/'
results_file = "MS007_3PreMed2_FOOOF.txt"

# Running functions
file_paths = edf_to_hdf5(in_file, path)
#create_reports(file_paths, f_range, sr, report_dir)
save_results(file_paths, f_range, sr, results_dir, results_file)
