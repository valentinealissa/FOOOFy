#!/usr/bin/env python
from neurodsp.spectral import compute_spectrum
import scipy.signal
from fooof import FOOOF
from neurodsp.plts.time_series import plot_time_series
from neurodsp.utils import create_times
from neurodsp.filt import filter_signal
import numpy as np
import h5py
import os

# sample values for file paths, sample rate, frequency range, and output directory
paths_to_files = ['/Users/alissavalentine/ephys_rotation/hdf5/PT-0017_STNLeft_Pass1_R1C_Snapshot-30.0sec2_R-3.71.hdf5',
                  '/Users/alissavalentine/ephys_rotation/hdf5/PT-0017_STNLeft_Pass1_R1C_Snapshot-30.0sec1_R-3.81.hdf5']
sr = 48000
f_range = (2, 50)
output_dir = '/Users/alissavalentine/ephys_rotation/hdf5/reports/'


def create_reports(file_paths, sample_rate, frequency_range, output_directory):
    # change directory to where you want the reports saved
    os.chdir(output_directory)

    # iterate through each file in the list of file paths
    for file in file_paths:
        # read the h5py file
        with h5py.File(file, "r") as f:
            # edit future file name to remove dots since it caused error
            file_name = file.split('.hdf5')[0].replace(".", "_")
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
            fm.save_report(file_name=file_name, file_path='/Users/alissavalentine/ephys_rotation/hdf5/reports/')


#create_reports(paths_to_files, sr, f_range, output_dir)
