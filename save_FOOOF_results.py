#!/usr/bin/env python
from neurodsp.spectral import compute_spectrum
import scipy.signal
from neurodsp.filt import filter_signal
import numpy as np
import h5py
import os
from fooof import FOOOF
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fm

# values included:
# file_name
# theta_CF,theta_PW,theta_BW
# alpha_CF,alpha_PW,alpha_BW
# beta_CF,beta_PW,beta_BW
# lgamma_CF,lgamma_PW,lgamma_BW
# Rsquare,error
# exponent,offset

# sample values for sample rate, frequency range, output directory, and file paths
paths_to_files = ['/Users/alissavalentine/ephys_rotation/hdf5/PT-0017_STNLeft_Pass1_R1C_Snapshot-30.0sec2_R-3.71.hdf5',
                  '/Users/alissavalentine/ephys_rotation/hdf5/PT-0017_STNLeft_Pass1_R1C_Snapshot-30.0sec1_R-3.81.hdf5']
sr = 48000
f_range = (2, 50)
output_dir = '/Users/alissavalentine/ephys_rotation/'
name = "test.txt"


def save_results(file_paths, sample_rate, frequency_range, output_directory, output_file_name):
    # empty list to hold future values from FOOOF
    final_values = []

    # iterate through files
    for file in file_paths:
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
            recording_name = file.rsplit("/")[-1].rstrip('.hdf5')
            # saving all extracted values into a row
            new_data = str(recording_name + "," + str(theta_CF) + "," + str(theta_PW) + "," + str(theta_BW) + "," + str(
                alpha_CF) + "," + str(alpha_PW) + "," + str(alpha_BW) + "," + str(beta_CF) + "," + str(
                beta_PW) + "," + str(beta_BW) + "," + str(low_gamma_CF) + "," + str(low_gamma_PW) + "," + str(
                low_gamma_BW) + "," + str(r_squared) + "," + str(error) + "," + str(exponent) + "," + str(offset))
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


#save_results(paths_to_files, sr, f_range, output_dir, name)
