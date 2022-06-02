#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np
import h5py
import os
import scipy.signal
import scipy.io as sio
from neurodsp.filt import filter_signal
from neurodsp.spectral import compute_spectrum
from neurodsp.plts.time_series import plot_time_series
from neurodsp.utils import create_times
from fooof import FOOOF
from fooof.bands import Bands
from fooof.analysis import get_band_peak_fm
# import modules from python scripts
from save_FOOOF_results import save_results
from create_FOOOF_reports import create_reports
from create_hdf5 import create_hdf5


def generate_data_reports(file_paths, input_dir, hdf5_dir, report_dir, results_dir,
                          results_file_name, sample_rate, frequency_range):
    """
    :return: new directory with hdf5 files, new directory with PDF FOOOF reports, FOOOF results file
    """
    # returns list of hdf5 file paths for input into FOOOF modeling and generates hdf5 files
    paths = create_hdf5(file_paths, input_dir, hdf5_dir)
    # saving reports as PDFs
    create_reports(paths, sample_rate, frequency_range, report_dir)
    # saving results from FOOOF to text file
    save_results(paths, sample_rate, frequency_range, results_dir, results_file_name)


if __name__ == '__main__':
    # Creating description and arguments for script use on the command line
    parser = argparse.ArgumentParser(description='Give sample rate, frequency range, '
                                                 'input and output file names/paths')
    parser.add_argument('--sr', '--sample_rate',
                        dest='s_rate',
                        type=int,
                        help='sample rate of recording',
                        required=True)
    parser.add_argument('--fr', '--frequency_range',
                        dest='f_range',
                        type=int,
                        help='frequency range for analysis',
                        required=True)
    parser.add_argument('--file_path', '--recording',
                        dest='file',
                        type=str,
                        help='file name or path of recording',
                        required=True)
    parser.add_argument('--in_dir', '--input_directory',
                        dest='in_dir',
                        type=str,
                        help='directory or folder with recording file',
                        required=True)
    parser.add_argument('--hdf5_dir', '--hdf5_directory',
                        dest='hdf5_dir',
                        type=str,
                        help='directory or folder to hold hdf5 files',
                        required=True)
    parser.add_argument('--report_dir', '--report_directory',
                        dest='report_dir',
                        type=str,
                        help='directory or folder to hold PDFs of FOOOF reports',
                        required=True)
    parser.add_argument('--results_dir', '--results_directory',
                        dest='results_dir',
                        type=str,
                        help='directory or folder to hold results from FOOOF report',
                        required=True)
    parser.add_argument('--results_file', '--results_file_name',
                        dest='results_file',
                        type=str,
                        help='file name for holding FOOOF results',
                        required=True)
    # parsing arguments
    args = parser.parse_args()
    generate_data_reports()
