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
    parser = argparse.ArgumentParser(description='Give # sentences added, % words replaced')
    parser.add_argument('--s', '--sentence',
                        dest='s_add',
                        type=int,
                        help='# sentences added',
                        required=False,
                        default=2)
    # parsing arguments
    args = parser.parse_args()
    generate_data_reports()
