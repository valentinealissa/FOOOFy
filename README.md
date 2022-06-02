# Parameter Extraction with FOOOF
> This folder contains scripts for generating values from neural 
> ephys recordings. The following steps should be followed when 
> using multiple matlab structured data files.
> 
## Implementation
### Modules required:
* `h5py`: https://docs.h5py.org/en/latest/build.html
* `neuroDSP`: https://neurodsp-tools.github.io/neurodsp/
* `scipy`: https://docs.scipy.org/doc/scipy/getting_started.html
* `FOOOF`: https://fooof-tools.github.io/fooof/index.html
### Files required:
These files should be run in the numerical order seen below
1. `file_agg.m`: aggregates file metadata (electrode depth, recording duration, file path) into .mat and .csv 
2. `file_qc.m`: quality controls .mat file for electrode depth and recording duration. creates new .mat and .csv file
3. **`FOOOF.py`: runs `create_hdf5.py`, `create_FOOOF_reports.py`, and `save_FOOOF_results.py` in sequential order
   1. `create_hdf5.py`: writes matlab structured data into hdf5 files
   2. `create_FOOOF_reports.py`: generates PDF reports of FOOOF model fitting for each recording
   3. `save_FOOOF_results.py`: generates txt file with extracted values from FOOOF model fitting

>** This is the final script that saves report PDFs and results data
### Running:
* In MATLAB:
  * Run `file_agg.m`
  * Run `file_qc.m`
* In command line:
  * Run `FOOOF.py` using text below:
> ```angular2html
> python ./FOOOF.py \
> --sr = 48000 \
> --fr = 2 50 \
> --in_file = qc_depth_duration_mini.csv \
> --in_dir = /Users/alissavalentine/ephys_rotation \
> --hdf5_dir = /Users/alissavalentine/ephys_rotation/hdf5 \
> --report_dir = /Users/alissavalentine/ephys_rotation/hdf5/reports/ \
> --results_dir = /Users/alissavalentine/ephys_rotation/ \
> --results_file = results_mini.txt
>```
