# DESI_SANDBOX
This is a package for:
- Create the input file based on the DESI Mock Spectra for the CNN DLA finder
- Evaluate the predicted result from the CNN DLA finder Algorithm

## Installation 
```
git clone https://github.com/JianiDing/DESI_SANDBOX.git
python setup.py develop
```

## Basic Usage 

DLA_meta_build: 
Used to create DLA catalog from the CNN DLA Finder predictions. 
run_command:
Used to run the CNN DLA Finder on the mock spectra
cpm: 
- match_rate: Used to compute the purity and completeness of the predicted results
- confusion_matrix: Used to compute the confusion matrix of the results
desi_meta_v05:
hdf5_writter: Used to write the hdf5 file storing the mock spectra (input of the DLA finder)
undetect_DLA_inspection:
- miss_id_check: Checking the false positive and negative spectra
- plotting: Plotting the false positive and negative spectra



The information for the CNN DLA finder Algorithm is listed here:
https://github.com/davidparks21/qso_lya_detection_pipeline
