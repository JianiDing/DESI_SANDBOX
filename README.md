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


### desi_meta_v05:
- hdf5_writter: Used to write the hdf5 file storing the mock spectra (input of the DLA finder)


### DLA_meta_build: 
- reading_data: Reading Mock DLA meta data from the DESI Mock DLA Catalog
- dla_output_fits: Creating output file for the Mock DLA meta data

### run_command:
- running_data: Run the CNN DLA Finder on the mock spectra


### cpm: 
- dlametadata: Create the DLA meta data from the CNN DLA Finder Predictions
- match_rate: Used to compute the purity and completeness of the predicted results
- confusion_matrix: Used to compute the confusion matrix of the results

### undetect_DLA_inspection:
- miss_id_check: Checking the false positive and false negative spectra
- plotting: Plotting the false positive and false negative spectra



The information for the CNN DLA finder Algorithm is listed here:
https://github.com/davidparks21/qso_lya_detection_pipeline
