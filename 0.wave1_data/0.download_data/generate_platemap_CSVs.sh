#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the main conda environment
conda activate nuclear_speckle_cp_env

# change directory into the metadata folder
cd ./metadata

# convert notebook to python file into the nbconverted folder
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run python script to analyze plates
python nbconverted/generate_platemap_csvs.py

# change directory back to original
cd ..
echo "Platemap CSV files have been generated!"
