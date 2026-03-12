#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the main conda environment
conda activate nuclear_speckle_cp_env

# convert notebook to python file into the scripts folder
jupyter nbconvert --to script --output-dir=scripts/ *.ipynb

# run python script to analyze plates
python scripts/cp_analysis.py
