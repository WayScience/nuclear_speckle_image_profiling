#!/usr/bin/env python
# coding: utf-8

# # Process TIFF files with illumination correction
# 
# Illumination correction (IC) is a step in the image analysis pipeline to remove any illumination errors that could be impacting the biology of the cells. 
# The function we generate in the CellProfiler pipeline are the "best" possible given the parameters available.
# We will process the 16-bit TIFF files with IC and save as the same bit-depth.

# ## Import libraries

# In[1]:


import pathlib
import pprint

import sys

sys.path.append("../utils")
import cp_parallel


# ## Set paths and variables

# In[2]:


# set the run type for the parallelization
run_name = "illum_correction"

# path to IC pipeline
path_to_pipeline = pathlib.Path("./pipelines/illumination_correction.cppipe").resolve(
    strict=True
)

# set main output dir for all plates
output_dir = pathlib.Path("./IC_corrected_images")
output_dir.mkdir(exist_ok=True)

# directory where images are located within folders
images_dir = pathlib.Path("../1.max_projection/max_projected_images")

# list for plate names based on folders to use to create dictionary
plate_names = []
# iterate through images directory and append plate names from folder names that contain image data from that plate
for file_path in images_dir.iterdir():
    if str(file_path.stem).startswith("slide"):
        plate_names.append(str(file_path.stem))

print(plate_names)


# ## Create dictionary with all info for each plate

# In[3]:


# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    name: {
        "path_to_images": pathlib.Path(list(images_dir.rglob(name))[0]).resolve(
            strict=True
        ),
        "path_to_output": pathlib.Path(f"{output_dir}/{name}"),
        "path_to_pipeline": path_to_pipeline,
    }
    for name in plate_names
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run illumination correction pipeline on each plate in parallel
# 
# In this notebook, we do not run the cells to completion as we prefer to run the notebooks as nbconverted python files due to better stability.

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)

