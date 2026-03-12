#!/usr/bin/env python
# coding: utf-8

# # Perform data preprocessing with pycytominer on single cell features
# 
# Note: Single cell is represented by only nuclei compartment features which was used to extract features across all three channels.

# ## Import libraries

# In[1]:


import pathlib
import pprint

import pandas as pd
from pycytominer import annotate, normalize, feature_select


# ## Set paths and variables

# In[2]:


# Path to dir with converted profiles per plate (each plate as a folder)
cleaned_dir = pathlib.Path("./data/cleaned_sc_profiles")

# path for plate map directory
platemap_dir = pathlib.Path(f"../0.download_data/metadata/platemaps")

# Output dir for the files to be saved to
output_dir = pathlib.Path("./data/single_cell_profiles")
output_dir.mkdir(exist_ok=True, parents=True)

# Extract the plate names from the file names
plate_names = [file.stem.split("_")[0] for file in platemap_dir.glob("*_platemap.csv")]

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
]

# print the plate names and how many plates there are (confirmation)
print(f"There are {len(plate_names)} plates in this dataset. Below are the names:")
for name in plate_names:
    print(name)


# ## Create dictionary with unique paths for each plate

# In[3]:


# create plate info dictionary
plate_info_dictionary = {
    name: {
        "profile_path": str(
            pathlib.Path(f"{cleaned_dir}/{name}_sc_cleaned.parquet").resolve(strict=True)
        ),
        "platemap_path": str(
            pathlib.Path(f"{platemap_dir}/{name}_platemap.csv").resolve(strict=True)
        ),
    }
    for name in plate_names
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Perform preprocessing on single cell features

# In[4]:


for plate, info in plate_info_dictionary.items():
    print(f"Performing pycytominer pipeline for {plate}")
    # Set output paths per preprocessing step
    output_annotated_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_annotated.parquet")
    )
    output_normalized_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_normalized.parquet")
    )
    output_feature_select_file = str(
        pathlib.Path(f"{output_dir}/{plate}_sc_feature_selected.parquet")
    )

    # Load in the converted profile to be used in the first step
    profile_df = pd.read_parquet(info["profile_path"])

    # Load in platemap file with most relevant columns for annotation
    platemap_df = pd.read_csv(info["platemap_path"], usecols=["Well", "CellLine", "Condition"])

    # Step 1: Annotation
    annotate(
        profiles=profile_df,
        platemap=platemap_df,
        join_on=["Metadata_Well", "Image_Metadata_Well"],
        output_file=output_annotated_file,
        output_type="parquet",
    )

    # Load the annotated parquet file to fix metadata columns names
    annotated_df = pd.read_parquet(output_annotated_file)

    # Rename columns using the rename() function
    column_name_mapping = {
        "Image_Metadata_Site": "Metadata_Site",
        "Image_Count_Nuclei": "Metadata_Nuclei_Site_Count",
        "Nuclei_Location_Center_X": "Metadata_Nuclei_Location_Center_X",
        "Nuclei_Location_Center_Y": "Metadata_Nuclei_Location_Center_Y",
    }

    annotated_df.rename(columns=column_name_mapping, inplace=True)

    # Save the modified DataFrame back to the same location
    annotated_df.to_parquet(output_annotated_file, index=False)

    # Step 2: Normalization
    normalized_df = normalize(
        profiles=output_annotated_file,
        method="standardize",
        output_file=output_normalized_file,
        output_type="parquet",
    )

    # Step 3: Feature selection
    feature_select(
        output_normalized_file,
        operation=feature_select_ops,
        output_file=output_feature_select_file,
        output_type="parquet",
    )
    print(
        f"Annotation, normalization, and feature selection have been performed for {plate}"
    )


# ## Check example output file to confirm that the process worked

# In[5]:


# Check output file
test_df = pd.read_parquet("./data/single_cell_profiles/slide1_sc_annotated.parquet")

print(test_df.shape)
test_df.head(2000)

