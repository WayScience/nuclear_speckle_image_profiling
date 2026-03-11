#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite file output from CellProfiler into parquet file using Cytotable

# ## Import libraries

# In[1]:


import logging
import pathlib

import pandas as pd

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)


# ## Set paths and variables

# In[2]:


# type of file output for CytoTable
dest_datatype = "parquet"

# set main output dir for all parquet files
output_dir = pathlib.Path("./data")
output_dir.mkdir(exist_ok=True)

# directory where SQLite files are located
sqlite_dir = pathlib.Path("../3.cp_analysis/analysis_output").resolve(strict=True)

# Set converted parquet dir
parquet_dir = pathlib.Path(f"{output_dir}/converted_profiles")
parquet_dir.mkdir(exist_ok=True)

# Set plate names as an empty list to append to
plate_names = []

# directory with plate maps
platemap_dir = pathlib.Path(f"../0.download_data/metadata/platemaps")

# list for plate names based on metadata files to use to create dictionary
plate_names = []
# iterate through metadata dir and append plate names from metadata files
for file_path in platemap_dir.iterdir():
    filename = file_path.stem
    first_index = filename.split("_")[0]
    plate_names.append(first_index)

# print the plate names and how many plates there are (confirmation)
print(f"There are {len(plate_names)} plates in this dataset. Below are the names:")
for name in plate_names:
    print(name)


# ## Run cytotable convert to output nuclei and image features separately for all plates

# In[3]:


# Iterate over directory with SQLite outputs
for plate_folder in sqlite_dir.iterdir():
    # Using the plate names list, only process files within that list
    if plate_folder.name in plate_names:
        # Construct output path for converted parquet file
        output_path = pathlib.Path(f"{parquet_dir}/{plate_folder.stem}/{plate_folder.stem}_converted.parquet")
        
        print("Starting conversion with cytotable for plate:", plate_folder.stem)

        # merge single cells and output as parquet file
        convert(
            source_path=str(plate_folder),
            dest_path=str(output_path),
            dest_datatype=dest_datatype,
            metadata=["image"],
            compartments=["nuclei"],
            identifying_columns=["ImageNumber"],
            joins="""
            SELECT
                *
            FROM
                read_parquet('per_image.parquet') as per_image
            INNER JOIN read_parquet('per_nuclei.parquet') AS per_nuclei ON
                per_nuclei.Metadata_ImageNumber = per_image.Metadata_ImageNumber
            """,
            chunk_size=10000,
        )

        print("Conversion finished for plate:", plate_folder.stem)


# ## Remove unwanted image + metadata columns and split the bulk and single-cell data from the main parquet file

# In[4]:


# path to unwanted image cols text file
unwanted_list_path = pathlib.Path("./unwanted_image_cols.txt")
# Load the list of columns to remove from the text file
with open(unwanted_list_path, "r") as file:
    columns_to_remove = [line.strip() for line in file]

# Iterate through directory with converted outputs
for plate_folder in parquet_dir.iterdir():
    # Only process the files that are in the plate names list
    if plate_folder.name in plate_names:
        # Read in each file as data frame
        plate_df = pd.read_parquet(
            pathlib.Path(f"{plate_folder}/{plate_folder.stem}_converted.parquet")
        )
        print(
            "Starting to edit image and nuclei data frames for plate:",
            plate_folder.stem,
        )

        # Drop the specified columns (ignore error if a column isn't there)
        plate_df = plate_df.drop(columns=columns_to_remove, errors="ignore")

        # Identify metadata columns for nuclei data frame
        metadata_columns = [
            "Metadata_ImageNumber",
            "Image_Metadata_Plate",
            "Image_Metadata_Site",
            "Image_Metadata_Well",
            "Image_Count_Nuclei",
            "Image_FileName_DAPI"
        ]

        # Create nuclei (single-cell) data frame
        nuclei_df = plate_df[
            metadata_columns
            + [col for col in plate_df.columns if col.startswith("Nuclei_")]
        ]

        # Create image (bulk) data frame
        image_df = plate_df[
            ["Metadata_ImageNumber"]
            + [col for col in plate_df.columns if col.startswith("Image_")]
        ]
        # Drop duplicate images in the image data frame since each image will have the same values even if the row is repeated
        image_df = image_df.drop_duplicates(subset="Metadata_ImageNumber")

        # Save nuclei and image data frames to the same folder as the plate
        nuclei_df.to_parquet(f"{plate_folder}/per_nuclei.parquet", index=False)
        image_df.to_parquet(f"{plate_folder}/per_image.parquet", index=False)

        # nuclei_df and image_df shape and one data frame to assess all looks correct
        print("Shape of nuclei data frame", nuclei_df.shape)
        print("Shape of image data frame", image_df.shape)

