#!/usr/bin/env python
# coding: utf-8

# # Generate platemap files from original position files
# 
# In this notebook, position files are used to generate plate map files to use for saving max projection images and annotating single-cell outputs downstream. Platemap files are saved which are only one row per Well with relevant perturbation and cell line info.
# 
# Intermediate files called `cellprofiler_csvs` are created to use for outputting TIFFs from the nd2 file in the next module (1.max_projection).

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd


# ## Generate intermediate CellProfiler CSV files for metadata module

# In[2]:


# Dir path for output of platemap CSV files
cp_csv_dir = pathlib.Path("./cellprofiler_csvs")
cp_csv_dir.mkdir(parents=True, exist_ok=True)

# Find all position txt files in the current directory starting with "slide" using glob
position_files = pathlib.Path().resolve().glob('slide*')

# Instantiate a empty list to append cellprofiler csvs to
cp_csv_dfs = []

# Iterate through each file to update "Point Name" and "Image" columns
for file in position_files:
    # Read the CSV file
    df = pd.read_csv(file, delimiter='\t', encoding='utf-16')
    
    # Remove '#' prefix from 'Point Name' column
    df['Point Name'] = df['Point Name'].str.lstrip('#')
    
    # Zero-index the 'Image' column
    df['Image'] = df['Image'] - 1
    
    # Save the processed DataFrame to the cellprofiler csvs directory
    output_file = pathlib.Path(f"{cp_csv_dir}/{file.stem}.csv")
    df.to_csv(output_file, index=False)
    
    # Append the processed DataFrame to the list
    cp_csv_dfs.append(df)

# Print the list of dataframes to verify that the process worked
for df in cp_csv_dfs:
    print(df.head())


# ## Generate platemap files

# In[5]:


# Dir path for output of platemap CSV files
platemap_dir = pathlib.Path("./platemaps")
platemap_dir.mkdir(parents=True, exist_ok=True)

# Find all position txt files in the current directory starting with "slide" using glob
position_files = pathlib.Path().resolve().glob('slide*')

# Instantiate an empty list to append platemaps to
platemap_dfs = []

# Iterate through each file to update and reduce the rows to one per well
for file in position_files:
    # Read the CSV file
    df = pd.read_csv(file, delimiter='\t', encoding='utf-16')
    
    # Only keep relevant columns to perturbation and cell line
    df = df[['Well', 'CellLine', 'Condition']]
    
    # Reduce rows down to one per well
    df = df.drop_duplicates(subset='Well')
    
    # Save the processed DataFrame to the platemap directory
    output_file = pathlib.Path(f"{platemap_dir}/{file.stem.split('.')[0]}_platemap.csv")
    df.to_csv(output_file, index=False)
    
    # Append the processed DataFrame to the list
    platemap_dfs.append(df)

# Print the list of dataframes to verify that the process worked
for df in platemap_dfs:
    print(df)


# In[ ]:




