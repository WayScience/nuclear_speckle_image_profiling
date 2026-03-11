#!/usr/bin/env python
# coding: utf-8

# # Generate UMAP coordinates for each plate

# ## Import libraries

# In[1]:


import glob
import pathlib
import pandas as pd
import umap

from pycytominer import feature_select
from pycytominer.cyto_utils import infer_cp_features


# ## Set constants

# In[2]:


# Set constants
umap_random_seed = 0
umap_n_components = 2

output_dir = pathlib.Path("results")
output_dir.mkdir(parents=True, exist_ok=True)


# ## Create list of paths to feature selected data per plate

# In[3]:


# Set input paths
data_dir = pathlib.Path("../4.preprocess_features/data/single_cell_profiles")

# Select only the feature selected files
file_suffix = "*sc_feature_selected.parquet"

# Obtain file paths for all feature selected plates
fs_files = glob.glob(f"{data_dir}/{file_suffix}")
fs_files


# In[4]:


# Load feature data into a dictionary, keyed on plate name
cp_dfs = {x.split("/")[-1]: pd.read_parquet(x) for x in fs_files}

# Print out useful information about each dataset
print(cp_dfs.keys())
[cp_dfs[x].shape for x in cp_dfs]


# ## Generate UMAP coordinates for each plate
# 
# **Note:** Only metadata that is common between plates are included in final data frame.

# In[5]:


desired_columns = [
    "Metadata_Plate",
    "Metadata_Well",
    "Metadata_Site",
    "Metadata_CellLine",
    "Metadata_Condition",
    "Metadata_Nuclei_Site_Count",
    "Metadata_Nuclei_Location_Center_X",
    "Metadata_Nuclei_Location_Center_Y"
]

# Fit UMAP features per dataset and save
for plate in cp_dfs:
    plate_name = pathlib.Path(plate).stem
    print("UMAP embeddings being generated for", plate_name)

    # Make sure to reinitialize UMAP instance per plate
    umap_fit = umap.UMAP(random_state=umap_random_seed, n_components=umap_n_components)

    # Make sure NA columns have been removed
    cp_df = cp_dfs[plate]
    cp_df = feature_select(cp_df, operation="drop_na_columns", na_cutoff=0)

    # Process cp_df to separate features and metadata
    cp_features = infer_cp_features(cp_df)
    meta_features = infer_cp_features(cp_df, metadata=True)
    filtered_meta_features = [
        feature for feature in meta_features if feature in desired_columns
    ]

    # Fit UMAP and convert to pandas DataFrame
    embeddings = pd.DataFrame(
        umap_fit.fit_transform(cp_df.loc[:, cp_features]),
        columns=[f"UMAP{x}" for x in range(0, umap_n_components)],
    )
    print(embeddings.shape)

    # Combine with metadata
    cp_umap_with_metadata_df = pd.concat(
        [cp_df.loc[:, filtered_meta_features].reset_index(drop=True), embeddings],
        axis=1,
    )

    # randomize the rows of the dataframe to plot the order of the data evenly
    cp_umap_with_metadata_df = cp_umap_with_metadata_df.sample(frac=1, random_state=0)

    # Generate output file and save
    output_umap_file = pathlib.Path(output_dir, f"UMAP_{plate_name}.tsv")
    cp_umap_with_metadata_df.to_csv(output_umap_file, index=False, sep="\t")


# In[6]:


# Print an example output file
print(cp_umap_with_metadata_df.shape)
cp_umap_with_metadata_df.head(10)


# In[7]:


# Sort the DataFrame by Metadata_Nuclei_Site_Count in ascending order
sorted_df = cp_umap_with_metadata_df.sort_values(by='Metadata_Nuclei_Site_Count')

# Print the shape of the sorted DataFrame
print(sorted_df.shape)

# Display the first 10 rows of the sorted DataFrame
sorted_df.head(10)


# In[8]:


# Filter rows where UMAP1 column has values above 10
filtered_df = cp_umap_with_metadata_df[cp_umap_with_metadata_df['UMAP1'] < 0]

# Print the shape of the filtered DataFrame
print(filtered_df.shape)

# Display the first 10 rows of the filtered DataFrame
filtered_df


# In[9]:


# Group by Metadata_Well and Metadata_Site and count the number of unique combinations
combo_counts = filtered_df.groupby(['Metadata_Well', 'Metadata_Site']).size().reset_index(name='counts')

# Display the resulting DataFrame with counts
print(combo_counts)


# In[10]:


# Filter rows where UMAP1 column has values above 10
filtered_df = cp_umap_with_metadata_df[cp_umap_with_metadata_df['UMAP0'] > 5]

# Print the shape of the filtered DataFrame
print(filtered_df.shape)

# Display the first 10 rows of the filtered DataFrame
filtered_df


# In[11]:


# Group by Metadata_Well and Metadata_Site and count the number of unique combinations
combo_counts = filtered_df.groupby(['Metadata_Well', 'Metadata_Site']).size().reset_index(name='counts')

# Display the resulting DataFrame with counts
print(combo_counts)


# In[ ]:




