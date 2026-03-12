#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using coSMicQC.
# To filter the single-cells, the default method is z-score to find outliers using the values from only one feature at a time. 
# We use features from the AreaShape module to assess the quality of the segmented single-cells:
# 
# ### Assessing poor nuclei segmentation
# 
# Segmentation parameters, though optimized, aren't perfect and may not segment all nuclei "correctly".
# When we say "correctly", the basic definition we use is that a segmentation is around a nuclei and it is mostly accurate (can be slightly under or over segmented).
# We are looking to identify technically "incorrect" segmentations, including clusters of nuclei, segmented background, multiple segmentation in one nuclei, etc.
# To identify nuclei experiencing mis-segmentation, we use:
# 
# - **Nuclei Area:** This metric quantifies the number of pixels in a nucleus segmentation. We detect nuclei that are abnormally large or small, which likely indicates poor nucleus segmentation.
# - **Nuclei FormFactor:** This metric quantifies how circular an object is. The equation used is `4*π*Area/Perimeter^2`. The range of values are 0 to 1 where 1 means a perfect circle and 0 meaning a very odd shaped nuclei. We are looking to remove nuclei that have low `roundness` as we have found these nuclei have rough edges due to the algorithm struggling to segment.
# - **Nuclei Eccentricity:** This metric quantifies how elongated an ellipse is. The range is 0 to 1 where 1 means it is a line shape and 0 is a circle. We are looking to remove nuclei that are very close to line shape as we have found these are related to mis-segmentations from mKate channel overlap into the Hoechst channel.
# 
# We use these features either in combination or alone when finding technical segmentation errors.
# To be specific, we use:
# 
# | Measurement(s)           | Target                                                           | Notes                                                                                                                                                        |
# |--------------------------|------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
# | Nuclei Area + FormFactor | Small and large area that are very non-circular (small formfactor) | We use two combinations of this where one finds the non-circular small cells and the other non-circular large area.                                           |
# | Nuclei Eccentricity      | Segmentations that look like a straight line (too elongated)      | We use this by itself since the poor segmentations can be any size when it is a mis-segmented straight line.                                                  |
# 

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import cosmicqc as cosm


# ## Set paths and variables

# In[2]:


# Directory with data
data_dir = pathlib.Path("./data/converted_profiles/")

# Directory to save cleaned data
cleaned_dir = pathlib.Path("./data/cleaned_sc_profiles/")
cleaned_dir.mkdir(exist_ok=True)

# Directory to save qc result CSV files
qc_results_dir = pathlib.Path("./qc_results")
qc_results_dir.mkdir(exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path("./qc_figures")
qc_fig_dir.mkdir(exist_ok=True)

# Create an empty dictionary to store data frames for each plate
all_qc_data_frames = {}

# metadata columns to include in output data frame
metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
    "Image_FileName_DAPI"
]

# Get a list of the plates to be used for single-cell QC
plates = [plate.name.split("_")[0] for plate in data_dir.iterdir() if plate.is_dir()]

# Print the list of plates
print("These are plates identified:")
for plate in plates:
    print(plate)


# ## Concat all plates together to assess quality control

# In[3]:


# Create an empty list to store data frames for each plate
all_plate_dfs = []

# Iterate through each plate and create the specified data frame
for plate in plates:
    # Only process the files that are in the plate names list
    plate_path = pathlib.Path(f"{data_dir}/{plate}/per_nuclei.parquet")

    # Read the parquet file into a DataFrame
    qc_df = pd.read_parquet(plate_path)

    # Append the data frame to the list
    all_plate_dfs.append(qc_df)

# Concatenate data frames for each plate
concat_df = pd.concat(all_plate_dfs)

# Reset the index so that the index doesn't reset per plate and cause issues downstream
concat_df.reset_index(drop=True, inplace=True)

print(concat_df.shape)
concat_df.head()


# ## Evaluate the distributions of the features across plates to determine best features to use for QC

# In[4]:


# Set the style
sns.set_style("whitegrid")

# Create subplots with 2 rows and 2 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 6)) # first is width and second is height

# Flatten the axes array for easy indexing
axes = axes.flatten()

# Create kdeplot for Nuclei_AreaShape_Area
sns.kdeplot(
    data=concat_df,
    x="Nuclei_AreaShape_Area",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[0],  # Set the first subplot
)
axes[0].set_title("Area")

# Create kdeplot for FormFactor
sns.kdeplot(
    data=concat_df,
    x="Nuclei_AreaShape_FormFactor",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[1],
    legend=False,
)
axes[1].set_title("Form Factor")

# Create kdeplot for Eccentricity
sns.kdeplot(
    data=concat_df,
    x="Nuclei_AreaShape_Eccentricity",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[2],
    legend=False,
)
axes[2].set_title("Eccentricity")

plt.suptitle(
    "Distribution of each quality control feature across all plates",
    fontsize=16,
)
plt.tight_layout(rect=[0, 0, 1.0, 1.0])  # Adjust subplot layout

# Save figure
plt.savefig(
    pathlib.Path(f"{qc_fig_dir}/all_plates_QC_features_dist_plot.png"), dpi=500
)

plt.show()


# ## Use Nuclei Area and Nuclei FormFactor to determine nuclei that are small and improperly shaped due to mis-segmentation

# In[5]:


# Set a negative threshold to identify both outlier small nuclei and low formfactor representing poor segmentations
outlier_threshold = -1

# find small nuclei and low intensity
feature_thresholds = {
    "Nuclei_AreaShape_Area": outlier_threshold,
    "Nuclei_AreaShape_FormFactor": outlier_threshold,
}

small_low_formfactor_outliers = cosm.find_outliers(
    df=concat_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

small_low_formfactor_outliers.sort_values(by="Nuclei_AreaShape_Area", ascending=True).head()


# ## Use Nuclei Area and Nuclei FormFactor to determine nuclei that are very large and improper segmentations where multiple nuclei are included in one segmentation

# In[6]:


# find large nuclei segmentations (above mean) and low formfactor
feature_thresholds = {"Nuclei_AreaShape_Area": 2, "Nuclei_AreaShape_FormFactor": -1}

# run function to identify outliers given conditions
large_area_formfactor_outliers_df = cosm.find_outliers(
    df=concat_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

# print out data frame
large_area_formfactor_outliers_df.sort_values(by="Nuclei_AreaShape_Area", ascending=False).head()


# ## Visualize the distribution of cells with Area and FormFactor

# ### Count

# In[7]:


# Hexbin plot with Matplotlib to show the normalized log scale counts across Area and FormFactor
plt.hexbin(concat_df['Nuclei_AreaShape_Area'], concat_df['Nuclei_AreaShape_FormFactor'], gridsize=50, cmap='plasma', norm=mcolors.LogNorm())
plt.colorbar(label='Count')
plt.xlabel('Nuclei Area')
plt.ylabel('Nuclei FormFactor')
plt.title('Distribution of single-cell nuclei area and formfactor')
plt.tight_layout()

# Save figure
plt.savefig(
    pathlib.Path(f"{qc_fig_dir}/nuclei_area_formfactor_count_hexbin.png"),
    dpi=500,
)

plt.show()


# ### Outlier status

# In[8]:


# Reset the default value to 'inlier'
concat_df["Outlier_Status"] = 0

# Update the 'Outlier_Status' column based on the outliers DataFrame using index
concat_df.loc[concat_df.index.isin(small_low_formfactor_outliers.index), "Outlier_Status"] = (
    1
)

# Update the 'Outlier_Status' column based on the outliers DataFrame using index
concat_df.loc[concat_df.index.isin(large_area_formfactor_outliers_df.index), "Outlier_Status"] = (
    1
)
# Ensure 'Outlier_Status' is numeric
concat_df['Outlier_Status'] = pd.to_numeric(concat_df['Outlier_Status'], errors='coerce')

# Define a custom colormap
cmap = mcolors.ListedColormap(['#006400', '#990090'])
bounds = [0, 0.5, 1]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Create hexbin plot with custom color mapping
hb = plt.hexbin(concat_df['Nuclei_AreaShape_Area'], concat_df['Nuclei_AreaShape_FormFactor'],
                C=concat_df['Outlier_Status'], gridsize=30, cmap=cmap, norm=norm)

plt.colorbar(hb, label='Outlier Status (0: Passing, 1: Failing)')
plt.xlabel('Nuclei_AreaShape_Area')
plt.ylabel('Nuclei_AreaShape_FormFactor')
plt.title('Outlier status across single-cell nuclei area and formfactor')
plt.tight_layout()

# Save figure
plt.savefig(
    pathlib.Path(f"{qc_fig_dir}/nuclei_area_formfactor_outlier_status_hexbin.png"),
    dpi=500,
)

plt.show()


# ## Use Nuclei Eccentricity to identify nuclei that are very long due to mis-segmentations

# In[9]:


# find very elongated nuclei segmentations (above mean)
feature_thresholds = {
    "Nuclei_AreaShape_Eccentricity": 2,
}

# run function to identify outliers given conditions
eccent_outliers_df = cosm.find_outliers(
    df=concat_df,
    metadata_columns=metadata_columns,
    feature_thresholds=feature_thresholds
)

# print out data frame
eccent_outliers_df.head()


# ### Visualize eccentricity outliers specifically in a histogram prior to merging all outliers together

# In[10]:


# Reset the default value to 'inlier'
concat_df["Outlier_Status"] = "Single-cell passed QC"

# Update the 'Outlier_Status' column based on the outliers DataFrame using index
concat_df.loc[concat_df.index.isin(eccent_outliers_df.index), "Outlier_Status"] = (
    "Single-cell failed QC"
)

# Create histogram
plt.figure(figsize=(10, 6))
sns.histplot(
    data=concat_df,
    x="Nuclei_AreaShape_Eccentricity",
    hue="Outlier_Status",
    multiple="stack",
    bins=50,  # Adjust the number of bins as needed
    palette={"Single-cell passed QC": "#006400", "Single-cell failed QC": "#990090"},
    legend=True,
)

plt.title(f"Distribution of single-cell nuclei eccentricity for all plates")
plt.xlabel("Nuclei Eccentricity")
plt.ylabel("Single-cell count")
plt.tight_layout()

# Show the legend
plt.legend(
    title="QC Status",
    loc="upper right",
    prop={"size": 10},
    labels=["Failed QC", "Passed QC"],
)

# Save figure
plt.savefig(
    pathlib.Path(f"{qc_fig_dir}/nuclei_eccentricity_histogram.png"),
    dpi=500,
)

plt.show()


# ## Combine all combinations of outlier data frames together to make one nuclei outlier data frame

# ### Concat on common columns to make one nuclei outlier data frame

# In[11]:


# Concat all outlier data frames together on common columns to make one nuclei outliers data frame
common_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Site",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
]

nuclei_outliers_df = pd.concat(
    [
        small_low_formfactor_outliers[common_columns],
        large_area_formfactor_outliers_df[common_columns],
        eccent_outliers_df[common_columns],
    ]
)

# Define column names to be added
columns_to_add = [
    "Small_Area_FormFactor_Failed",
    "Large_Area_FormFactor_Failed",
    "Eccentricity_Failed",
]

# Iterate over each column name
for column_name in columns_to_add:
    # Assign 1 to rows where the column exists in the corresponding dataframe, else 0
    nuclei_outliers_df[column_name] = 0
    nuclei_outliers_df.loc[
        nuclei_outliers_df.index.isin(small_low_formfactor_outliers.index)
        & (column_name == "Small_Area_FormFactor_Failed"),
        column_name,
    ] = 1
    nuclei_outliers_df.loc[
        nuclei_outliers_df.index.isin(large_area_formfactor_outliers_df.index)
        & (column_name == "Large_Area_FormFactor_Failed"),
        column_name,
    ] = 1
    nuclei_outliers_df.loc[
        nuclei_outliers_df.index.isin(eccent_outliers_df.index)
        & (column_name == "Eccentricity_Failed"),
        column_name,
    ] = 1

# drop any duplicates based on index
nuclei_outliers_df = nuclei_outliers_df.drop_duplicates(subset=None)

print(nuclei_outliers_df.shape)
nuclei_outliers_df.head()


# ### Plot a histogram of the number of cells detected across conditions

# In[12]:


# Define the palette
palette = {
    "Small_Area_FormFactor_Failed": "#ff6600",
    "Large_Area_FormFactor_Failed": "#ff0000",
    "Eccentricity_Failed": "#990090",
}

# Count the number failing cells in each condition
column_counts = nuclei_outliers_df[columns_to_add].sum()

# Get the colors from the palette corresponding to each condition
colors = [palette.get(condition, "skyblue") for condition in column_counts.index]

# Plot histogram with colors matching the conditions
plt.figure(figsize=(10, 6))
plt.bar(column_counts.index, column_counts.values, color=colors)
plt.yscale('log')
plt.title("Count of failed single-cells for each condition")
plt.xlabel("Condition")
plt.ylabel("Log of the count of failed single-cells")
plt.xticks(rotation=45)
plt.tight_layout()

# Save figure
plt.savefig(
    pathlib.Path(f"{qc_fig_dir}/qc_conditions_sc_count_histogram.png"),
    dpi=500,
)

plt.show()


# ## Output data frame to go over examples of passing or failing single-cells

# In[13]:


concat_df[
    [
        "Outlier_Status",
        "Nuclei_AreaShape_Area",
        "Nuclei_AreaShape_FormFactor",
        "Image_Metadata_Plate",
        "Image_Metadata_Well",
        "Image_Metadata_Site",
        "Nuclei_Location_Center_X",
        "Nuclei_Location_Center_Y",
    ]
].sort_values(by="Nuclei_AreaShape_Area", ascending=False).head()


# ## Remove any single-cell that was identified as failing QC and save cleaned data per plate

# ### Remove single-cells identified as failing QC from the dataset

# In[14]:


# Identify failing QC single-cells based on index values
outlier_indices = nuclei_outliers_df.index

# Remove rows with outlier indices (identified above) from the concat data frame
concat_df_cleaned = concat_df.drop(outlier_indices)

# Remove columns from z-scoring or assigning outliers (not included for downstream analysis)
columns_to_keep = [
    col
    for col in concat_df_cleaned.columns
    if not col.startswith("Z_Score") and col != "Outlier_Status"
]

# Filter columns based on conditions
concat_df_cleaned = concat_df_cleaned[columns_to_keep]

# Get the number of indices removed during cleaning
num_indices_removed = len(concat_df.index) - len(concat_df_cleaned.index)

# Save cleaned data for each plate and show the number of single-cells removed per plate
for plate, plate_df in concat_df_cleaned.groupby("Image_Metadata_Plate"):
    num_removed_per_plate = len(
        concat_df[concat_df["Image_Metadata_Plate"] == plate]
    ) - len(plate_df)
    plate_df.to_parquet(f"{cleaned_dir}/{plate}_sc_cleaned.parquet")
    print(f"Plate {plate}: Number of single-cells dropped: {num_removed_per_plate}")

# Verify the result
print(f"Number of single-cells dropped: {num_indices_removed}")
print(plate_df.shape)
plate_df.head()

