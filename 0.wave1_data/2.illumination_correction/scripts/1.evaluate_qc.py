#!/usr/bin/env python
# coding: utf-8

# # Whole image quality control metric evaluation
# 
# In this notebook, we will use the outputted QC metrics to start working on developing thresholds using z-score to flag and skip images during CellProfiler illumination correction.
# We merge all 4 plates together to find the most optimal thresholds across all of them.
# Since the A647 and GOLD channels tend to be less bright and have less objects within the image, we will not use them for thresholding and only use the DAPI channel in this notebook.
# 
# **Blur metric to detect out of focus images** -> PowerLogLogSlope
# 
# **Saturation metric to detect large smudges** -> PercentMaximal

# ## Import libraries

# In[1]:


import pathlib
import pandas as pd

from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns


# ## Set paths and load in data frame

# In[2]:


# Set the threshold for identifying outliers with z-scoring for all metrics (# of standard deviations away from mean)
threshold_z = 2

# Directory for figures to be outputted
figure_dir = pathlib.Path("./qc_figures")
figure_dir.mkdir(exist_ok=True)

# Directory containing the QC results
qc_results_dir = pathlib.Path("./qc_results")

# List to store DataFrames
dfs = []

# Iterate over each folder in qc_results and read the Image.csv file
for folder in qc_results_dir.iterdir():
    image_csv_path = folder / "Image.csv"
    if image_csv_path.exists():
        df = pd.read_csv(image_csv_path)
        dfs.append(df)

# Concatenate all DataFrames into a single DataFrame
qc_df = pd.concat(dfs, ignore_index=True)

print(qc_df.shape)
qc_df.head()


# ## Create concat dataframe combining blur and saturation metrics from all channels

# In[3]:


# Define the channel of interest
channel = "DAPI"

# Create a DataFrame for the channel with all Metadata columns (excluding Series and Frame)
df = (
    qc_df.filter(like="Metadata_")
    .drop(columns=[
        "Metadata_Series",
        "Metadata_Frame",
        "Metadata_Channel",
        "Metadata_FileLocation",
    ])
    .copy()
)

# Add ImageQuality columns for the channel
df[f"ImageQuality_PowerLogLogSlope"] = qc_df[f"ImageQuality_PowerLogLogSlope_{channel}"]
df[f"ImageQuality_PercentMaximal"] = qc_df[f"ImageQuality_PercentMaximal_{channel}"]

# Add "Channel" column
df["Channel"] = channel

print(df.shape)
df.head()


# ## Visualize blur metric
# 
# We use PowerLogLogSlope as the metric for blur (most recommended). We are looking for the most negative or blurry images as we have found in this dataset that anything very positive is still a good quality image.

# In[4]:


sns.set_style('whitegrid')
sns.kdeplot(data=df, x='ImageQuality_PowerLogLogSlope', hue='Channel', palette=['b'],fill=True, common_norm=False)
plt.title(f'Density of PowerLogLogSlope in DAPI channel for all plates')
plt.xlabel('ImageQuality_PowerLogLogSlope')
plt.ylabel('Density')

plt.tight_layout()
plt.savefig(pathlib.Path(f"{figure_dir}/all_plates_channels_blur_density.png"), dpi=500)
plt.show()


# In[5]:


summary_statistics = df["ImageQuality_PowerLogLogSlope"].describe()
print(summary_statistics)


# In[6]:


# Calculate Z-scores for the column
z_scores = zscore(df['ImageQuality_PowerLogLogSlope'])

# Identify outlier rows based on Z-scores below the mean since we are looking for the blurriest images (more negative)
blur_outliers = df[z_scores < -3]

print(blur_outliers.shape)
print(blur_outliers['Channel'].value_counts())
blur_outliers.sort_values(by='ImageQuality_PowerLogLogSlope', ascending=False)


# In[7]:


# Calculate the mean and standard deviation
mean_value = df["ImageQuality_PowerLogLogSlope"].mean()
std_dev = df["ImageQuality_PowerLogLogSlope"].std()

# Calculate the threshold values
threshold_value_below_mean = mean_value + -3 * std_dev

# Print the calculated threshold values
print("Threshold for outliers below the mean:", threshold_value_below_mean)


# ## Saturation metric
# 
# We use Percent Maxmimal as the saturation metric. We are looking for images with a large, bright smudge which tends to only occur in the DAPI channel.

# In[8]:


summary_statistics = df["ImageQuality_PercentMaximal"].describe()
print(summary_statistics)


# In[9]:


# Create a histogram plot
plt.figure(figsize=(10, 6))
sns.histplot(df['ImageQuality_PercentMaximal'], color='skyblue', alpha=0.7)

# Set labels and title
plt.ylabel('log(count)')
plt.xlabel('Percent Maximal')
plt.title('Distribution of PercentMaximal metric in the DAPI channel for all plates')
plt.yscale('log')  # Set y-axis to logarithmic scale
plt.tight_layout()

plt.savefig(pathlib.Path(f"{figure_dir}/all_plates_percent_maximal.png"), dpi=500)

# Show the plot
plt.show()


# In[10]:


# Calculate Z-scores for the column
z_scores = zscore(df['ImageQuality_PercentMaximal'])

# Identify outlier rows based on Z-scores greater than as to identify whole images with abnormally high saturated pixels
sat_outliers = df[z_scores > 2]

print(sat_outliers.shape)
print(sat_outliers['Channel'].value_counts())
sat_outliers.sort_values(by='ImageQuality_PercentMaximal', ascending=True).head()


# In[11]:


# Calculate the mean and standard deviation
mean_value = df["ImageQuality_PercentMaximal"].mean()
std_dev = df["ImageQuality_PercentMaximal"].std()

# Calculate the threshold values
threshold_value_above_mean = mean_value + 2 * std_dev

# Print the calculated threshold values
print("Threshold for outliers above the mean:", threshold_value_above_mean)

